import pickle
import sys
import matplotlib.pyplot as plt
from scipy import stats
import numpy as np
import triangle
from matplotlib.backends.backend_pdf import PdfPages
import simple_model


def opt_bin(A,B):

    bounds = [A.min(), B.min(), A.max(), B.max()]
    bounds.sort()
    sizes = [np.sqrt(A.size), np.sqrt(B.size)]
    sizes.sort()

    return np.linspace(bounds[0], bounds[3], sizes[1])


def lookatresults(data, name):
    plots, thetas, modes = [], [], []
    P = data[-1][0]
    for i in xrange(len(P)):
        x = P[i]
        theta = r'$\theta_{3:}$ {1:.2f} +{2:.2f}/-{0:.2f}'.format(
            stats.mstats.mode(x)[0][0]-stats.scoreatpercentile(x, 16),
            stats.mstats.mode(x)[0][0],
            stats.scoreatpercentile(x, 84)-stats.mstats.mode(x)[0][0], i+1)

        modes.append(stats.mstats.mode(x)[0][0])

        thetas.append(r'$\theta_{}$'.format(i+1))
        f = plt.figure()
        plt.suptitle(name)
        plt.subplot(111)
        plt.title(theta)
        ker = stats.gaussian_kde(x)
        plt.hist(x, normed=True, alpha=0.2)
        X = np.linspace(0.0, max(x) + .1*max(x), 1000)
        plt.plot(X,ker(X))
        plt.xlabel(r"$\theta_{}$".format(i+1))
        #plt.savefig('theta_{}.png'.format(i))
        plots.append(f)


    f = plt.figure()
    plt.subplot(211)
    plt.plot(data['epsilon'], 'o-')
    plt.title(r'$\epsilon$')
    plt.subplot(212)
    plt.plot(data['n total'], 'o-')
    plt.title('N Trials')
    plots.append(f)


    alphas = np.linspace(0, 1, data.size)

    for j in xrange(len(data[0][0])):
        f = plt.figure()
        for i, D in enumerate(data):
            F = stats.gaussian_kde(D[0][j])
            x = np.linspace(D[0][j].min(), D[0][j].max(), 300)
            plt.plot(x, F(x), alpha=alphas[i])
            plt.xlabel(r"$\theta_{}$".format(j+1))
            if i == data.size - 1:
                plt.plot(x, F(x), c='m', ls='--', lw=2, zorder=1)

        plots.append(f)


    plt.figure()
    f = triangle.corner(P.T, labels=thetas)
    #plt.savefig('trianle.png'.format(i))
    plots.append(f)

    return plots, modes

def plot_modes(obs, modes, stars, model, data):

        plots = []
        synth = model.generate_data(modes)
        synth_stats = model.summary_stats(synth)
        obs_stats = model.summary_stats(obs)




        f = plt.figure()
        plt.suptitle('Obs Cand.:{}; Sim Cand.:{}'.format(obs.size, synth.size))
        plt.rc('legend', fontsize='xx-small', frameon=False)
        plt.subplot(121)
        bins = opt_bin(obs_stats[0],synth_stats[0])
        plt.hist(obs_stats[0], bins=bins, histtype='step', label='Data')
        plt.hist(synth_stats[0], bins=bins, histtype='step', label='Simulation')
        plt.xlabel(r'$\xi$')
        plt.legend()

        plt.subplot(122)
        bins = opt_bin(obs_stats[1],synth_stats[1])
        plt.hist(obs_stats[1], bins=np.arange(bins.min()-0.5, bins.max()+1.5,
                                              1),
                 histtype='step', label='Data', log=True)
        plt.hist(synth_stats[1], bins=np.arange(bins.min()-0.5, bins.max()+1.5,
                                                1),
                 histtype='step', label='Simulation', log=True)
        plt.xlabel(r'$N_p$')
        plt.legend()

        plots.append(f)

        f = plt.figure()
        plt.suptitle('Normalized Planet Counts')
        plt.subplot(121)
        plt.hist(obs_stats[1], bins=np.arange(bins.min()-0.5, bins.max()+1.5,
                                              1),
                 histtype='step', label='Data', normed=True)
        plt.hist(synth_stats[1], bins=np.arange(bins.min()-0.5, bins.max()+1.5,
                                                1),
                 histtype='step', label='Simulation', normed=True)
        plt.xlabel(r'$N_p$')
        plt.legend()


        plt.subplot(122)
        plt.hist(obs_stats[1], bins=np.arange(bins.min()-0.5, bins.max()+1.5,
                                              1),
                 histtype='step', label='Data', log=True, normed=True)
        plt.hist(synth_stats[1], bins=np.arange(bins.min()-0.5, bins.max()+1.5,
                                                1),
                 histtype='step', label='Simulation', log=True, normed=True)
        plt.xlabel(r'$N_p$')
        plt.legend()

        plots.append(f)

        D = np.zeros(100)
        for i in xrange(100):
            S = model.generate_data(modes)
            SS = model.summary_stats(S)
            D[i] = model.distance_function(obs_stats, SS)

        f = plt.figure()
        bins = opt_bin(D, np.array(data[-1]['D accepted']))
        plt.hist(D, bins=bins, histtype='step', label='Modes')
        plt.hist(data[-1]['D accepted'], bins=bins,
                     histtype='step', label='PMC')
        plt.axvline(data[-1]['epsilon'], lw=3, label=r'$\epsilon$')
        plt.legend()
        plots.append(f)


        for i in synth.dtype.names:
            f = plt.figure()
            bins = np.linspace(synth[i].min(), synth[i].max(),
                            np.sqrt(synth[i].size))
            plt.hist(synth[i], bins=bins, histtype='step', label='Sim')
            plt.xlabel(i)
            plt.title('min={:.4f}; max={:.4f}'.format(synth[i].min(),
                                                      synth[i].max()))
            if i in obs.dtype.names:
                plt.hist(obs[i], bins=bins, histtype='step', label='Data')
                plt.legend()
            plots.append(f)


        return plots




def main():
    data = pickle.load(file(sys.argv[1]))
    obs = pickle.load(file("/".join(sys.argv[1].split('/')[:-1] +
                            ['obs_data.pkl'])))

    stars = pickle.load(file('stars.pkl'))
    model = simple_model.MyModel(stars)

    results, modes = lookatresults(data, sys.argv[1])

    results = results + plot_modes(obs, modes, stars, model, data)

    report_plot = PdfPages(sys.argv[1].replace('.pkl', '_testplots.pdf'))

    for f in results:
        report_plot.savefig(f)

    report_plot.close()

if __name__ == "__main__":
    main()