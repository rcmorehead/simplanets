import pickle
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats
import numpy as np
import triangle
from matplotlib.backends.backend_pdf import PdfPages
import simple_model
from scipy.optimize import minimize

def opt_bin(A,B):

    bounds = [A.min(), B.min(), A.max(), B.max()]
    bounds.sort()
    sizes = [np.sqrt(A.size), np.sqrt(B.size)]
    sizes.sort()

    return np.linspace(bounds[0], bounds[3], sizes[1])

def neg(x, function=None):
    return -function(x)


def lookatresults(data, name, modes):
    plots, thetas = [], []
    P = data[-1][0]

    for i in xrange(len(P)):
        x = P[i]
        theta = r'$\theta_{3:}$ {1:.2f} +{2:.2f}/-{0:.2f}'.format(
            modes[i]-stats.scoreatpercentile(x, 16),
            modes[i],
            stats.scoreatpercentile(x, 84)-modes[i], i+1)



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
    #f = triangle.corner(P.T, labels=thetas)
    #plt.savefig('trianle.png'.format(i))
    #plots.append(f)

    return plots

def plot_modes(obs, modes, stars, model, data):

        plots = []
        synth = model.generate_data(modes)
        synth_stats = model.summary_stats(synth)
        obs_stats = model.summary_stats(obs)




        

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
                pass
                #plt.hist(obs[i], bins=bins, histtype='step', label='Data')
                #plt.legend()
            plots.append(f)


        return plots




def main():
    data = pickle.load(file(sys.argv[1]))
    obs = pickle.load(file("/".join(sys.argv[1].split('/')[:-1] +
                            ['obs_data.pkl'])))

    stars = pickle.load(file('stars.pkl'))
    model = simple_model.MyModel(stars)
    model.set_epsilon(data[-1]['epsilon'])

    f = stats.gaussian_kde(data[-1][0])
    int_guess = np.mean(data[-1][0], axis=1)

    modes = minimize(neg, int_guess, args=(f)).x

    results = lookatresults(data, sys.argv[1], modes)

    results = results + plot_modes(obs, modes, stars, model, data)

    report_plot = PdfPages(sys.argv[1].replace('.pkl', '_testplots.pdf'))

    for f in results:
        report_plot.savefig(f)

    report_plot.close()

if __name__ == "__main__":
    main()
