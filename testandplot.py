import pickle
import sys
import matplotlib.pyplot as plt
from scipy import stats
import numpy as np
import triangle
from matplotlib.backends.backend_pdf import PdfPages


def lookatresults(data, name):
    plots, thetas = [], []
    P = data[-1][0]
    for i in xrange(len(P)):
        x = P[i]
        theta = r'$\theta_{3:}$ {1:.2f} +{2:.2f}/-{0:.2f}'.format(
            stats.mstats.mode(x)[0][0]-stats.scoreatpercentile(x, 16),
            stats.mstats.mode(x)[0][0],
            stats.scoreatpercentile(x, 84)-stats.mstats.mode(x)[0][0], i+1)
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

    return plots

def plot_modes(data):
    pass


def main():
    data = pickle.load(file(sys.argv[1]))
    obs = pickle.load(file("/".join(sys.argv[1].split('/')[:-1] +
                            ['obs_data.pkl'])))
    results = lookatresults(data, sys.argv[1])


    report_plot = PdfPages(sys.argv[1].replace('.pkl', '_testplots.pdf'))

    for f in results:
        report_plot.savefig(f)

    report_plot.close()

if __name__ == "__main__":
    main()