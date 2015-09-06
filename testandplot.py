import pickle
import sys
import pylab as plt
from scipy import stats
import numpy as np
import triangle
from matplotlib.backends.backend_pdf import PdfPages


def lookatresults(data):
    plots, thetas = [], []
    P = data[-1][0].T
    for i in xrange(len(P[0])):
        x = [x[i] for x in P]
        theta = r'$\theta_{3:}$ {1:.2f} +{2:.2f}/-{0:.2f}'.format(
            stats.mstats.mode(x)[0][0]-stats.scoreatpercentile(x, 16),
            stats.mstats.mode(x)[0][0],
            stats.scoreatpercentile(x, 84)-stats.mstats.mode(x)[0][0], i+1)
        thetas.append(r'$\theta_{}$'.format(i+1))
        f = plt.figure()
        plt.subplot(111)
        plt.title(theta)
        ker = stats.gaussian_kde(x)
        plt.hist(x, normed=True, alpha=0.2)
        X = np.linspace(0.0, max(x) + .1*max(x), 1000)
        plt.plot(X,ker(X))
        plt.xlabel(r"$\theta_{}$".format(i+1))
        #plt.savefig('theta_{}.png'.format(i))
        plots.append(f)



    plt.figure()
    f = triangle.corner(P, labels=thetas)
    #plt.savefig('trianle.png'.format(i))
    plots.append(f)

    return plots


def main():
    data = pickle.load(file(sys.argv[1]))
    results = lookatresults(data)

    report_plot = PdfPages(sys.argv[1].replace('.pkl', '_testplots.pdf'))

    for f in results:
        report_plot.savefig(f)

    report_plot.close()

if __name__ == "__main__":
    main()