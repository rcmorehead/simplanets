from simpleabc import simple_abc
from astropy.io import ascii
import simple_model
import numpy as np
import pickle
import pylab as plt
from scipy import stats
import time
import triangle
import matplotlib
from scipy.optimize import minimize

plt.style.use('ggplot')
matplotlib.rcParams['legend.numpoints'] = 1
matplotlib.rcParams['lines.linewidth'] = 2.0
matplotlib.rcParams['axes.labelsize'] = 'xx-large'
matplotlib.rcParams['xtick.labelsize'] = 'xx-large'
matplotlib.rcParams['ytick.labelsize'] = 'xx-large'
matplotlib.rcParams['legend.fontsize'] = 'x-large'
matplotlib.rcParams['xtick.major.size'] = 0
matplotlib.rcParams[ 'xtick.major.width'] = 2
matplotlib.rcParams['xtick.minor.size'] = 3
matplotlib.rcParams['xtick.minor.width'] = 1
matplotlib.rcParams['ytick.major.size'] = 7
matplotlib.rcParams['ytick.major.width'] = 2
matplotlib.rcParams['ytick.minor.size'] = 3
matplotlib.rcParams['ytick.minor.width'] = 1
matplotlib.rcParams['figure.subplot.bottom'] = 0.15
matplotlib.rcParams['savefig.dpi'] = 300

def loadabc(filename):
    data = pickle.load(file(filename))
    obs = pickle.load(file("/".join(filename.split('/')[:-1] +
                            ['obs_data.pkl'])))

    stars = pickle.load(file('stars.pkl'))
    model = simple_model.MyModel(stars)

    f = stats.gaussian_kde(data[-1][0])
    int_guess = np.mean(data[-1][0], axis=1)

    modes = minimize(neg, int_guess, args=(f)).x

    return data, obs, stars, model, modes

def neg(x, function=None):
    return -function(x)


def opt_bin(A,B):

    bounds = [A.min(), B.min(), A.max(), B.max()]
    bounds.sort()
    sizes = [np.sqrt(A.size), np.sqrt(B.size)]
    sizes.sort()

    return np.linspace(bounds[0], bounds[3], sizes[1])

def lookatresults(data, modes, theta=None, vert=False):

    if vert == True:
        subplots = [311, 312, 313]
        figsize = (6,18)
    else:
        subplots = [131, 132, 133]
        figsize = (15, 3)

    f = stats.gaussian_kde(data[-1][0])
    int_guess = np.mean(data[-1][0], axis=1)
    modes = minimize(neg, int_guess, args=(f)).x

    thetas = []
    P = data[-1][0]

    for i in xrange(len(P)):
        x = P[i]
        t = r'$\theta_{3:}$ {1:.2f} +{2:.2f}/-{0:.2f}'.format(
            modes[i]-stats.scoreatpercentile(x, 16),
            modes[i],
            stats.scoreatpercentile(x, 84)-modes[i], i+1)

        thetas.append(t)

    if P.shape[1] > 10:
        bins = np.sqrt(P.shape[1])
    else:
        bins=10
    fig = plt.figure(figsize=figsize)
    plt.subplot(subplots[0])
    plt.title(thetas[0])
    ker = stats.gaussian_kde(P[0])
    plt.hist(P[0], bins=bins, normed=True, alpha=0.2)
    x = np.linspace(0.0,90,1000)
    plt.plot(x,ker(x))
    plt.xlabel(r"$\sigma_{mi}$")
    if theta != None:
        plt.axvline(theta[0])

    plt.subplot(subplots[1])
    plt.title(thetas[1])
    ker = stats.gaussian_kde(P[1])
    plt.hist(P[1], bins=bins, normed=True, alpha=0.2)
    x = np.linspace(0.0,1.0,1000)
    plt.plot(x,ker(x))
    plt.xlabel(r"$\sigma_{e}$")
    if theta != None:
        plt.axvline(theta[1])

    plt.subplot(subplots[2])
    plt.title(thetas[2])
    ker = stats.gaussian_kde(P[2])
    plt.hist(P[2], bins=bins, normed=True, alpha=0.2)
    x = np.linspace(0,20,1000)
    plt.plot(x,ker(x))
    plt.xlabel(r"$\lambda$")
    if theta != None:
        plt.axvline(theta[2])

    return fig


def sim_results(obs, modes, stars, model, data):


    synth = model.generate_data(modes)
    synth_stats = model.summary_stats(synth)
    obs_stats = model.summary_stats(obs)


    f = plt.figure(figsize=(15,3))
    plt.suptitle('Obs Cand.:{}; Sim Cand.:{}'.format(obs.size, synth.size))
    plt.rc('legend', fontsize='xx-small', frameon=False)
    plt.subplot(121)
    bins = opt_bin(obs_stats[0],synth_stats[0])
    plt.hist(obs_stats[0], bins=bins, histtype='step', label='Data', lw=2)
    plt.hist(synth_stats[0], bins=bins, histtype='step', label='Simulation', lw=2)
    plt.xlabel(r'$\xi$')
    plt.legend()

    plt.subplot(122)
    bins = opt_bin(obs_stats[1],synth_stats[1])
    plt.hist(obs_stats[1], bins=np.arange(bins.min()-0.5, bins.max()+1.5,
                                          1),
             histtype='step', label='Data', log=True, lw=2)
    plt.hist(synth_stats[1], bins=np.arange(bins.min()-0.5, bins.max()+1.5,
                                            1),
             histtype='step', label='Simulation', log=True, lw=2)
    plt.xlabel(r'$N_p$')
    plt.legend()
