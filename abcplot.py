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

# plt.style.use('ggplot')
# matplotlib.rcParams['legend.numpoints'] = 1
# matplotlib.rcParams['lines.linewidth'] = 2.0
# matplotlib.rcParams['axes.labelsize'] = 'xx-large'
# matplotlib.rcParams['xtick.labelsize'] = 'xx-large'
# matplotlib.rcParams['ytick.labelsize'] = 'xx-large'
# matplotlib.rcParams['legend.fontsize'] = 'x-large'
# matplotlib.rcParams['xtick.major.size'] = 0
# matplotlib.rcParams[ 'xtick.major.width'] = 2
# matplotlib.rcParams['xtick.minor.size'] = 3
# matplotlib.rcParams['xtick.minor.width'] = 1
# matplotlib.rcParams['ytick.major.size'] = 7
# matplotlib.rcParams['ytick.major.width'] = 2
# matplotlib.rcParams['ytick.minor.size'] = 3
# matplotlib.rcParams['ytick.minor.width'] = 1
# matplotlib.rcParams['figure.subplot.bottom'] = 0.15
# matplotlib.rcParams['savefig.dpi'] = 300

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

def lookatresults(data, modes, theta=None, vert=False, labels=None):


    P = data[-1][0]
    n = P.shape[0]

    if labels == None:
        labels = [""] * n
    else:
        pass 

    if vert == True:
        subplots = range(n*100+11,n*100+n+11,1)
        figsize = (6, 3*n)
    elif vert == 'four':
        subplots = [221, 222, 223, 224]
        figsize = (10, 10)
    else:
        subplots = range(100+n*10+1,100+n*10+1+n,1)
        figsize = (5*n, 3)

    f = stats.gaussian_kde(data[-1][0])
    int_guess = np.mean(data[-1][0], axis=1)
    modes = minimize(neg, int_guess, args=(f)).x

    thetas = []
    P = data[-1][0]
    labelpad = 20

    for i in xrange(n):
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
    
    for i in xrange(n):
        print subplots[i]
        plt.subplot(int(subplots[i]))
        #plt.title(thetas[0])
        ker = stats.gaussian_kde(P[i])
        h = plt.hist(P[i], bins=bins, normed=True, alpha=1)
        x = np.linspace(h[1][0],h[1][-1],1000)
        plt.plot(x,ker(x))
        plt.xlabel(labels[i], labelpad=labelpad, fontsize=24)
        if theta != None:
            plt.axvline(theta[0])

    for t in thetas: 
        print t

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
