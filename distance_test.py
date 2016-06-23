from simpleabc import simple_abc
import matplotlib
matplotlib.use('Agg')
import simple_model
import numpy as np
import pickle
from scipy import stats
import pylab as plt
import sys
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

N = 100

stars = pickle.load(file('stars.pkl'))
model = simple_model.MyModel(stars)

prior_bounds = [(0, 90.0), (0, 10), (0, 20), (0, 1)]

model.set_prior([stats.uniform(prior_bounds[0][0],prior_bounds[0][1]),
                 stats.uniform(prior_bounds[1][0],prior_bounds[1][1]),
                 stats.uniform(prior_bounds[2][0],prior_bounds[2][1]),
                  stats.uniform(prior_bounds[3][0],prior_bounds[3][1])])

model.set_epsilon(0.05)
theta_0 = (10, 1.0, 10, 0.5)
obs = model.generate_data(theta_0)

model.set_data(obs)

sum_stat = model.summary_stats(obs)

thetas, distances, X, D = [],[],[],[]

for n in xrange(N):

    theta = model.draw_theta()
    synth = model.generate_data(theta)
    synth_sum = model.summary_stats(synth)
    distance = model.distance_function(sum_stat, synth_sum)

    thetas.append(theta)
    distances.append(distance)

for i in xrange(len(model.prior)):
    T = list(theta_0) 
    Xs = np.linspace(prior_bounds[i][0]+0.01,prior_bounds[i][1]-0.01, 100)
    d = []
    for x in Xs:
        T[i] = x
        synth = model.generate_data(T)
        synth_sum = model.summary_stats(synth)
        distance = model.distance_function(sum_stat, synth_sum)
        d.append(distance)
    X.append(Xs)
    D.append(d)


thetas = np.asarray(thetas).T

out = file('{}_distance_test.pkl'.format(sys.argv[1]), 'w')
pickle.dump([thetas, distances, X, D], out)
out.close()

report_plot = PdfPages('{:}_testplots.pdf'.format(sys.argv[1]))


for i in xrange(thetas.shape[0]):
    f = plt.figure()
    plt.plot(thetas[i], distances, 'o')
    plt.subplots_adjust(hspace=.45)
    plt.axvline(theta_0[i], color='gray')
    plt.xlabel(r'$\theta_{}$'.format(i+1))
    plt.ylabel('Distance')
    report_plot.savefig(f)

for i in xrange(thetas.shape[0]):
    f = plt.figure()
    plt.semilogy(thetas[i], distances, 'o')
    plt.subplots_adjust(hspace=.45)
    plt.axvline(theta_0[i], color='gray')
    plt.xlabel(r'$\theta_{}$'.format(i+1))
    plt.ylabel('Distance')
    report_plot.savefig(f)

for i in xrange(thetas.shape[0]):
    f = plt.figure()
    plt.plot(X[i], D[i], 'o')
    plt.subplots_adjust(hspace=.45)
    plt.axvline(theta_0[i], color='gray')
    plt.xlabel(r'$\theta_{}$'.format(i+1))
    plt.ylabel('Distance')
    report_plot.savefig(f)

for i in xrange(thetas.shape[0]):
    f = plt.figure()
    plt.semilogy(X[i], D[i], 'o')
    plt.subplots_adjust(hspace=.45)
    plt.axvline(theta_0[i], color='gray')
    plt.xlabel(r'$\theta_{}$'.format(i+1))
    plt.ylabel('Distance')
    report_plot.savefig(f)

report_plot.close()
