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

N = 1000

stars = pickle.load(file('stars.pkl'))
model = simple_model.MyModel(stars)

model.set_prior([stats.uniform(0, 90.0),
                 stats.uniform(0, 1),
                 stats.uniform(0, 20),
                  stats.uniform(0, 1)])

theta_0 = (10, 0.1, 10, 0.5)
obs = model.generate_data(theta_0)

model.set_data(obs)

sum_stat = model.summary_stats(obs)

thetas, distances = [],[]

for n in xrange(N):

    theta = model.draw_theta()
    synth = model.generate_data(theta)
    synth_sum = model.summary_stats(synth)
    distance = model.distance_function(sum_stat, synth_sum)

    thetas.append(theta)
    distances.append(distance)



thetas = np.asarray(thetas).T

out = file('{}_distance_test.pkl'.format(sys.argv[1]), 'w')
pickle.dump([thetas, distances], out)
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

report_plot.close()
