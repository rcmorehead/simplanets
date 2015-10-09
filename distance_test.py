from simpleabc import simple_abc
import simple_model
import numpy as np
import pickle
from scipy import stats
import pylab as plt
import sys

N = 1000

stars = pickle.load(file('stars.pkl'))
model = simple_model.MyModel(stars)

model.set_prior([stats.uniform(0, 90.0),
                 stats.uniform(0, 1),
                 stats.uniform(0, 20)])

theta_0 = (10, 0.1, 10)
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


f = plt.figure()
for i in xrange(thetas.shape[0]):
    plt.subplot(thetas.shape[0], 1, i+1)
    plt.plot(thetas[i], distances, 'o')
    plt.subplots_adjust(hspace=.45)
    plt.axvline(theta_0[i], color='gray')
    plt.xlabel(r'$\theta_{}$'.format(i+1))
    plt.ylabel('Distance')

f.savefig('{}_distance_test.pdf'.format(sys.argv[1]), dpi=300)
