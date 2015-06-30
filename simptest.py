import matplotlib
matplotlib.use('Agg')
import simple_abc
import simple_model
#from astropy.io import ascii
import numpy as np
import pickle 
import pylab as plt
import time 
from scipy import stats


def main():
    #TODO Choose new random seed after testing
    np.random.seed(917)

    steps = 5
    eps = 0.25
    min_part = 10

    #stars = pickle.load(file('stars.pkl'))
    stars = pickle.load(file('stars_trimmed.pkl'))
    #obs = pickle.load(file('data.pkl'))

    model = simple_model.MyModel(stars)
    model.set_prior([stats.uniform(0.5, 1.0),
                    stats.uniform(0, 1.0)])

    theta = (0.513265306122, 0.1)

    obs = model.generate_data(theta)
    model.set_data(obs)





    n_procs = [1, 2, 3, 4, 5, 6, 7, 8]

    start = time.time()
    OT = simple_abc.pmc_abc(model, obs, epsilon_0=eps, min_particles=min_part, steps=steps,
    						 target_epsilon=eps, parallel=False, plot=True)
    end = time.time()
    print 'Serial took {}s'.format(end - start)
    out_pickle = file('simptest.pkl', 'w')
    pickle.dump(OT, out_pickle)
    out_pickle.close()



if __name__ == '__main__':
	main()

