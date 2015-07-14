#import sys
#sys.path.append("/Users/Robert/CODE")
from simpleabc import simple_abc
import simple_model
import numpy as np
import pickle
from scipy import stats
import time


np.random.seed(914)

steps = 15
eps = 0.25
min_part = 100

#stars = pickle.load(file('stars.pkl'))
stars = pickle.load(file('stars_trimmed.pkl'))
#obs = pickle.load(file('data.pkl'))

model = simple_model.MyModel(stars)

#theta = (mututal inclination, eccentricity, planet number)

model.set_prior([stats.uniform(0, 90.0),
                stats.uniform(0, 1.0),
                stats.uniform(0, 20)])


#theta = (0.513265306122, 0.1)
theta = (2.0, 0.1, 5)

obs = model.generate_data(theta)
model.set_data(obs)



n_procs = [1, 2, 3, 4, 5, 6, 7, 8]

start = time.time()
OT = simple_abc.pmc_abc(model, obs, epsilon_0=eps, min_particles=min_part, steps=steps,
                        target_epsilon=eps, parallel=False)
end = time.time()
print 'Serial took {}s'.format(end - start)
out_pickle = file('demo.pkl', 'w')
pickle.dump(OT, out_pickle)
out_pickle.close()


