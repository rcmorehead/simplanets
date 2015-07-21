from simpleabc import simple_abc
import simple_model
import numpy as np
import pickle
from scipy import stats
import time


steps = 10
eps = 0.25
min_part = 25

stars = pickle.load(file('stars.pkl'))

model = simple_model.MyModel(stars)

obs = np.recfromcsv('04012015_trimmed.csv',usecols=(1,4,14,32),delimiter=",")
obs = obs[obs['koi_disposition'] != "FALSE POSITIVE"]
obs.dtype.names = 'ktc_kepler_id','koi_disposition','period', 'T'

model.set_prior([stats.uniform(0, 90.0),
                 stats.uniform(0, 1),
                 stats.uniform(0, 20)])

model.set_data(obs)

start = time.time()
OT = simple_abc.pmc_abc(model, obs, epsilon_0=eps, min_particles=min_part,
                        steps=steps, target_epsilon=eps, parallel=False)
end = time.time()
print 'Serial took {}s'.format(end - start)
out_pickle = file('kepler_pmc_xi_based.pkl', 'w')
pickle.dump(OT, out_pickle)
out_pickle.close()