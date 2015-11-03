from simpleabc import simple_abc
import simple_model
import numpy as np
import pickle
from scipy import stats
import time
import sys

name = sys.argv[1]
steps = int(sys.argv[2])
eps = float(sys.argv[3])
min_part = int(sys.argv[4])
known = sys.argv[5]
#print known, type(known)

if known == "True":
    known = True
else:
    known = False
#print known, type(known)

stars = pickle.load(file('stars.pkl'))

model = simple_model.MyModel(stars)

if known:
    theta_0 = (10, 0.1, 10)
    obs = model.generate_data(theta_0)
    print obs[0:3]
else:
    obs = np.recfromcsv('q1_q17_dr24_koi_mcmc_only.csv',
                        usecols=(1,4,14,32,48), delimiter=",")
    obs = obs[obs['koi_disposition'] != "FALSE POSITIVE"]
    obs.dtype.names = ('ktc_kepler_id','koi_disposition','period',
                       'T', 'koi_prad')
    obs = obs[(obs['period'] >= 10.0) & (obs['period'] <= 320.0) &
                (obs['koi_prad'] <=20.0)]
    print obs[0:3]




model.set_prior([stats.uniform(0, 90.0),
                 stats.uniform(0, 1),
                 stats.uniform(0, 20)])

model.set_data(obs)

start = time.time()
OT = simple_abc.pmc_abc(model, obs, epsilon_0=eps, min_samples=min_part,
                        steps=1, parallel=True, n_procs='all')
if known:
    out_pickle = file('RUNS/{0}/KNOWN/{0}_{1}samples_0.pkl'.format(name,
                                                                min_part), 'w')
    observed = file('RUNS/{0}/KNOWN/obs_data.pkl'.format(name), 'w')
    pickle.dump(obs, observed)
    observed.close()

else:
    out_pickle = file('RUNS/{0}/SCIENCE/{0}_{1}samples_0.pkl'.format(name,
                                                                min_part), 'w')
    observed = file('RUNS/{0}/SCIENCE/obs_data.pkl'.format(name), 'w')
    pickle.dump(obs, observed)
    observed.close()

pickle.dump(OT, out_pickle)
out_pickle.close()

for i in range(1, steps):
    PT = OT
    OT = simple_abc.pmc_abc(model, obs, epsilon_0=eps, min_samples=min_part,
                        resume=PT, steps=1, parallel=True, n_procs='all')
    if known:
        out_pickle = file(
            'RUNS/{0}/KNOWN/{0}_{1}samples_{2}.pkl'.format(name, min_part, i),
            'w')
    else:
        out_pickle = file(
            'RUNS/{0}/SCIENCE/{0}_{1}samples_{2}.pkl'.format(name, min_part, i),
            'w')

    pickle.dump(OT, out_pickle)
    out_pickle.close()

end = time.time()
print 'This run took {}s'.format(end - start)
