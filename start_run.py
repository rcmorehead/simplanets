import os
import sys

print "*** Simple Planets - As easy as ABC! ***"
name = str(raw_input('Enter run name: '))
steps = int(raw_input('Enter number of pmc steps: '))
eps = float(raw_input('Enter initial tolerance size: '))
min_part = int(raw_input('Enter sample size: '))
n_procs = int(raw_input('Enter number of cores (max 20)'))

print "Running {:}, steps = {:}, epsilon = {:}, samples = {:}".format(name,
                                                                      steps,
                                                                      eps,
                                                                      min_part)

os.makedirs('RUNS/{:}'.format(name))
os.makedirs('RUNS/{:}/KNOWN'.format(name))
os.makedirs('RUNS/{:}/SCIENCE'.format(name))

pbs_head = '''
#!/bin/bash
#
#PBS -N {0:}
#PBS -M abc-sim@psu.edu
#PBS -m abe
#PBS -A ebf11_collab
#PBS -l pmem={1:}gb
#PBS -l nodes=1:ppn={1:}
#PBS -l walltime=048:00:00
#PBS -o runs/
#PBS -e runs/
#PBS -j oe
#
cd $PBS_O_WORKDIR'''.format(name, n_procs)

science = 'pbs_scripts/{:}_science.pbs'.format(name)
known = 'pbs_scripts/{:}_known.pbs'.format(name)

sci_file = file(science, 'w')
print >> sci_file, pbs_head
print >> sci_file, """python simpleplanets_kepler.py {:} {:} {:} {:} {:} False
                   """.format(name, steps, eps, min_part, n_procs)
print >> sci_file, """python testandplot.py RUNS/{0}/SCIENCE/{0}_{1}samples_{2}.pkl
                   """.format(name, min_part, steps-1)

sci_file.close()

known_file = file(known, 'w')
print >> known_file, pbs_head
print >> known_file, """python simpleplanets_kepler.py {:} {:} {:} {:} {:} True
                   """.format(name, steps, eps, min_part, n_procs)
print >> known_file, """python testandplot.py RUNS/{0}/KNOWN/{0}_{1}samples_{2}.pkl
                   """.format(name, min_part, steps-1)
known_file.close()

os.system('echo "steps = {1:}, epsilon = {2:}, samples = {3:}" > RUNS/{0:}/{0:}_log.txt'.format(name,
                                                                      steps,
                                                                      eps,
                                                                      min_part))
os.system(' git status -u none >> RUNS/{0:}/{0:}_log.txt'.format(name))
os.system(' git rev-parse HEAD >> RUNS/{0:}/{0:}_log.txt'.format(name))

if len(sys.argv) > 1 and sys.argv[1] == 'local':
    os.system("python simpleplanets_kepler.py {:} {:} {:} {:} False".format(name, steps, eps, min_part))
    os.system("python simpleplanets_kepler.py {:} {:} {:} {:} True".format(name, steps, eps, min_part))
    os.system("python testandplot.py RUNS/{0}/SCIENCE/{0}_{1}samples_{2}.pkl".format(name, min_part, steps-1))
    os.system("python testandplot.py RUNS/{0}/KNOWN/{0}_{1}samples_{2}.pkl".format(name, min_part, steps-1))

else:
    os.system("qsub {:}".format(known))
    os.system("qsub {:}".format(science))
