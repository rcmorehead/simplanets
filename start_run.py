import os

print "*** Simple Planets - As easy as ABC! ***"
name = str(raw_input('Enter run name: '))
steps = int(raw_input('Enter number of pmc steps: '))
eps = float(raw_input('Enter initial tolerance size: '))
min_part = int(raw_input('Enter sample size: '))

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
#PBS -N drow
#PBS -M rcm242@psu.edu
#PBS -m abe
#PBS -A ebf11_collab
#PBS -l pmem=4gb
#PBS -l nodes=1:ppn=1
#PBS -l walltime=048:00:00
#PBS -o runs/
#PBS -e runs/
#PBS -j oe
#
cd $PBS_O_WORKDIR'''

science = 'pbs_scripts/{:}_science.pbs'.format(name)
known = 'pbs_scripts/{:}_known.pbs'.format(name)

sci_file = file(science, 'w')
print >> sci_file, pbs_head
print >> sci_file, """python simpleplanets_kepler.py {:} {:} {:} {:} False
                   """.format(name, steps, eps, min_part)
sci_file.close()

known_file = file(known, 'w')
print >> known_file, pbs_head
print >> known_file, """python simpleplanets_kepler.py {:} {:} {:} {:} True
                   """.format(name, steps, eps, min_part)
known_file.close()

os.system('echo "steps = {1:}, epsilon = {2:}, samples = {3:}" > RUNS/{0:}/{0:}_log.txt'.format(name,
                                                                      steps,
                                                                      eps,
                                                                      min_part))

os.system("qsub {:}".format(known))
os.system("qsub {:}".format(science))