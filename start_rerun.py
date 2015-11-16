import os
import sys

print "*** Simple Planets - Almost as easy as ABC! ***"
name = str(raw_input('Enter run name to restart: '))
start = int(raw_input('Enter number of last step: '))
min_part = int(raw_input('Enter number of particles '))
steps = int(raw_input('Enter number of pmc steps: '))
n_procs = int(raw_input('Enter number of cores (max 20)'))

print "Rerunning {:}, steps = {:}".format(name, steps,)


pbs_head = '''
#!/bin/bash
#
#PBS -N {0:}_restart
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

science = 'pbs_scripts/{:}rerun_science.pbs'.format(name)
known = 'pbs_scripts/{:}rerun_known.pbs'.format(name)

sci_file = file(science, 'w')
print >> sci_file, pbs_head
print >> sci_file, """python simpleplanets_rerun.py {:} {:} {:} {:} {:} False
                   """.format(name, steps, start, min_part, n_procs)
print >> sci_file, """python testandplot.py RUNS/{0}/SCIENCE/{0}_{1}samples_{2}.pkl
                   """.format(name, min_part, start+steps)

sci_file.close()

known_file = file(known, 'w')
print >> known_file, pbs_head
print >> known_file, """python simpleplanets_kepler.py {:} {:} {:} {:} {:} True
                   """.format(name, steps, start, min_part, n_procs)
print >> known_file, """python testandplot.py RUNS/{0}/KNOWN/{0}_{1}samples_{2}.pkl
                   """.format(name, min_part, start+steps)
known_file.close()



if len(sys.argv) > 1 and sys.argv[1] == 'local':
    os.system("python simpleplanets_rerun.py {:} {:} {:} {:} {:} False".format(name, steps, start, min_part, n_procs))
    os.system("python simpleplanets_rerun.py {:} {:} {:} {:} {:} True".format(name, steps, start, min_part, n_procs))
    #os.system("python testandplot.py RUNS/{0}/SCIENCE/{0}_{1}samples_{2}.pkl".format(name, min_part, steps-1))
    #os.system("python testandplot.py RUNS/{0}/KNOWN/{0}_{1}samples_{2}.pkl".format(name, min_part, steps-1))

else:
    os.system("qsub {:}".format(known))
    os.system("qsub {:}".format(science))
