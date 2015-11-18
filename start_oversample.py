import os
import sys

print "*** Simple Planets - Almost as easy as ABC! ***"
name = str(raw_input('Enter run name to oversample: '))
science = str(raw_input('Type "True" if KNOWN: '))
start = int(raw_input('Enter number step to sample from: '))
min_part = int(raw_input('Enter number of particles pbs job:'))
jobs = int(raw_input('Enter number of pbs jobs: '))


print "Oversampling {}".format(name)





for j in xrange(jobs):

  pbs_head = '''
#!/bin/bash
#
#PBS -N {0:}_osamp_{2:}
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
cd $PBS_O_WORKDIR'''.format(name, n_procs, j)

  pbs_name = 'pbs_scripts/{:}oversample_{:}.pbs'.format(name,j)


  pbs_file = file(pbs_name, 'w')
  print >> pbs_file, pbs_head
  print >> pbs_file, """python simpleplanets_oversample.py {:} {:} {:} {:} {:} {:} {:}
                   """.format(name, 1, start-1, min_part, 1, science, j)

  pbs_file.close()
  
  os.system("qsub {:}".format(pbs_name))

