import os
import sys

print "*** Simple Planets - Distance Test ***"
name = str(raw_input('Enter run name: '))

print "Running distance test {:}".format(name)

pbs_head = '''
#!/bin/bash
#
#PBS -N {0:}
#PBS -M abc-sim@psu.edu
#PBS -m abe
#PBS -A ebf11_collab
#PBS -l pmem=4gb
#PBS -l nodes=1:ppn=1
#PBS -l walltime=04:00:00
#PBS -o runs/
#PBS -e runs/
#PBS -j oe
#
cd $PBS_O_WORKDIR
python distance_test.py {0:}'''.format(name)

science = 'dist_test.pbs'
sci_file = file(science, 'w')
print >> sci_file, pbs_head
sci_file.close()

os.system("qsub dist_test.pbs")
