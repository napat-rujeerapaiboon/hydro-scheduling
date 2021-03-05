import os
import time
from datetime import datetime

# list of samples that we're running
samples = range(25)
logfiles = {}

# ensure that queue is empty initially
os.system ("qstat | grep \"predict\" | wc -l > output.txt")

f = open("output.txt", "r")
count = int(f.read())

if count != 0:
    print ("Aborting: job queue is not empty.")

# remove any potential status files
os.system ("rm success_*.txt")

# submit all jobs
for i in samples:
    fout = open("predict_{}.pbs".format(i), "w")
    fout.write("#!/bin/sh\n")
    fout.write("#PBS -l select=1:ncpus=4:mem=16GB\n")
    fout.write("#PBS -l walltime=1:00:00\n\n")
    fout.write("module load gurobi/9.0.1\n")
    fout.write("cd $HOME/energy/Sourcecode/prediction\n")
    fout.write("export PBS_ARRAY_INDEX={}\n".format(i))
    fout.write("./a.out\n")
    fout.close()
    
    os.system ("qsub predict_{}.pbs > output.txt".format(i))
    f = open("output.txt", "r")
    logfiles[i] = f.read().strip()

