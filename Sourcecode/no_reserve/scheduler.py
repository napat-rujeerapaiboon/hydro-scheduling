import os
import time
from datetime import datetime

# list of samples that we're running
samples = range(25)
logfiles = {}

# ensure that queue is empty initially
os.system ("qstat | grep \"noreserve\" | wc -l > output.txt")

f = open("output.txt", "r")
count = int(f.read())

if count != 0:
    print ("Aborting: job queue is not empty.")

# remove any potential status files
os.system ("rm success_*.txt")

# submit all jobs
for i in samples:
    fout = open("noreserve_{}.pbs".format(i), "w")
    fout.write("#!/bin/sh\n")
    fout.write("#PBS -l select=1:ncpus=4:mem=16GB\n")
    fout.write("#PBS -l walltime=1:00:00\n\n")
    fout.write("module load gurobi/9.0.1\n")
    fout.write("cd $HOME/energy/Sourcecode/no_reserve\n")
    fout.write("export PBS_ARRAY_INDEX={}\n".format(i))
    fout.write("./a.out 10\n")
    fout.close()
    
    os.system ("qsub noreserve_{}.pbs > output.txt".format(i))
    f = open("output.txt", "r")
    logfiles[i] = f.read().strip()

time.sleep (60)

# periodic scheduler
while True:
    print ("{} -- Periodic scheduler check.".format(datetime.now()))

    # go through each sample
    flog = open("logfile.txt", "a+")
    for i in samples:
        os.system ("qstat | grep \"noreserve_{}.pbs\" | wc -l > output.txt".format(i))

        f = open("output.txt", "r")
        count = int(f.read())

        if count == 0:
            print ("  - Sample {} ({}) completed execution.".format(i, logfiles[i]))
            flog.write("{} -- Sample {} ({}) completed execution.\n".format(datetime.now(), i, logfiles[i]))

            os.system ("qsub noreserve_{}.pbs > output.txt".format(i))
            f = open("output.txt", "r")
            logfiles[i] = f.read().strip()
    flog.close()

    # wait
    time.sleep (60)
