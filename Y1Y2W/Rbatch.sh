#!/bin/bash
#SBATCH  -J bmmclainMPC_2_02       # job name
#SBATCH  -N 1                     # Number of nodes
#SBATCH  -n 1                    # number of task set to 1 for serial job, for parallel increase to 28, 48 for Sabine, or Carya respectively
#SBATCH  -o JobOutputDump/Misc/job_name.o%j          # output and error file name (%j expands to jobID)
#SBATCH  -t 00:10:00               # run time (hh:mm:ss) - 1.5 hours
#SBATCH  --mem-per-cpu=3750mb 

module add R

Rscript Graph_Y1Y2W.R 940603

