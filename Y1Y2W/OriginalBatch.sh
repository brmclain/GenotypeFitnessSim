#!/bin/bash
#SBATCH  -J bmmclainMPC_2_01       # job name
#SBATCH  -N 1                     # Number of nodes
#SBATCH  -n 48                    # number of task set to 1 for serial job, for parallel increase to 28, 48 for Sabine, or Carya respectively
#SBATCH  -o job_name.o%j          # output and error file name (%j expands to jobID)
#SBATCH  -t 3:00:00               # run time (hh:mm:ss) - 1.5 hours
#SBATCH  --mem-per-cpu=3750mb 





# Add the necessary programs to your working environment (optional, only use if you need 
module load python/3.9
cd /project/meisel/users/bmmclain/
python MultiProcessingCode_V2_01.py 1000 1000 48 dominant --control_file control_file.txt output3.txt

