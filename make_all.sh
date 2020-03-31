#!/bin/sh
#$ -cwd           # Set the working directory for the job to the current directory
#$ -pe smp 1      # Request 1 core
#$ -l h_rt=240:0:0 # Request 24 hour runtime
#$ -l h_vmem=2G   # Request 1GB RAM
#$ -m ea
#$ -M m.j.haughey@qmul.ac.uk
#$ -t 1-50
#$ -N fractalDimensions

bash run_codes.sh $SGE_TASK_ID $SGE_TASK_LAST $1 $2


