#!/bin/sh
#$ -cwd           # Set the working directory for the job to the current directory
#$ -pe smp 1      # Request 1 core
#$ -l h_rt=240:0:0 # Request 24 hour runtime
#$ -l h_vmem=1G   # Request 1GB RAM
#$ -m ea
#$ -M m.j.haughey@qmul.ac.uk
#$ -t 1-1
#$ -N collect_allFractalDimensionValues

inPath=$1
outPath=$1'/all_fractalDimensions.dat'

bash get_all_fractalDimensions.sh $inPath 50 >> $outPath
