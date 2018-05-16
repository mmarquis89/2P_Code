#!/bin/bash

#SBATCH -c 1                    			# Number of cores requested
#SBATCH -t 00:30:00                   		# Runtime in minutes
#SBATCH -p short                			# Partition (queue) to submit to
#SBATCH --mem=32G               			# memory needed
#SBATCH --mail-user=mmarquis89@gmail.com
#SBATCH --mail-type=END         			# Mail when the job ends  

echo initial_analysis_processing

parentDir=$1
sessionDataFile=$2

module load matlab/2017a
myVar=$(matlab -nodesktop -r "initial_analysis_processing('$parentDir', '$sessionDataFile')")