#!/bin/bash

#SBATCH -c 4                    			# Number of cores requested
#SBATCH -t 00:50:00                   		# Runtime in minutes
#SBATCH -p short                			# Partition (queue) to submit to
#SBATCH --mem=128G               			# memory needed
#SBATCH --mail-user=mmarquis89@gmail.com
#SBATCH --mail-type=END         			# Mail when the job ends  

echo behavioral_state_dff_calc
parentDir=$1
sessionDataFile=$2

module load matlab/2017a
myVar=$(matlab -nodesktop -r "behavioral_state_dff_calc('$parentDir', '$sessionDataFile')")