#!/bin/bash

#SBATCH -c 6                    			# Number of cores requested
#SBATCH -t 4:00:00                   			# Runtime in minutes
#SBATCH -p short                			# Partition (queue) to submit to
#SBATCH --mem=100G               			# memory needed
#SBATCH --mail-user=mmarquis89@gmail.com
#SBATCH --mail-type=END         			# Mail when the job ends  

parentDir=$1
sid=$2
expDate=$3
imgSaveDir=$4

module load matlab/2017a
matlab -nodesktop -r "preReg_routine_MM('$parentDir', '$sid', '$expDate', 'OutputDir', '$imgSaveDir')"