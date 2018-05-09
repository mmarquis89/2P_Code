#!/bin/bash

#SBATCH -c 2                    			# Number of cores requested
#SBATCH -t 10                   			# Runtime in minutes
											# Or use HH:MM:SS or D-HH:MM:SS, instead of just number of minutes
#SBATCH -p priority                			# Partition (queue) to submit to
#SBATCH --mem=16G               			# memory needed (memory PER CORE)
#SBATCH --mail-user=mmarquis89@gmail.com
#SBATCH --mail-type=END         			# Mail when the job ends  


module load matlab/2017a
matlab -nodesktop -r


create_anatomy_stack('/n/scratch2/mjm60/2018_04_14_exp_1')