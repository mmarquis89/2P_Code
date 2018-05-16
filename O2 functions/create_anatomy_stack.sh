#!/bin/sh

#SBATCH -c 1                    			# Number of cores requested
#SBATCH -t 2:00:00                   	    # Runtime in minutes
#SBATCH -p short               			    # Partition (queue) to submit to
#SBATCH --mem=16G               			# memory needed (memory PER CORE)
#SBATCH --mail-user=mmarquis89@gmail.com
#SBATCH --mail-type=END         			# Mail when the job ends  


parentDir=$1
imgSaveDir=$2

echo "$parentDir $imgSaveDir" >> test.txt

module load matlab/2017a
matlab -nodesktop -r "create_anatomy_stack('$parentDir', 'OutputDir', '$imgSaveDir')"