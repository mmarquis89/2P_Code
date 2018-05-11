#!/bin/bash

#SBATCH -c 4                   			    # Number of cores requested
#SBATCH -t 15                   			# Runtime in minutes
#SBATCH -p short                			# Partition (queue) to submit to
#SBATCH --mem=8G               			    # memory needed (memory PER CORE)
#SBATCH --mail-user=mmarquis89@gmail.com
#SBATCH --mail-type=END         			# Mail when the job ends  
#SBATCH --job-name="makeVid"

vidDataDir=$1
vidSaveDir=$2
sid=$3
tid=$4

module load matlab/2017a
matlab -nodesktop -r "make_vid('$vidDataDir', $sid, $tid, 'OutputDir', '$vidSaveDir')"