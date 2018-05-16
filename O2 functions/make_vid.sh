#!/bin/bash

#SBATCH -c 1                   			    # Number of cores requested
#SBATCH -t 10                   			# Runtime in minutes
#SBATCH -p short                			# Partition (queue) to submit to
#SBATCH --mem=2G               			    # memory needed (memory PER CORE)
#SBATCH --job-name=makeVid

echo make_vid

vidDataDir=$1
vidSaveDir=$2
sid=$3
tid=$4

module load matlab/2017a
matlab -nodesktop -r "make_vid('$vidDataDir', $sid, $tid, 'OutputDir', '$vidSaveDir')"