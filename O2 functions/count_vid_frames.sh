#!/bin/bash

#SBATCH -c 2                   			    # Number of cores requested
#SBATCH -t 30                   			# Runtime in minutes
#SBATCH -p short                			# Partition (queue) to submit to

vidDir=$1
sid=$2

module load matlab/2017a
myVar=$(matlab -nodesktop -r "count_vid_frames('$vidDir', $sid)")
