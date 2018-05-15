#!/bin/bash

#SBATCH -c 1                   			    # Number of cores requested
#SBATCH -t 2:00:00                			# Runtime in minutes
#SBATCH -p priority                			# Partition (queue) to submit to
#SBATCH --mem=16G               			# memory needed (memory PER CORE)

echo concat_vids

vidDir=$1
sid=$2
echo $sid >> test.txt

fileStr="*sid_${sid}_tid*.avi"
outputFileName="sid_${sid}_AllTrials"
echo $fileStr >> test.txt
echo $outputFileName >> test.txt
module load matlab/2017a
matlab -nodesktop -r "concat_vids('$vidDir', '$fileStr', 'OutputFile', '$outputFileName')"