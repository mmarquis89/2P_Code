#!/bin/bash

#SBATCH -c 4                   			    # Number of cores requested
#SBATCH -t 30                   			# Runtime in minutes
#SBATCH -p priority                			# Partition (queue) to submit to
#SBATCH --mem=32G               			# memory needed (memory PER CORE)
#SBATCH --mail-user=mmarquis89@gmail.com
#SBATCH --mail-type=END         			# Mail when the job ends  


vidDir=$1
sid=$2
echo $sid >> test.txt

fileStr="*sid_${sid}_tid*.avi"
outputFileName="sid_${sid}_AllTrials"
echo $fileStr >> test.txt
echo $outputFileName >> test.txt
module load matlab/2017a
matlab -nodesktop -r "concat_vids('$vidDir', '$fileStr', 'OutputFile', '$outputFileName')"