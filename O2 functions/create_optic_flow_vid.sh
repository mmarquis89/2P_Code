#!/bin/bash

#SBATCH -c 6                   			    # Number of cores requested
#SBATCH -t 60                				# Runtime in minutes
#SBATCH -p short                			# Partition (queue) to submit to
#SBATCH --mem=8G               			    # memory needed (memory PER CORE)
#SBATCH --mail-user=mmarquis89@gmail.com
#SBATCH --mail-type=END         			# Mail when the job ends  


parentDir=$1
sid=$2

inputVidName="sid_${sid}_AllTrials"
roiDataFile="Behavior_Vid_ROI_Data.mat"
frameCountFile="sid_${sid}_frameCounts.txt"

module load matlab/2017a
matlab -nodesktop -r "create_optic_flow_vid('$parentDir', '$inputVidName', '$roiDataFile', '$frameCountFile')"
