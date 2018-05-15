#!/bin/bash

#SBATCH -c 1                   			    # Number of cores requested
#SBATCH -t 02:00:00                			# Runtime in minutes
#SBATCH -p short                			# Partition (queue) to submit to
#SBATCH --mem=16G             			    # memory needed (memory PER CORE

echo create_optic_flow_vid
parentDir=$1
sid=$2

inputVidName="sid_${sid}_AllTrials"
roiDataFile="Behavior_Vid_ROI_Data.mat"
frameCountFile="sid_${sid}_frameCounts.txt"

module load matlab/2017a
matlab -nodesktop -r "create_optic_flow_vid('$parentDir', '$inputVidName', '$roiDataFile', '$frameCountFile')"
