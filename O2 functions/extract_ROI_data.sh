#!/bin/bash

#SBATCH -c 1                    			# Number of cores requested
#SBATCH -t 00:30:00                   		# Runtime in minutes
#SBATCH -p short                			# Partition (queue) to submit to
#SBATCH --mem=64G               			# memory needed
#SBATCH --mail-user=mmarquis89@gmail.com
#SBATCH --mail-type=END         			# Mail when the job ends  

echo extract_ROI_data

expDate=$1
sid=$2

sessionDataFile="rigid_sid_${sid}_sessionFile.mat"
ROIfile="ROI_metadata.mat"
imgSaveDir="/home/mjm60/${expDate}/sid_${sid}/ImagingData"

module load matlab/2017a
matlab -nodesktop -r "extract_ROI_data('$imgSaveDir', '$sessionDataFile', '$ROIfile')"
