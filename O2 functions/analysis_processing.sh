#!/bin/bash

#SBATCH -c 1                    			# Number of cores requested
#SBATCH -t 00:02:00                   			# Runtime in minutes
#SBATCH -p priority              			# Partition (queue) to submit to
#SBATCH --mail-user=mmarquis89@gmail.com
#SBATCH --mail-type=END         			# Mail when the job ends  

echo analysis_processing

expDate=$1
sid=$2

imgDataDir="/n/scratch2/mjm60/${expDate}/ImagingData"
vidDataDir="/n/scratch2/mjm60/${expDate}/BehaviorVideo"
imgSaveDir="/home/mjm60/${expDate}/sid_${sid}/ImagingData"
vidSaveDir="/home/mjm60/${expDate}/sid_${sid}/BehaviorVideo"

sessionDataFile="rigid_sid_${sid}_sessionFile.mat"
echo $sessionDataFile

#------------- Create and save an analysis metadata file ---------------
jid1=$(sbatch initial_analysis_processing.sh $imgSaveDir $sessionDataFile)
jid1="${jid1//[!0-9]/}"

#------------------ Save plane data ---------------------------
sbatch  --dependency=afterok:$jid1 save_plane_data.sh $imgSaveDir $sessionDataFile

#------------------ Calculate overall behavior state dF/F ---------------------------
sbatch  --dependency=afterok:$jid1 behavioral_state_dff_calc.sh $imgSaveDir $sessionDataFile
