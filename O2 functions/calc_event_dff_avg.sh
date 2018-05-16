#!/bin/bash

#SBATCH -c 1                    			# Number of cores requested
#SBATCH -t 01:00:00                   			# Runtime in minutes
#SBATCH -p short              			    # Partition (queue) to submit to
#SBATCH --mem=64G             			    # memory needed (memory PER CORE
#SBATCH --mail-user=mmarquis89@gmail.com
#SBATCH --mail-type=END         			# Mail when the job ends  

echo calc_event_dff

expDate=$1
sid=$2

imgSaveDir="/home/mjm60/${expDate}/sid_${sid}/ImagingData"
fileStr="EventData*.mat"
sessionDataFile="rigid_sid_${sid}_sessionFile.mat"
echo $sessionDataFile

module load matlab/2017a
matlab -nodesktop -r "calc_event_dff_avg('$imgSaveDir', '$fileStr', '$sessionDataFile')"