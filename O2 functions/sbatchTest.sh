#!/bin/sh

#SBATCH -t 5                  		# Runtime in minutes
#SBATCH -p priority              			# Partition (queue) to submit to
#SBATCH --mail-user=mmarquis89@gmail.com
#SBATCH --mail-type=END         			# Mail when the job ends  

expDate=$1
sid=$2

#expDate="2018_04_14_exp_1"
#sid=0

imgDataDir="/n/scratch2/mjm60/$expDate/ImagingData"
vidDataDir="/n/scratch2/mjm60/$expDate/BehaviorVideo"
imgSaveDir="/home/mjm60/$expDate/ImagingData"
vidSaveDir="/home/mjm60/$expDate/BehaviorVideo"

# Make average fluorescence vids
#sbatch make_fluorescence_vid.sh "$imgDataDir" "$sid" "$imgSaveDir"

# Make anatomy stack
#sbatch create_anatomy_stack.sh "$imgDataDir" "$imgSaveDir"

#-----------------------------------------------------------------
# Imaging Data
#-----------------------------------------------------------------

#jid1=$(sbatch pre_reg_processing.sh $imgDataDir $sid $expDate $imgSaveDir)
sbatch pre_reg_processing.sh "$imgDataDir" "$sid" "$expDate" "$imgSaveDir"

#sbatch --dependency=afterok:$jid1 NoRMCorre_registration.sh $imgSaveDir $sid $expDate

#-----------------------------------------------------------------
# Behavior Vids
#-----------------------------------------------------------------


## Create behavior vids
#jid1=$(sbatch make_behavior_vids.sh $expDate)

# Concatenate vids
#jid2=$(sbatch concatenate_vids.sh)

