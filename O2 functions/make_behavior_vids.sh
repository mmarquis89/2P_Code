#!/bin/bash

#SBATCH -c 1                    			# Number of cores requested
#SBATCH -t 00:10:00                   		# Runtime in minutes
#SBATCH -p priority              			# Partition (queue) to submit to
#SBATCH --mem=2G             			    # memory needed (memory PER CORE
#SBATCH --mail-user=mmarquis89@gmail.com
#SBATCH --mail-type=END         			# Mail when the job ends  

echo make_behavior_vids

vidDataDir=$1
vidSaveDir=$2
sid=$3

# Create text file with list of trials/sids/tids
module load matlab/2017a
matlab -nodesktop -r "create_trial_list('$vidDataDir', $sid)"
sleep 5

#------------ Start a job to create each video -------------------------------------
file="$vidDataDir/dirList.txt"
while IFS=',' read -r dirName tid
do
	sbatch make_vid.sh $dirName $vidSaveDir $sid $tid 
done <"$file"
sleep 5

# ------------ Concatenate vids afterwards -----------------------------------------

jid1=$(sbatch --dependency=singleton --job-name=makeVid concat_vids.sh $vidSaveDir $sid)
jid1="${jid1//[!0-9]/}"

# ---- Count the number of frames in each of the individual trial videos -----------
jid2=$(sbatch --dependency=singleton --job-name=makeVid count_vid_frames.sh $vidSaveDir $sid)
jid2="${jid2//[!0-9]/}"

# ------------- Create optic flow vid for anvil annotation -------------------------
sbatch --dependency=afterok:$jid1:$jid2 create_optic_flow_vid.sh $vidSaveDir $sid

