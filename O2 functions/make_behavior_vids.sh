#!/bin/bash

#SBATCH -c 1                    			# Number of cores requested
#SBATCH -t 10                   			# Runtime in minutes
											# Or use HH:MM:SS or D-HH:MM:SS, instead of just number of minutes
#SBATCH -p priority              			# Partition (queue) to submit to
#SBATCH --mem=4G               			    # memory needed (memory PER CORE)
#SBATCH --mail-user=mmarquis89@gmail.com
#SBATCH --mail-type=END         			# Mail when the job ends  


# Create text file with list of trials/sids/tids
module load matlab/2017a
vidDataDir=$1
vidSaveDir=$2
matlab -nodesktop -r "create_trial_list('$vidDataDir')"
sleep 1

# Start a job for each video
file="$vidDataDir/dirList.txt"
while IFS=',' read -r dirName sid tid
do
	sbatch make_vid.sh $dirName $sid $tid $vidSaveDir
done <"$file"
