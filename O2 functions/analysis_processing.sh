#!/bin/bash

#SBATCH -c 1                    			# Number of cores requested
#SBATCH -t 00:02:00                   			# Runtime in minutes
#SBATCH -p priority              			# Partition (queue) to submit to
#SBATCH --mail-user=mmarquis89@gmail.com
#SBATCH --mail-type=END         			# Mail when the job ends  

#-----------------------------------------------------------------
# Initial analysis processing
#-----------------------------------------------------------------