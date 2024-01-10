#!/bin/bash

## script to run fine-mapping for cis-regions

#!#############################################################
#!#### Modify the options in this section as appropriate ######
#!#############################################################

#! Name of the job:
#SBATCH -J cis_fm_pQTLs

#! Which project should be charged:
#SBATCH -A MRC-EPID-SL0-CPU

#! How many whole nodes should be allocated?
# In your case you should leave this at 1
#SBATCH --nodes=1

#SBATCH --ntasks=1

#SBATCH --cpus-per-task 6

#! define how much memory for each node
#SBATCH --mem-per-cpu=5G

##SBATCH --exclusive

#SBATCH --array=1-6%20

#! Specify required run time
#SBATCH --time=72:00:00

#! What types of email messages do you wish to receive?
#SBATCH --mail-type=FAIL

#SBATCH --output=slurm-%x-%j.out

#! Uncomment this to prevent the job from being requeued (e.g. if
#! interrupted by node failure or system downtime):
##SBATCH --no-requeue

#! Do not change:
#SBATCH -p epid

#! ############################################################
#! Modify the settings below to specify the application's environment, location
#! and launch method:

#! Optionally modify the environment seen by the application
#! (note that SLURM reproduces the environment at submission irrespective of ~/.bashrc):
. /etc/profile.d/modules.sh                # Leave this line (enables the module command)
module purge                               # Removes all modules still loaded
module load rhel7/default-peta4            # REQUIRED - loads the basic environment

#! Insert additional module load commands after this line if needed:

module load gcc/5
module load r-3.6.0-gcc-5.4.0-bzuuksv

## assign directories used for the analysis
DIR='path'

cd ${DIR}

## get file name
export FL="input/${1}"

echo ${FL}

## run as array job
echo "Job ID: $SLURM_ARRAY_TASK_ID"
olink="$(awk -v var="$SLURM_ARRAY_TASK_ID" -F '\t' 'NR == var {print $1}' ${FL})"
id="$(awk -v var="$SLURM_ARRAY_TASK_ID" -F '\t' 'NR == var {print $2}' ${FL})"
chr="$(awk -v var="$SLURM_ARRAY_TASK_ID" -F '\t' 'NR == var {print $3}' ${FL})"
poss="$(awk -v var="$SLURM_ARRAY_TASK_ID" -F '\t' 'NR == var {print $4}' ${FL})"
pose="$(awk -v var="$SLURM_ARRAY_TASK_ID" -F '\t' 'NR == var {print $5}' ${FL})"

echo "Phenotype ${olink} : Chromosome ${chr} : Locus start ${poss} : Locus end ${pose}"

#------------------------------#
## -->      run coloc     <-- ##
#------------------------------#

## run simple coloc with original association statistics
scripts/03_fine_mapping.R "${olink}" "${id}" "${chr}" "${poss}" "${pose}"


