#!/bin/bash
#SBATCH --job-name=gunzip	                ### Job name
#SBATCH --output=Output_Log-gunzip.log	### File to store output
#SBATCH --error=Error_Log-gunzip.err	    ### File to store error messages
#SBATCH --qos=normal				    ### Quality of service queue
#SBATCH --time=12:00:00					### Time limit in [days]-[hh]:[mm]:[ss]
#SBATCH --nodes=1						### Nodes requested for job
#SBATCH --ntasks-per-node=1				### Number of tasks to be launched per node
#SBATCH --cpus-per-task=1				### Number of cores per node requested. Max 20/node
#SBATCH --mem=64G						### Memory requested for job. Max 64GB/node
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lwilli24@tulane.edu

echo Start Job
echo nodes: $SLURM_JOB_NODELIST
echo job id: $SLURM_JOB_ID
echo Number of tasks: $SLURM_NTASKS

################################################################################
# JOB DESCRIPTION
################################################################################
# batch gunzip script calls this script for each single sample


################################################################################
# ARGUMENTS
################################################################################
GZ_FILE=$1
OUTPUT_DIR=$2

################################################################################
# EXECUTION CODE
################################################################################

gunzip -c $GZ_FILE > "${OUTPUT_DIR}/${GZ_FILE%.gz}"

################################################################################
# END JOB
################################################################################
echo End Job
