#!/bin/bash
#SBATCH --job-name=mut-sigs	                ### Job name
#SBATCH --output=Output_Log-mut-sigs.log	### File to store output
#SBATCH --error=Error_Log-mut-sigs.err	    ### File to store error messages
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

# Job info ##############################
# This script runs an R script which analyzes a set of VCFs for mutation spectra patterns using the MutationalPatterns R package.

# Arguments ##############################
ANNOTATED_VCFS_DIR=

# Code ##############################
# Load R module
module load R/3.5.2-intel

# Execute Rscript
Rscript 01--MutationalPatterns.R $args?

echo End Job