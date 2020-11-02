#!/bin/bash
#SBATCH --job-name=mergeChrFa	                ### Job name
#SBATCH --output=Output_Log-mergeChrFa.log	### File to store output
#SBATCH --error=Error_Log-mergeChrFa.err	    ### File to store error messages
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
# ARGUMENTS
OUT_GENOME_FASTA_HEADER=$1

# Combine individual chromosome FASTA files into a single genome FASTA file
cat mm10_chr_fa_UCSC/*.fa > "$OUT_GENOME_FASTA_HEADER".fa


################################################################################
# END JOB
################################################################################
echo End Job

