#!/bin/bash
#SBATCH --job-name refFastaIndex.sh 			### Job name
#SBATCH -o Output_Log-refFastaIndex.sh.log	### File to store output
#SBATCH -e Error_Log-refFastaIndex.sh.err	### File to store error messages
#SBATCH --qos normal						### Quality of service queue
#SBATCH -t 12:00:00							### Time limit in [days]-[hh]:[mm]:[ss]
#SBATCH --nodes=1							### Nodes requested for job
#SBATCH --ntasks-per-node=1					### Number of tasks to be launched per node
#SBATCH --cpus-per-task=1					### Number of cores per node requested. Max 20/node
#SBATCH --mem=1G							### Memory requested for job. Max 64GB/node
#SBATCH --mail-type ALL
#SBATCH --mail-user lwilli24@tulane.edu

# NOTE: Benchmarked: 1 CPU core, 1GB memory (12s runtime)

echo Start Job
echo nodes: $SLURM_JOB_NODELIST
echo job id: $SLURM_JOB_ID
echo Number of tasks: $SLURM_NTASKS

################################################################################
# JOB DESCRIPTION
################################################################################
# Make an index for the reference genome FASTA with SAMtools. This is needed for use by various GATK tools such as HaplotypeCaller


################################################################################
# ARGUMENTS
################################################################################
REF_GENOME_FASTA=$1

################################################################################
# EXECUTION CODE
################################################################################
module load samtools

samtools faidx "$REF_GENOME_FASTA"

################################################################################
# END JOB
################################################################################
echo End Job
