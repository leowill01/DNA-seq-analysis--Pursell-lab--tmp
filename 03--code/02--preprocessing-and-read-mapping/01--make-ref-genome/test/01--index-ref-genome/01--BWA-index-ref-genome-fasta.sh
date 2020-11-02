#!/bin/bash
#SBATCH --job-name bwaRefGenomeIndexforRaw.sh 			### Job name
#SBATCH -o Output_Log-bwaRefGenomeIndexforRaw.sh.log	### File to store output
#SBATCH -e Error_Log-bwaRefGenomeIndexforRaw.sh.err	### File to store error messages
#SBATCH --qos normal						### Quality of service queue
#SBATCH -t 12:00:00							### Time limit in [days]-[hh]:[mm]:[ss]
#SBATCH --nodes=1							### Nodes requested for job
#SBATCH --ntasks-per-node=1					### Number of tasks to be launched per node
#SBATCH --cpus-per-task=10					### Number of cores per node requested. Max 20/node
#SBATCH --mem=5G							### Memory requested for job. Max 64GB/node
#SBATCH --mail-type ALL
#SBATCH --mail-user lwilli24@tulane.edu

# NOTE: Benchmark: 10 CPU cores, 5GB memory (50min)

echo Start Job
echo nodes: $SLURM_JOB_NODELIST
echo job id: $SLURM_JOB_ID
echo Number of tasks: $SLURM_NTASKS

################################################################################
# JOB DESCRIPTION
################################################################################
# Uses BWA to index reference genome FASTA file.
# Usage: sbatch <script> <genome.fa>

################################################################################
# ARGUMENTS
################################################################################
REF_GENOME_FASTA=$1


################################################################################
# EXECUTION CODE
################################################################################
module load bwa

bwa index \
-a bwtsw \
"$REF_GENOME_FASTA"


################################################################################
# END JOB
################################################################################
echo End Job
