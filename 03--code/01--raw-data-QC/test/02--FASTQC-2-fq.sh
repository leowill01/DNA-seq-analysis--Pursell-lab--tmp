#!/bin/bash
#SBATCH --job-name FASTQC.sh 			### Job name
#SBATCH -o Output_Log-FASTQC.sh.log	### File to store output
#SBATCH -e Error_Log-FASTQC.sh.err	### File to store error messages
#SBATCH --qos normal						### Quality of service queue
#SBATCH -t 24:00:00							### Time limit in [days]-[hh]:[mm]:[ss]
#SBATCH --nodes=1							### Nodes requested for job
#SBATCH --ntasks-per-node=1					### Number of tasks to be launched per node
#SBATCH --cpus-per-task=1					### Number of cores per node requested. Max 20/node
#SBATCH --mem=64G							### Memory requested for job. Max 64GB/node
#SBATCH --mail-type ALL
#SBATCH --mail-user lwilli24@tulane.edu

echo Start Job
echo nodes: $SLURM_JOB_NODELIST
echo job id: $SLURM_JOB_ID
echo Number of tasks: $SLURM_NTASKS

################################################################################
# JOB DESCRIPTION
################################################################################
# Quality control check with FASTQC - uses 4 input .fastq files for samples with tumor and normal F+R reads each


################################################################################
# ARGUMENTS
################################################################################
IN_FASTQ_1=$1
IN_FASTQ_2=$2
OUT_DIRECTORY=$3


################################################################################
# EXECUTION CODE
################################################################################
module load fastqc

fastqc \
--format fastq \
--threads 1 \
"$IN_FASTQ_1" \
--outdir "$OUT_DIRECTORY"

fastqc \
--format fastq \
--threads 1 \
"$IN_FASTQ_2" \
--outdir "$OUT_DIRECTORY"


################################################################################
# END JOB
################################################################################
echo End Job
