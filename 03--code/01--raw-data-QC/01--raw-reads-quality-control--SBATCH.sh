#!/bin/bash
#SBATCH --job-name FQ_QC.sh 			### Job name
#SBATCH -o Output_Log-FQ_QC.sh.log	### File to store output
#SBATCH -e Error_Log-FQ_QC.sh.err	### File to store error messages
#SBATCH --qos normal						### Quality of service queue
#SBATCH -t 24:00:00							### Time limit in [days]-[hh]:[mm]:[ss]
#SBATCH --nodes=1							### Nodes requested for job
#SBATCH --ntasks-per-node=1					### Number of tasks to be launched per node
#SBATCH --cpus-per-task=4					### Number of cores per node requested. Max 20/node
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
# Analyzes quality of raw FASTQ reads with FASTQC


################################################################################
# ARGUMENTS
################################################################################
IN_DIR=$1 # dir of .fq files
OUT_DIR=$2 # dir where FASTQC results will go


################################################################################
# EXECUTION CODE
################################################################################
# FIXME: how to not overlap error and log files so that I can automate the making of the sample's experiment dir?

# Load FastQC module on Cypress
module load fastqc

# Run FastQC on all the .fq files in the sample FASTQ folder
fastqc "${IN_DIR}"/**/*.fq \
--format fastq \
--threads 4 \
--extract \
--outdir "$OUT_DIR"

################################################################################
# END JOB
################################################################################
echo End Job
