#!/bin/bash
#SBATCH --job-name mpileup.sh 			### Job name
#SBATCH -o Output_Log-mpileup.sh.log	### File to store output
#SBATCH -e Error_Log-mpileup.sh.err	### File to store error messages
#SBATCH --qos long						### Quality of service queue
#SBATCH -t 7-00:00:00							### Time limit in [days]-[hh]:[mm]:[ss]
#SBATCH --nodes=1							### Nodes requested for job
#SBATCH --ntasks-per-node=1					### Number of tasks to be launched per node
#SBATCH --cpus-per-task=1					### Number of cores per node requested. Max 20/node
#SBATCH --mem=64G							### Memory requested for job. Max 64GB/node
#SBATCH --mail-type ALL
#SBATCH --mail-user lwilli24@tulane.edu

# NOTE: benchmarked

echo Start Job
echo nodes: $SLURM_JOB_NODELIST
echo job id: $SLURM_JOB_ID
echo Number of tasks: $SLURM_NTASKS

################################################################################
# JOB DESCRIPTION
################################################################################
# Step 1 of VarScan2 Basic Protocol 2: Somatic Mutation Calling in Tumor-Normal Pairs


################################################################################
# ARGUMENTS
################################################################################
NORMAL_BAM=$1 # Clean BAM
TUMOR_BAM=$2 # Clean BAM
OUT_PREFIX=$3


################################################################################
# EXECUTION CODE
################################################################################
# 1. Run SAMtools mpileup on the normal and tumor BAM files:
module load samtools

samtools mpileup -B -q 1 -f \
"/lustre/project/zpursell/leo/projects/project-02--DNA-seq-analysis-of-mouse-POLE-exo-mutant-tumors/02--raw-data/02--standard-reference/genome/mm10--UCSC/mm10_UCSC.fa" \
"$NORMAL_BAM" \
"$TUMOR_BAM" \
> "$OUT_PREFIX".mpileup


################################################################################
# END JOB
################################################################################
echo End Job
