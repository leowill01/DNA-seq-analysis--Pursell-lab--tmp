#!/bin/bash
#SBATCH --job-name markDuplicates.sh 			### Job name
#SBATCH -o Output_Log-markDuplicates.sh.log	### File to store output
#SBATCH -e Error_Log-markDuplicates.sh.err	### File to store error messages
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
# Mark read duplicates in BAM files that are artifacts from sequencing with MarkDuplicates (Picard).
# It takes a previously-made adapter-removed, sorted, merged, and indexed BAM file as input.

################################################################################
# ARGUMENTS
################################################################################
RG_RMADAP_MERGED_SORTED_BAM=$1
OUT_CLEAN_BAM_HEADER=$2

################################################################################
# EXECUTION CODE
################################################################################
module load gatk

gatk MarkDuplicates \
-I="$RG_RMADAP_MERGED_SORTED_BAM" \
-O="$OUT_CLEAN_BAM_HEADER".clean.bam \
-M="$OUT_CLEAN_BAM_HEADER".clean.bam.metrics.txt \
--CREATE_INDEX=true \
-TMP_DIR=$TMPDIR

################################################################################
# END JOB
################################################################################
echo End Job
