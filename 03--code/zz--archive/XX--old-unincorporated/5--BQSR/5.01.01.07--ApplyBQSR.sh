#!/bin/bash
#SBATCH --job-name applyBQSR.sh 			### Job name
#SBATCH -o Output_Log-applyBQSR.sh.log	### File to store output
#SBATCH -e Error_Log-applyBQSR.sh.err	### File to store error messages
#SBATCH --qos normal						### Quality of service queue
#SBATCH -t 12:00:00							### Time limit in [days]-[hh]:[mm]:[ss]
#SBATCH --nodes=1							### Nodes requested for job
#SBATCH --ntasks-per-node=1					### Number of tasks to be launched per node
#SBATCH --cpus-per-task=1					### Number of cores per node requested. Max 20/node
#SBATCH --mem=64G							### Memory requested for job. Max 64GB/node
#SBATCH --mail-type ALL
#SBATCH --mail-user lwilli24@tulane.edu

# NOTE: Benchmark: 1 core

echo Start Job
echo nodes: $SLURM_JOB_NODELIST
echo job id: $SLURM_JOB_ID
echo Number of tasks: $SLURM_NTASKS

################################################################################
# JOB DESCRIPTION
################################################################################
# Uses GATK ApplyBQSR to take a recalibration table (made from HaplotypeCaller) and recalibrate original clean BAM files to correct for systematic sequencing biases in base quality scores

################################################################################
# ARGUMENTS
################################################################################
REF_GENOME_FASTA=$1
IN_UNRECAL_BAM=$2
IN_RECAL_TABLE=$3
OUT_RECAL_BAM_HEADER=$4


################################################################################
# EXECUTION CODE
################################################################################
module load gatk

gatk ApplyBQSR \
-R="$REF_GENOME_FASTA" \
-I="$IN_UNRECAL_BAM" \
--bqsr-recal-file="$IN_RECAL_TABLE" \
-O="$OUT_RECAL_BAM_HEADER".BQSR_R1_recal.bam \
-TMP_DIR=$TMPDIR

################################################################################
# END JOB
################################################################################
echo End Job
