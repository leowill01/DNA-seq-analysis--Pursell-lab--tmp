#!/bin/bash
#SBATCH --job-name baseRecalibrator.sh 			### Job name
#SBATCH -o Output_Log-baseRecalibrator.sh.log	### File to store output
#SBATCH -e Error_Log-baseRecalibrator.sh.err	### File to store error messages
#SBATCH --qos normal						### Quality of service queue
#SBATCH -t 12:00:00							### Time limit in [days]-[hh]:[mm]:[ss]
#SBATCH --nodes=1							### Nodes requested for job
#SBATCH --ntasks-per-node=1					### Number of tasks to be launched per node
#SBATCH --cpus-per-task=1					### Number of cores per node requested. Max 20/node
#SBATCH --mem=64G							### Memory requested for job. Max 64GB/node
#SBATCH --mail-type ALL
#SBATCH --mail-user lwilli24@tulane.edu

# NOTE benchmarked


echo Start Job
echo nodes: $SLURM_JOB_NODELIST
echo job id: $SLURM_JOB_ID
echo Number of tasks: $SLURM_NTASKS

################################################################################
# JOB DESCRIPTION
################################################################################
# Uses GATK BaseRecalibrator and a known germline variant file to make a recalibration table from a raw BAM file

################################################################################
# ARGUMENTS
################################################################################
REF_GENOME_FASTA=$1
IN_BQSR1_RECAL_BAM=$2
KNOWNSITES_SNP_VCF=$3
KNOWNSITES_INDEL_VCF=$4
OUT_RECAL_TABLE_HEADER=$5


################################################################################
# EXECUTION CODE
################################################################################
module load gatk

gatk BaseRecalibrator \
-R="$REF_GENOME_FASTA" \
-I="$IN_BQSR1_RECAL_BAM" \
--known-sites="$KNOWNSITES_SNP_VCF" \
--known-sites="$KNOWNSITES_INDEL_VCF" \
-O="$OUT_RECAL_TABLE_HEADER".BQSR_R2.recal.table \
-TMP_DIR=$TMPDIR


################################################################################
# END JOB
################################################################################
echo End Job
