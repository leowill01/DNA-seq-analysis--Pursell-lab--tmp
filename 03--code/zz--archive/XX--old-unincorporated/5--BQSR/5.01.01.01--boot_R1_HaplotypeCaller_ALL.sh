#!/bin/bash
#SBATCH --job-name HaplotypeCaller.sh 				### Job name
#SBATCH -o Output_Log-HaplotypeCaller.sh.log			### File to store output
#SBATCH -e Error_Log-HaplotypeCaller.sh.err			### File to store error messages
#SBATCH --qos normal												### Quality of service queue
#SBATCH -t 24:00:00													### Time limit in [days]-[hh]:[mm]:[ss]
#SBATCH --nodes=1													### Nodes requested for job
#SBATCH --ntasks-per-node=1											### Number of tasks to be launched per node
#SBATCH --cpus-per-task=1											### Number of cores per node requested. Max 20/node
#SBATCH --mem=64G													### Memory requested for job. Max 64GB/node
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
# Use GATK HaplotypeCaller to make an initial raw germline SNP VCF when no dbSNP known variant sites file is available for organism

################################################################################
# ARGUMENTS
################################################################################
CLEAN_NORMAL_BAM_1=$1
CLEAN_NORMAL_BAM_2=$2
CLEAN_NORMAL_BAM_3=$3
CLEAN_NORMAL_BAM_4=$4
CLEAN_NORMAL_BAM_5=$5
CLEAN_NORMAL_BAM_6=$6
OUT_PREFIX=$7
REF_GENOME_FASTA=$8


################################################################################
# EXECUTION CODE
################################################################################
module load gatk

gatk HaplotypeCaller \
--input "$CLEAN_NORMAL_BAM_1" \
--input "$CLEAN_NORMAL_BAM_2" \
--input "$CLEAN_NORMAL_BAM_3" \
--input "$CLEAN_NORMAL_BAM_4" \
--input "$CLEAN_NORMAL_BAM_5" \
--input "$CLEAN_NORMAL_BAM_6" \
--output "$SAMPLE_ID"_bootGermVar_R1.raw.g.vcf \
--reference "$REF_GENOME_FASTA" \
--emit-ref-confidence GVCF \
-TMP_DIR=$TMPDIR

################################################################################
# END JOB
################################################################################
echo End Job
