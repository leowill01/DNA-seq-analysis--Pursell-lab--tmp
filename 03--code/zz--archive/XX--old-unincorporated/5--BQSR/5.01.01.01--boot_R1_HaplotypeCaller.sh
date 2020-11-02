#!/bin/bash
#SBATCH --job-name HaplotypeCaller.sh 				### Job name
#SBATCH -o Output_Log-HaplotypeCaller.sh.log			### File to store output
#SBATCH -e Error_Log-HaplotypeCaller.sh.err			### File to store error messages
#SBATCH --qos normal												### Quality of service queue
#SBATCH -t 12:00:00													### Time limit in [days]-[hh]:[mm]:[ss]
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
CLEAN_NORMAL_BAM=$1
SAMPLE_ID=$2
REF_GENOME_FASTA=$3


################################################################################
# EXECUTION CODE
################################################################################
module load gatk

gatk HaplotypeCaller \
--input="$CLEAN_NORMAL_BAM" \
--output="$SAMPLE_ID"_bootGermVar_R1.raw.vcf \
--reference="$REF_GENOME_FASTA" \
-TMP_DIR=$TMPDIR

################################################################################
# END JOB
################################################################################
echo End Job
