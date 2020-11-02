#!/bin/bash
#SBATCH --job-name selectVariantsIndel.sh 			### Job name
#SBATCH -o Output_Log-selectVariantsIndel.sh.log	### File to store output
#SBATCH -e Error_Log-selectVariantsIndel.sh.err	### File to store error messages
#SBATCH --qos normal						### Quality of service queue
#SBATCH -t 12:00:00							### Time limit in [days]-[hh]:[mm]:[ss]
#SBATCH --nodes=1							### Nodes requested for job
#SBATCH --ntasks-per-node=1					### Number of tasks to be launched per node
#SBATCH --cpus-per-task=1					### Number of cores per node requested. Max 20/node
#SBATCH --mem=64G							### Memory requested for job. Max 64GB/node
#SBATCH --mail-type ALL
#SBATCH --mail-user lwilli24@tulane.edu

# NOTE: Benchmark: 1 core
# NOTE: Run interactively to avoid queue time


echo Start Job
echo nodes: $SLURM_JOB_NODELIST
echo job id: $SLURM_JOB_ID
echo Number of tasks: $SLURM_NTASKS

################################################################################
# JOB DESCRIPTION
################################################################################
# Uses GATK SelectVariants to select only Indel variants from a raw variant call VCF created with HaplotypeCaller

################################################################################
# ARGUMENTS
################################################################################
REF_GENOME_FASTA=$1
RAW_VAR_VCF=$2
RAW_INDEL_VCF_HEADER=$3


################################################################################
# EXECUTION CODE
################################################################################
module load gatk

gatk SelectVariants \
-R="$REF_GENOME_FASTA" \
-V="$RAW_VAR_VCF" \
-select-type="INDEL" \
-O="$RAW_INDEL_VCF_HEADER".indel.vcf

################################################################################
# END JOB
################################################################################
echo End Job
