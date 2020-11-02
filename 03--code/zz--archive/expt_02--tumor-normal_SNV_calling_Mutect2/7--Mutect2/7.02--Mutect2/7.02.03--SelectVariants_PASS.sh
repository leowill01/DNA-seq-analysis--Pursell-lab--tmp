#!/bin/bash
#SBATCH --job-name selectVariantsSNP.sh 			### Job name
#SBATCH -o Output_Log-selectVariantsSNP.sh.log	### File to store output
#SBATCH -e Error_Log-selectVariantsSNP.sh.err	### File to store error messages
#SBATCH --qos normal						### Quality of service queue
#SBATCH -t 12:00:00							### Time limit in [days]-[hh]:[mm]:[ss]
#SBATCH --nodes=1							### Nodes requested for job
#SBATCH --ntasks-per-node=1					### Number of tasks to be launched per node
#SBATCH --cpus-per-task=1					### Number of cores per node requested. Max 20/node
#SBATCH --mem=64G							### Memory requested for job. Max 64GB/node
#SBATCH --mail-type ALL
#SBATCH --mail-user lwilli24@tulane.edu

# NOTE: Benchmark: 1 core, 64GB memory

echo Start Job
echo nodes: $SLURM_JOB_NODELIST
echo job id: $SLURM_JOB_ID
echo Number of tasks: $SLURM_NTASKS

################################################################################
# JOB DESCRIPTION
################################################################################
# Uses GATK SelectVariants to select only SNP variants from a raw variant call VCF created with HaplotypeCaller

################################################################################
# ARGUMENTS
################################################################################
REF_GENOME_FASTA=$1
FILTERED_MUTECT_VCF=$2
HC_SNP_VCF_PREFIX=$3
VARIANT_TYPE=$4 # options: {SNP|INDEL}


################################################################################
# EXECUTION CODE
################################################################################
module load gatk

gatk SelectVariants \
-R "$REF_GENOME_FASTA" \
-V "$FILTERED_MUTECT_VCF" \
-O "$HC_SNP_VCF_PREFIX".hc."$VARIANT_TYPE".vcf \
-select-type "$VARIANT_TYPE" \
--exclude-filtered true

################################################################################
# END JOB
################################################################################
echo End Job
