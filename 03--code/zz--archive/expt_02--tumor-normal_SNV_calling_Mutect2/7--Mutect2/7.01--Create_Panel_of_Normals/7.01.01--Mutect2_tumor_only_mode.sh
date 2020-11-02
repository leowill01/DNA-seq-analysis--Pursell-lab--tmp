#!/bin/bash
#SBATCH --job-name tumorOnlyModeMutect2.sh 			### Job name
#SBATCH -o Output_Log-tumorOnlyModeMutect2.sh.log	### File to store output
#SBATCH -e Error_Log-tumorOnlyModeMutect2.sh.err	### File to store error messages
#SBATCH --qos normal						### Quality of service queue
#SBATCH -t 12:00:00							### Time limit in [days]-[hh]:[mm]:[ss]
#SBATCH --nodes=1							### Nodes requested for job
#SBATCH --ntasks-per-node=1					### Number of tasks to be launched per node
#SBATCH --cpus-per-task=20					### Number of cores per node requested. Max 20/node
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
# Run individually on a set of normal samples to create files to use with CreateSomaticPanelOfNormals


################################################################################
# ARGUMENTS
################################################################################
REF_GENOME_FASTA=$1
IN_CLEAN_NORMAL_BAM=$2
SAMPLE_NAME=$3
OUT_PREFIX=$4


################################################################################
# EXECUTION CODE
################################################################################
module load gatk

# $SAMPLE_NAME is not actually the tumor sample name--put the normal sample name because it is a given that you are running for the purposes of making a PoN
# Using a germline variant reference allows for upfront filtering of common germline variant alleles, which effectively omits common germline variant *alleles* from the PoN. I would like to use the previously created bootstrapped knownSites.vcf as a germline resource, but I've not confirmed that this is accurate enough to use as the germline resource, so I won't use it here yet.

gatk Mutect2 \
--reference="$REF_GENOME_FASTA" \
--input="$IN_CLEAN_NORMAL_BAM" \
--tumor="$SAMPLE_NAME" \
-O="$OUT_PREFIX".vcf


################################################################################
# END JOB
################################################################################

echo End Job
