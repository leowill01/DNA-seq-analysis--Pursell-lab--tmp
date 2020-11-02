#!/bin/bash
#SBATCH --job-name mutect2BqsrPon.sh 			### Job name
#SBATCH -o Output_Log-mutect2BqsrPon.sh.log	### File to store output
#SBATCH -e Error_Log-mutect2BqsrPon.sh.err	### File to store error messages
#SBATCH --qos normal						### Quality of service queue
#SBATCH -t 12:00:00							### Time limit in [days]-[hh]:[mm]:[ss]
#SBATCH --nodes=1							### Nodes requested for job
#SBATCH --ntasks-per-node=1					### Number of tasks to be launched per node
#SBATCH --cpus-per-task=5					### Number of cores per node requested. Max 20/node
#SBATCH --mem=64G							### Memory requested for job. Max 64GB/node
#SBATCH --mail-type ALL
#SBATCH --mail-user lwilli24@tulane.edu

# NOTE: Benchmarked (7h runtime)
# WARNING: Running with 1 core has timed out at 12h one some samples.

echo Start Job
echo nodes: $SLURM_JOB_NODELIST
echo job id: $SLURM_JOB_ID
echo Number of tasks: $SLURM_NTASKS


################################################################################
# JOB DESCRIPTION
################################################################################
# Call tumor mutations in a tumor-normal pair using GATK4 MuTect2

# A `--germline-resource` is not used because it was not available.

# The -bamout option was added for testing the feature. Description: Generate the reassembled alignments file with -bamout. The bamout alignments contain the artificial haplotypes and reassembled alignments for the normal and tumor and enable manual review of calls. The parameter is not required by the tool but is recommended as adding it costs only a small fraction of the total run time.


################################################################################
# ARGUMENTS
################################################################################
REF_GENOME_FASTA=$1
CLEAN_BQSR_TUMOR_BAM=$2
TUMOR_BAM_READ_GROUP_SAMPLE_NAME=$3
CLEAN_BQSR_NORMAL_BAM=$4
NORMAL_BAM_READ_GROUP_SAMPLE_NAME=$5
PANEL_OF_NORMALS_VCF=$6
BOOT_GERM_SNP_VCF=$7
OUT_VCF_BAM_PREFIX=$8


################################################################################
# EXECUTION CODE
################################################################################
module load gatk

gatk Mutect2 \
-R "$REF_GENOME_FASTA" \
-I "$CLEAN_BQSR_TUMOR_BAM" \
-tumor "$TUMOR_BAM_READ_GROUP_SAMPLE_NAME" \
-I "$CLEAN_BQSR_NORMAL_BAM" \
-normal "$NORMAL_BAM_READ_GROUP_SAMPLE_NAME" \
--panel-of-normals "$PANEL_OF_NORMALS_VCF" \
--germline-resource "$BOOT_GERM_SNP_VCF" \
-O "$OUT_VCF_BAM_PREFIX".mutect2.maskPonBootGerm.vcf \
-bamout "$OUT_VCF_BAM_PREFIX".mutect2.maskPonBootGerm.reassem.bam \
-TMP_DIR=$TMPDIR



################################################################################
# END JOB
################################################################################

echo End Job
