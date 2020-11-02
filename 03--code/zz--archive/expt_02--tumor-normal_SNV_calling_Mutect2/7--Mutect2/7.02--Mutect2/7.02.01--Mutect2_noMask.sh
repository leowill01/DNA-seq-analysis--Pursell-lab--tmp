#!/bin/bash
#SBATCH --job-name mutect2NoMask.sh 			### Job name
#SBATCH -o Output_Log-mutect2NoMask.sh.log	### File to store output
#SBATCH -e Error_Log-mutect2NoMask.sh.err	### File to store error messages
#SBATCH --qos long						### Quality of service queue
#SBATCH -t 3-00:00:00						### Time limit in [days]-[hh]:[mm]:[ss]
#SBATCH --nodes=1							### Nodes requested for job
#SBATCH --ntasks-per-node=1					### Number of tasks to be launched per node
#SBATCH --cpus-per-task=5					### Number of cores per node requested. Max 20/node
#SBATCH --mem=64G							### Memory requested for job. Max 64GB/node
#SBATCH --mail-type ALL
#SBATCH --mail-user lwilli24@tulane.edu

# NOTE: benchmarked for 5 cores (2 day runtime). 10 offers minimal improvement.

echo Start Job
echo nodes: $SLURM_JOB_NODELIST
echo job id: $SLURM_JOB_ID
echo Number of tasks: $SLURM_NTASKS


################################################################################
# JOB DESCRIPTION
################################################################################
# Call tumor mutations in a tumor-normal pair using GATK4 MuTect2

# Neither a germline SNP resource nor a Panel of Normals are used bc this is a comparison run to analyze a raw (i.e. no preprocessing)-aligned BAM.

# The -bamout option was added for testing the feature. Description: Generate the reassembled alignments file with -bamout. The bamout alignments contain the artificial haplotypes and reassembled alignments for the normal and tumor and enable manual review of calls. The parameter is not required by the tool but is recommended as adding it costs only a small fraction of the total run time.

# When running raw-aligned BAM files through Mutect2, use the "long" QOS queue because it takes a lot longer to run than pre-procssed BAMs.


################################################################################
# ARGUMENTS
################################################################################
REF_GENOME_FASTA=$1
TUMOR_BAM=$2
TUMOR_BAM_READ_GROUP_SAMPLE_NAME=$3
NORMAL_BAM=$4
NORMAL_BAM_READ_GROUP_SAMPLE_NAME=$5
OUT_VCF_BAM_PREFIX=$6


################################################################################
# EXECUTION CODE
################################################################################
module load gatk

gatk Mutect2 \
-R "$REF_GENOME_FASTA" \
-I "$TUMOR_BAM" \
-tumor "$TUMOR_BAM_READ_GROUP_SAMPLE_NAME" \
-I "$NORMAL_BAM" \
-normal "$NORMAL_BAM_READ_GROUP_SAMPLE_NAME" \
-O "$OUT_VCF_BAM_PREFIX".mutect2.noMask.vcf \
-bamout "$OUT_VCF_BAM_PREFIX".mutect2.noMask.reassem.bam \
-TMP_DIR=$TMPDIR



################################################################################
# END JOB
################################################################################

echo End Job
