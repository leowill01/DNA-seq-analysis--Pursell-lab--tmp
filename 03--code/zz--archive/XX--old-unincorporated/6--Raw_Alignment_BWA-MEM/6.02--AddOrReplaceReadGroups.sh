#!/bin/bash
#SBATCH --job-name addOrReplaceReadGroups.sh 			### Job name
#SBATCH -o Output_Log-addOrReplaceReadGroups.sh.log		### File to store output
#SBATCH -e Error_Log-addOrReplaceReadGroups.sh.err		### File to store error messages
#SBATCH --qos normal									### Quality of service queue
#SBATCH -t 12:00:00										### Time limit in [days]-[hh]:[mm]:[ss]
#SBATCH --nodes=1										### Nodes requested for job
#SBATCH --ntasks-per-node=1								### Number of tasks to be launched per node
#SBATCH --cpus-per-task=1								### Number of cores per node requested. Max 20/node
#SBATCH --mem=40G										### Memory requested for job. Max 64GB/node
#SBATCH --mail-type ALL
#SBATCH --mail-user lwilli24@tulane.edu

# NOTE: Benchmark: 1 core & 40G memory (whole exome) 

echo Start Job
echo nodes: $SLURM_JOB_NODELIST
echo job id: $SLURM_JOB_ID
echo Number of tasks: $SLURM_NTASKS

###############################################################################
# JOB DESCRIPTION
###############################################################################
# Uses Picard AddOrReplaceReadGroups to fix read group information on BAM files.
# !!!NOTE!!! I only wrote this to add platform ID (RGPL) which was missing for my BAM files. It will need to be modified to include other read group fields if needed.

# RGID: Must be unique. Usually <flowcellID.lane> for Illumina reads
# RGPU: <flowcellID.lane.sampleName>
# RGLB: Library preparation identifier i.e. group of samples you prepped at the same time to send for sequencing? Just put "lib1"
# RGPL: Valid values = {ILLUMINA, SOLID, LS454, HELICOS, PACBIO}


###############################################################################
# ARGUMENTS
################################################################################
IN_BAM=$1
OUT_BAM_HEADER=$2
SAMPLE_NAME=$3
# READ_GROUP_ID=$4
# Same as sample name if not otherwise specified
# READ_GROUP_PLATFORM_UNIT=$5
# Same as sample name if not otherwise specified
READ_GROUP_PLATFORM_NAME=$4
# Valid values = {ILLUMINA, SOLID, LS454, HELICOS, PACBIO}
# READ_GROUP_LIBRARY=$5
# Same as sample name if not otherwise specified


###############################################################################
# EXECUTION CODE
###############################################################################
module load gatk

gatk AddOrReplaceReadGroups \
--INPUT="$IN_BAM" \
--OUTPUT="$OUT_BAM_HEADER".RG.bam \
--RGSM="$SAMPLE_NAME" \
--RGID="$SAMPLE_NAME" \
--RGPU="$SAMPLE_NAME" \
--RGPL="$READ_GROUP_PLATFORM_NAME" \
--RGLB="$SAMPLE_NAME" \
--SORT_ORDER=coordinate \
--CREATE_INDEX=true \
-TMP_DIR=$TMPDIR


###############################################################################
# END JOB
################################################################################
echo End Job
