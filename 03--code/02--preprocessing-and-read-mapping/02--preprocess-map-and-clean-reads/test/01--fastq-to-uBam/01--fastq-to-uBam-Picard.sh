#!/bin/bash
#SBATCH --job-name fastqToUbam.sh 			### Job name
#SBATCH -o Output_Log-fastqToUbam.sh.log	### File to store output
#SBATCH -e Error_Log-fastqToUbam.sh.err	### File to store error messages
#SBATCH --qos normal						### Quality of service queue
#SBATCH -t 24:00:00							### Time limit in [days]-[hh]:[mm]:[ss]
#SBATCH --nodes=1							### Nodes requested for job
#SBATCH --ntasks-per-node=1					### Number of tasks to be launched per node
#SBATCH --cpus-per-task=1					### Number of cores per node requested. Max 20/node
#SBATCH --mem=64G							### Memory requested for job. Max 64GB/node
#SBATCH --mail-type ALL
#SBATCH --mail-user lwilli24@tulane.edu

echo Start Job
echo nodes: $SLURM_JOB_NODELIST
echo job id: $SLURM_JOB_ID
echo Number of tasks: $SLURM_NTASKS

###############################################################################
# JOB DESCRIPTION
###############################################################################
# Generate an unmapped BAM file (uBAM) from FASTQ read files to combine with raw BAM alignment, and add read group information
#
# Converts a FASTQ file to an unaligned BAM or SAM file. This tool extracts read sequences and base qualities from the input FASTQ file and writes them out to a new file in unaligned BAM (uBAM) format. Read group information can be provided on the command line.
#
# Three versions of FASTQ quality scales are supported: FastqSanger, FastqSolexa and FastqIllumina (see http://maq.sourceforge.net/fastq.shtml for details). Input FASTQ files can be in GZip format (with .gz extension).
#
# RGID: Must be unique. Usually <flowcellID.lane> for Illumina reads
# RGPU: <flowcellID.lane.sampleName>
# RGLB: Library preparation identifier i.e. group of samples you prepped at the same time to send for sequencing? Just put "lib1"
# RGPL: Valid values = {ILLUMINA, SOLID, LS454, HELICOS, PACBIO}


###############################################################################
# ARGUMENTS
###############################################################################
FWD_READ_FASTQ=$1
REV_READ_FASTQ=$2
OUT_BAM_HEADER=$3
READ_GROUP_SAMPLE_NAME=$4

# READ_GROUP_PLATFORM_NAME=$5 # Valid values = {ILLUMINA, SOLID, LS454, HELICOS, PACBIO}

###############################################################################
# JOB
################################################################################
# Use GATK Picard FastqToSam to make a uBAM from FASTQ files
# Only required options are shown
# <SAMPLE_READ_GROUP_NAME> goes on uBAM read group header
module load gatk

gatk FastqToSam \
--FASTQ="$FWD_READ_FASTQ" \
--FASTQ2="$REV_READ_FASTQ" \
--OUTPUT="$OUT_BAM_HEADER".RG.u.bam \
--SAMPLE_NAME="$READ_GROUP_SAMPLE_NAME" \
--READ_GROUP_NAME="$READ_GROUP_SAMPLE_NAME" \
--PLATFORM="illumina" \
-TMP_DIR=$TMPDIR

# --PLATFORM_UNIT="$READ_GROUP_SAMPLE_NAME" \
# --LIBRARY_NAME="$READ_GROUP_SAMPLE_NAME" \

################################################################################
# END JOB
################################################################################


echo End Job
