#!/bin/bash
#SBATCH --job-name=var-calling	                ### Job name
#SBATCH --output=Output_Log-var-calling.log	### File to store output
#SBATCH --error=Error_Log-var-calling.err	    ### File to store error messages
#SBATCH --qos=normal				    ### Quality of service queue
#SBATCH --time=24:00:00					### Time limit in [days]-[hh]:[mm]:[ss]
#SBATCH --nodes=1						### Nodes requested for job
#SBATCH --ntasks-per-node=1 			### Number of tasks to be launched per node
#SBATCH --cpus-per-task=1				### Number of cores per node requested. Max 20/node
#SBATCH --mem=64G						### Memory requested for job. Max 64GB/node
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lwilli24@tulane.edu

# NOTE: CPU benchmarked

echo Start Job
echo nodes: $SLURM_JOB_NODELIST
echo job id: $SLURM_JOB_ID
echo Number of tasks: $SLURM_NTASKS

# About ########################################

# This script runs VarScan2 Basic Protocol 2: somatic mutation calling in tumor-normal pairs

# Arguments ########################################

VARSCAN2_JARFILE="/lustre/project/zpursell/leo/code-packages-external/varscan-2.4.3/VarScan.v2.4.3.jar"
TUMOR_BAM=$1
NORMAL_BAM=$2
OUT_PREFIX=$3
OUTPUT_DIR=$4
REF_GENOME_FASTA=$5

# Make output dir
	# Make timestamped run results dir inside output sample run dir
	RESULTS_DIR="${OUTPUT_DIR}/results--$(date +\%Y-\%m-\%d-%H%M%S)"
	mkdir $RESULTS_DIR

echo Reference genome FASTA: "$REF_GENOME_FASTA"
echo Normal BAM: "$NORMAL_BAM"
echo Tumor BAM: "$TUMOR_BAM"
echo VarScan2 jarfile: "$VARSCAN2_JARFILE"

# Code ########################################

# 1. Run SAMtools mpileup on the normal and tumor BAM files
# Load SAMtools module on Cypress
module load samtools

samtools \
mpileup -B -q 1 \
-f "$REF_GENOME_FASTA" \
"$NORMAL_BAM" \
"$TUMOR_BAM" \
> "${RESULTS_DIR}/${OUT_PREFIX}.mpileup"
# TODO: Pipe mpileup directly into VarScan2 command

# 2. Run VarScan2 in somatic mode
java -jar "$VARSCAN2_JARFILE" \
somatic "${RESULTS_DIR}/${OUT_PREFIX}.mpileup" \
"${RESULTS_DIR}/${OUT_PREFIX}" \
--mpileup 1 \
--min-coverage 10 \
--min-var-freq 0.005 \
--somatic-p-value 0.10 \
--output-vcf 1
# IDEA: --strand-filter 1 ?
# NOTE: --min-var-freq is changed from default of 0.08

rm "${RESULTS_DIR}/${OUT_PREFIX}.mpileup"

# 3. Run processSomatic subcommand to divide the output into separate files based on somatic status and confidence
java -jar "$VARSCAN2_JARFILE" \
processSomatic "${RESULTS_DIR}/${OUT_PREFIX}.snp.vcf" \
--min-var-freq 0.005

java -jar "$VARSCAN2_JARFILE" \
processSomatic "${RESULTS_DIR}/${OUT_PREFIX}.indel.vcf" \
--min-var-freq 0.005

# 4. Run additional filter to identify and remove somatic mutations that are likely false positives due to alignment problems around indels.
java -jar "$VARSCAN2_JARFILE" \
somaticFilter "${RESULTS_DIR}/${OUT_PREFIX}.snp.Somatic.hc.vcf" \
--indel-file "${RESULTS_DIR}/${OUT_PREFIX}.indel.vcf" \
--output-file "${RESULTS_DIR}/${OUT_PREFIX}.snp.Somatic.hc.filter.vcf" \
--min-var-freq 0.005

# End job ########################################
echo End Job

