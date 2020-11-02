#!/bin/bash
#SBATCH --job-name=var-filter-only	                ### Job name
#SBATCH --output=Output_Log-var-filter-only.log	### File to store output
#SBATCH --error=Error_Log-var-filter-only.err	    ### File to store error messages
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

echo Start job

# About ########################################

# This script runs only the processSomatic and somaticFilter filtering steps of the VarScan2 Basic Protocol 2: somatic mutation calling in tumor-normal pairs

# Arguments ########################################

REF_GENOME_FASTA="/lustre/project/zpursell/leo/projects/project-02--DNA-seq-analysis-of-mouse-POLE-exo-mutant-tumors/02--raw-data/02--reference/genome/mm10--UCSC/mm10_UCSC.fa"
# NORMAL_BAM=$1
# TUMOR_BAM=$2
OUT_PREFIX=$1
VARSCAN2_JARFILE="/lustre/project/zpursell/leo/code-packages-external/varscan-2.4.3/VarScan.v2.4.3.jar"

echo Reference genome FASTA: $REF_GENOME_FASTA
echo Normal BAM: $NORMAL_BAM
echo Tumor BAM: $TUMOR_BAM
echo VarScan2 jarfile: $VARSCAN2_JARFILE


# Code ########################################
# NOTE: Only runs the filtering steps after variant calling
# mkdir results

# # 1. Run SAMtools mpileup on the normal and tumor BAM files
# # Load SAMtools module on Cypress
# module load samtools

# samtools \
# mpileup -B -q 1 \
# -f $REF_GENOME_FASTA \
# $NORMAL_BAM \
# $TUMOR_BAM \
# > ${OUT_PREFIX}.mpileup
# # TODO: Pipe mpileup directly into VarScan2 command

# # 2. Run VarScan2 in somatic mode
# java -jar $VARSCAN2_JARFILE \
# somatic ${OUT_PREFIX}.mpileup \
# "results/${OUT_PREFIX}" \
# --mpileup 1 \
# --min-coverage 10 \
# --min-var-freq 0.05 \
# --somatic-p-value 0.05 \
# --output-vcf 1
# # IDEA: --strand-filter 1 ?
# # NOTE: --min-var-freq is changed from default of 0.08

# rm ${OUT_PREFIX}.mpileup

# 3. Run processSomatic subcommand to divide the output into separate files based on somatic status and confidence
java -jar $VARSCAN2_JARFILE \
processSomatic "results/${OUT_PREFIX}.snp.vcf" \
--min-var-freq 0.05

java -jar $VARSCAN2_JARFILE \
processSomatic "results/${OUT_PREFIX}.indel.vcf" \
--min-var-freq 0.05

# 4. Run additional filter to identify and remove somatic mutations that are likely false positives due to alignment problems around indels.
java -jar $VARSCAN2_JARFILE \
somaticFilter "results/${OUT_PREFIX}.snp.Somatic.hc.vcf" \
--indel-file "results/${OUT_PREFIX}.indel.vcf" \
--output-file "results/${OUT_PREFIX}.snp.Somatic.hc.filter.vcf" \
--min-var-freq 0.05


# End job ########################################
echo End Job

