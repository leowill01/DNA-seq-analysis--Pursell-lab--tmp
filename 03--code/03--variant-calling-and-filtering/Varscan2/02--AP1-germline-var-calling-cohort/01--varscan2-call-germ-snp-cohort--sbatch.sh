#!/bin/bash
#SBATCH --job-name=vs-germ-var-cohort	                ### Job name
#SBATCH --output=Output_Log-vs-germ-var.log	### File to store output
#SBATCH --error=Error_Log-vs-germ-var.err	    ### File to store error messages
#SBATCH --qos=normal				    ### Quality of service queue
#SBATCH --time=12:00:00					### Time limit in [days]-[hh]:[mm]:[ss]
#SBATCH --nodes=1						### Nodes requested for job
#SBATCH --ntasks-per-node=1				### Number of tasks to be launched per node
#SBATCH --cpus-per-task=1				### Number of cores per node requested. Max 20/node
#SBATCH --mem=64G						### Memory requested for job. Max 64GB/node
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lwilli24@tulane.edu

echo Start Job
echo nodes: $SLURM_JOB_NODELIST
echo job id: $SLURM_JOB_ID
echo Number of tasks: $SLURM_NTASKS

# title: 01--varscan2-call-germ-snp-cohort--sbatch.sh
# author: Leo Williams
# date: 2020-07-07 2:34 PM

# USAGE ########################################

# sbatch 01--varscan2-call-germ-snp-cohort--sbatch.sh $

# ABOUT ########################################

# This script follows Varscan2 Alternate Protocol 1: Germline Variant Calling in a Cohort of Individuals (Koboldt, et al., Curr Protoc Bioinformatics., 2013) to produce a single VCF with consensus germline SNP calls in all samples in a cohort.

# ARGUMENTS ########################################

ORDERED_BAM_LIST=$1
ORDERED_SAMPLE_LIST_TXT=$2
OUT_PREFIX=$3
OUT_DIR=$4
REF_FASTA="/lustre/project/zpursell/leo/projects/project_02--DNA-seq-analysis/02--data/02--reference-data/genome/mm10--UCSC/mm10_UCSC.fa"
VARSCAN2_JARFILE="/lustre/project/zpursell/leo/code-packages-external/varscan-2.4.3/VarScan.v2.4.3.jar"

# SETUP ########################################

# load lustre modules
module load samtools/1.10

# Make timestamped dir for results output
RESULTS_DIR="${OUT_DIR}/results--$(date +\%Y-\%m-\%d-%H%M%S)"
mkdir "$RESULTS_DIR"

# Make a copy of the script run and put it into the results dir
SCRIPT_PATH=$(echo "$0")
cp "$SCRIPT_PATH" "$RESULTS_DIR"

# MAIN CODE ########################################

# Alternate Protocol 1: Germline Variant Calling in a Cohort of Individuals

# The ability of SAMtools to generate multiple-sample pileup (mpileup) output makes it possible to call variants across a set of many samples simultaneously. This approach, called cross-sample variant calling, is advantageous because it identifies all variant positions in a group of samples, and provides genotype calls for all samples at each one. Most users will prefer to obtain output in variant call format (VCF), which is compatible with many downstream analysis tools (Danecek, Auton et al. 2011).

# Notably, the size of the output file increases exponentially with the number of samples included, since each sample must be genotyped at every variant position (and variant positions unique to that sample must be genotyped in all other samples). Nevertheless, a single master VCF file containing all variant calls in a cohort offers a convenient format for data sharing and downstream analysis.

# The command and parameter settings for cross-sample variant calling are the same as for Basic Protocol 1. At each position, variant calling requirements are evaluated on each individual's mpileup data. Any positions at which at least one sample contains a qualifying variant are reported for the entire cohort.

# Necessary Resources

# Software: SAMtools (http://samtools.sourceforge.net) for generating the input mpileup

# Data: Aligned sequencing data in binary alignment map (BAM) format, with one BAM file for each sample in the cohort. A sample list for use in VCF column headers (e.g. samples.txt), should also be provided. This is a text file, one sample per line, and the order of the samples should match the order of the BAM files given to the SAMtools command in step 1.

# Perform cross-sample variant calling
	# 1. Run SAMtools mpileup on all BAM files in a single command
	MPILEUP_FILE="${RESULTS_DIR}/${OUT_PREFIX}.mpileup"

	samtools mpileup \
	--no-BAQ \
	--min-MQ 1 \
	--fasta-ref "$REF_FASTA" \
	$(echo "$ORDERED_BAM_LIST") > "$MPILEUP_FILE"

	# 2. Run Varscan mpileup2snp to call SNVs (meant to say SNPs since germline?)
	OUT_COHORT_GERM_SNP_VCF="${RESULTS_DIR}/${OUT_PREFIX}.cohort.varScan.snp.vcf"

	# TODO: fix offset of sample name columns by 1 - it puts the title of the .txt file as the 1st sample name
	java -jar "$VARSCAN2_JARFILE" mpileup2snp \
	"$MPILEUP_FILE" \
	--vcf-sample-list "$ORDERED_SAMPLE_LIST_TXT" \
	--min-coverage 10 \
	--min-var-freq 0.20 \
	--p-value 0.05 \
	--output-vcf 1 \
	> "$OUT_COHORT_GERM_SNP_VCF"

	# 3. Run Varscan mpileup2indel to call indels
	OUT_COHORT_GERM_INDEL_VCF="${RESULTS_DIR}/${OUT_PREFIX}.cohort.varScan.indel.vcf"

	java -jar "$VARSCAN2_JARFILE" mpileup2indel \
	"$MPILEUP_FILE" \
	--vcf-sample-list "$ORDERED_SAMPLE_LIST_TXT" \
	--min-coverage 10 \
	--min-var-freq 0.10 \
	--p-value 0.10 \
	--output-vcf 1 \
	> "$OUT_COHORT_GERM_INDEL_VCF"

	# Remove mpileup file
	rm "$MPILEUP_FILE"