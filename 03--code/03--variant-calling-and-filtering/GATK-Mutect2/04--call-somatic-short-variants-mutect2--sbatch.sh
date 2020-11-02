#!/usr/bin/env bash

# title: Somatic variant calling in tumor-normal pairs with Mutect2
# author: Leo Williams
# date: 2020-09-18-17:18:54--Friday, Sep 18, 2020 05:18 PM

#SBATCH --job-name=mutect2TN    					### Job name
#SBATCH --output=Output_Log-mutect2TN.log		### File to store output
#SBATCH --error=Error_Log-mutect2TN.err			### File to store error messages
#SBATCH --qos=normal				    	### Quality of service queue
#SBATCH --time=24:00:00						### Time limit in [days]-[hh]:[mm]:[ss]
#SBATCH --nodes=1							### Nodes requested for job
#SBATCH --ntasks-per-node=1 				### Number of tasks to be launched per node
#SBATCH --cpus-per-task=1					### Number of cores per node requested. Max 20/node
#SBATCH --mem=64G							### Memory requested for job. Max 64GB/node
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lwilli24@tulane.edu

# IS BENCHMARKED: NO

echo Start Job
echo nodes: $SLURM_JOB_NODELIST
echo job id: $SLURM_JOB_ID
echo Number of tasks: $SLURM_NTASKS

# USAGE ########################################

#

# ABOUT ########################################

# This script follows the GATK pipeline for calling somatic short variants in matched tumor-normal pair samples with Mutect2.

# INPUT:
# - tumor BAM
# - normal BAM
# - reference genome FASTA
# - panel of normals


# OUTPUT:
# - raw unfiltered somatic short variant callset VCF
# - reassembled BAM
# - Mutect2 .stats file

# ARGUMENTS ########################################

in_ref_fasta=$1 # required
in_tumor_bam=$2 # required
in_normal_bam=$3 # required
normal_sample_name=$4 # required
out_dir=$5 # required
in_tumor_bam_basename=$(basename "$in_tumor_bam")
in_normal_bam_basename=$(basename "$in_normal_bam")
out_prefix="${in_tumor_bam_basename%%.*}_T--${in_normal_bam_basename%%.*}_N"
# in_panel_of_normals=$6 # optional
# in_germline_resource=$7 # optional

# SETUP ########################################

# Make timestamped dir for results output
dir_results="${out_dir}/results--$(date +\%Y-\%m-\%d-%H%M%S)"
mkdir "$dir_results"

# Make a copy of the script run and put it into the results dir
path_script=$(echo "$0")
cp "$path_script" "$dir_results"

# MAIN CODE ########################################

## LOAD MODULES
module load gatk/4.1.8.1

## RUN ANALYSIS

# define output filepaths
out_raw_vcf_gz="${dir_results}/${out_prefix}.mutect2.raw.vcf.gz"
out_filtered_vcf_gz="${dir_results}/${out_prefix}.mutect2.filtered.vcf.gz"
out_filtered_pass_vcf_gz="${dir_results}/${out_prefix}.mutect2.filtered.pass.vcf.gz"
# out_f1r2_tar_gz="${dir_results}/${out_prefix}.f1r2.tar.gz"
# out_read_orientation_model="${dir_results}/${out_prefix}.read-orientation-model.tar.gz"


# call somatic variants in tumors with matched normal samples
	echo "Running Mutect2 without Panel of Normals..."

	gatk Mutect2 \
	--reference "$in_ref_fasta" \
	--input "$in_tumor_bam" \
	--input "$in_normal_bam" \
	--normal-sample "$normal_sample_name" \
	--output "$out_raw_vcf_gz"

	# FIXME: disabled to see if running was problem with PON arg
	# --germline-resource "$in_germline_resource" \
	# # if PON is supplied, run with, else run without
	# if [[ -z "$in_panel_of_normals" ]] ; then
	# 	echo "Running Mutect2 with Panel of Normals..."

	# 	gatk Mutect2 \
	# 	--reference "$in_ref_fasta" \
	# 	--input "$in_tumor_bam" \
	# 	--input "$in_normal_bam" \
	# 	--normal-sample "$normal_sample_name" \
	# 	--panel-of-normals "$in_panel_of_normals" \
	# 	--output "$out_raw_vcf_gz"
	# 	# --germline-resource "$in_germline_resource" \
	# 	# --f1r2-tar-gz "$out_f1r2_tar_gz" \
	# else
	# 	echo "Running Mutect2 without Panel of Normals..."

	# 	gatk Mutect2 \
	# 	--reference "$in_ref_fasta" \
	# 	--input "$in_tumor_bam" \
	# 	--input "$in_normal_bam" \
	# 	--normal-sample "$normal_sample_name" \
	# 	--output "$out_raw_vcf_gz"
	# 	# --germline-resource "$in_germline_resource" \
	# fi

# # pass raw data to LearnReadOrientationModel
# gatk LearnReadOrientationModel \
# -I "$out_f1r2_tar_gz" \
# -O "$out_read_orientation_model"

# # run GetPileupSummaries to summarize read support for a set number of known variant sites

# filter the raw Mutect2 calls
gatk FilterMutectCalls \
--reference "$in_ref_fasta" \
--variant "$out_raw_vcf_gz" \
--output "$out_filtered_vcf_gz"
# --contamination-table contamination.table \
# --tumor-segmentation segments.tsv \

# make separate VCF that has only the calls with '.' or 'PASS'
bcftools view -f '.,PASS' "$out_filtered_vcf_gz" > "$out_filtered_pass_vcf_gz"

# END ########################################
echo "End"