#!/usr/bin/env bash

# title: Make panel of normals for GATK Mutect2 tumor-normal variant calling
# author: Leo Williams
# date: 2020-09-18-17:18:54--Friday, Sep 18, 2020 05:18 PM

#SBATCH --job-name=makePON    					### Job name
#SBATCH --output=Output_Log-makePON.log		### File to store output
#SBATCH --error=Error_Log-makePON.err			### File to store error messages
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



# ARGUMENTS ########################################

in_dir_with_normal_pon_vcfs=$1
in_ref_fasta=$2
in_interval_list=$3
out_pon_rootname=$4
out_dir=$5
# in_list_of_normal_bams=$
# in_mutect_tumor_only_script=$

echo "Arguments:"
echo "in_dir_with_normal_pon_vcfs: $in_dir_with_normal_pon_vcfs"
echo "in_ref_fasta: $in_ref_fasta"
echo "in_refout_pon_rootname_fasta: $out_pon_rootname"
echo "out_dir: $out_dir"

# SETUP ########################################

# Set Bash options
	# set globstar option
	shopt -s globstar

# Make timestamped dir for results output
dir_results="${out_dir}/results--$(date +\%Y-\%m-\%d-%H%M%S)"
mkdir "$dir_results"

	# make subdir for GenomicsDB
	dir_genomicsDB="${dir_results}/genomicsDB"
	mkdir "$dir_genomicsDB"

# Make a copy of the script run and put it into the results dir
path_script=$(echo "$0")
cp "$path_script" "$dir_results"

# define function to get absolute filepath
absfilepath () {
	echo "$(cd "$(dirname "$1")"; pwd -P)/$(basename "$1")"
}

# MAIN CODE ########################################

## LOAD MODULES
module load gatk/4.1.8.1

## RUN ANALYSIS

# 2. Create a GenomicsDB for the normal Mutect2 calls
	# Make an arguments file of a list of every tumor-only normal VCF filepath preceded by a `-V` to add to the `GenomcicsDB` command
		# create the empty arguments file
		genomicsDB_arg_file="${dir_genomicsDB}/genomicsDB-arguments.txt"
		touch "$genomicsDB_arg_file"

		# loop through all the newly-made normal VCFs called in tumor-only mode and add them to the arguments file preceded with "-V"
		for normal_vcf in "${in_dir_with_normal_pon_vcfs}/"**/*"normalForPON.vcf.gz"; do
			echo -n "-V $(absfilepath "$normal_vcf") " >> "$genomicsDB_arg_file"
		done

	# create GenomicsDB
	genomicsDB_workspace="${dir_genomicsDB}/pon-DB" # define workspace path
	gatk --java-options "-Xmx4g -Xms4g" GenomicsDBImport \
	-R "$in_ref_fasta" \
	--intervals "$in_interval_list" \
	--genomicsdb-workspace-path "$genomicsDB_workspace" \
	--arguments_file "$genomicsDB_arg_file"
	# --tmp-dir

# 3. Create the Panel of Normals by combining the normal calls
	gatk CreateSomaticPanelOfNormals \
	--reference "$in_ref_fasta" \
	--variant gendb://"$genomicsDB_workspace" \
	--output "${dir_results}/${out_pon_rootname}.PON.vcf.gz"

# END ########################################
echo "End"