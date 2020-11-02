#!/bin/bash
#SBATCH --job-name=CNVs-varscan2	                ### Job name
#SBATCH --output=Output_Log-CNVs-varscan2.log	### File to store output
#SBATCH --error=Error_Log-CNVs-varscan2.err	    ### File to store error messages
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

################################################################################
# ABOUT
################################################################################

# This script is adapted from the "Recommended Workflow" under the current VarScan2 webpage for copy number calling [http://dkoboldt.github.io/varscan/copy-number-calling.html]
# Other information can be found in "Alternate Protocol 2: Somatic Copy Number Alteration Detection in Tumor-Normal Pairs" from [Koboldt, et. al, 2013, "Using VarScan 2 for Germline Variant Calling and Somatic Mutation Detection"].

# This script 1) creates an mpileup file with Samtools from a tumor and matched normal BAM files, 2) report all contiguous regions that meet the coverage requirement, and 3) adjusts copy number change for GC content and extracts potential homozygous deletions.

# The output is a list of regions with coordinates and log2 ratio of copy number changes.

# A subsequent script uses the R package DNAcopy to perform circular binary segementation (CBS) and plots the copy number data.

################################################################################
# USAGE
################################################################################

# sbatch script.sh

# Run from within /experiment_00/sample_name/ folder using the following command:

################################################################################
# SETUP
################################################################################

# Arguments
normal_bam="/lustre/project/zpursell/leo/projects/project_02--DNA-seq/rsync-from-local/04--analysis/expt_02--read-mapping--postprocessing/mm10-reference-genome/02--read-mapping/m079/normal/results/m079-normal-tail-WES.clean.bam"
tumor_bam="/lustre/project/zpursell/leo/projects/project_02--DNA-seq/rsync-from-local/04--analysis/expt_02--read-mapping--postprocessing/mm10-reference-genome/02--read-mapping/m079/tumor/results/m079-tumor-spleen-WES.clean.bam"
reference_genome_fasta="/lustre/project/zpursell/leo/projects/project_02--DNA-seq/rsync-from-local/02--data/02--reference-data/genome/mm10--UCSC/mm10_UCSC.fa"
varscan_jar="/lustre/project/zpursell/leo/code-packages-external/varscan-2.4.3/VarScan.v2.4.3.jar"
out_basename="m079-TN-CNV"

# Make results directory
results_dir_name="results--VarScan2-CNVs--$(date '+%Y-%m-%d-%H%M%S')"
mkdir $results_dir_name

################################################################################
# ANALYSIS
################################################################################

# TODO: auto calculate $data_ratio?

# 1) Make a tumor-normal .mpileup file and pipe to VarScan2 copynumber command
#     NOTE: More recent webpage doesn't mention calculation of data ratio
#     NOTE: (Koboldt, 2013) says to calculate the "data ratio" beforehand, which is the ratio of uniquely mapped normal reads to uniquely mapped tumor reads, or (unique_mapped_normal_reads)/(unique_mapped_tumor_reads)
    # module load samtools
    # samtools mpileup -q 1 -f $reference_genome_fasta \
    # $normal_bam \
    # $tumor_bam > "$out_basename".mpileup

    java -jar $varscan_jar copynumber "${out_basename}.mpileup" "${out_basename}" --mpileup 1
    
    # Other options:
    # -B
    # --min-coverage 10 \
    # --data-ratio $data_ratio
    # --min-segment-size 20 \
    # --max-segment-size 100

# # 2) Run VarScan2 copyCaller to adjust for GC content and make preliminary calls
#     java -jar $varscan_jar copyCaller \
#     "${results_dir_name}/${out_basename}.copynumber" \
#     --output-file "${results_dir_name}/${out_basename}.copynumber.called" \
#     --output-homdel-file "${results_dir_name}/${out_basename}.copynumber.called.homdel"

#     # See webpage for info on other options:
#     # --min-coverage [8]
#     # --amp-threshold [0.25]
#     # --del-threshold [0.25]
#     # --min-region-size [10]
#     # --recenter-up [0]
#     # --recenter-down [0]

# # 3) Apply circular binary segmentation in R using the library DNAcopy from Bioconductor
#     # TODO: pass $results_dir_name?
#     # Rscript 02--CBS-DNAcopy.r <copy_number_file> <output_basename>
#     Rscript 02--CBS-DNAcopy.r "${out_basename}.copynumber.called" "$out_basename"

# 4) Visualize CBS results and recenter if necessary

echo End Job