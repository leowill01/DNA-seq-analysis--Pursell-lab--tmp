#!/usr/bin/Rscript

# About -------------------------------------------------------------------

# This script takes copy number output from VarScan2 and uses DNAcopy to perform circular binary segmentation (CBS)
# Workflow is adapted from the current VarScan2 webpage: http://dkoboldt.github.io/varscan/copy-number-calling.html#copy-number-segmentation and also the maybe-outdated VarScan2 paper [Koboldt, et. al, 2013, "Using VarScan 2 for Germline Variant Calling and Somatic Mutation Detection"]

# Usage -------------------------------------------------------------------

# Rscript <copy_number_file> <outfile_basename>

# This script is referenced from a shell script that runs samtools mpileup --> varscan copynumber --> varscan copyCaller --> this_script.r
# If run standalone:

# Setup -------------------------------------------------------------------

# Load packages
library('DNAcopy')

# TESTING
    # setwd("/Volumes/Pursell-Lab-HD/Pursell-lab/projects/project_02--DNA-seq-analysis-POLE-mutant-tumors/04--analysis/expt_07--copy-number-variation-calling/test-m079") # to run locally. leave blank to run on server.
    cli_args = c("m079-TN-CNV.copynumber.called", "m079-TN-CNV")
    cn_file = cli_args[1]
    out_basename = cli_args[2]

# Get args passed from command line
# cli_args = commandArgs(trailingOnly = T)
# cn_file = cli_args[1]
# out_basename = cli_args[2]

# Set up results/output folder
results_dir = paste0("results--", as.character(basename(cli_args[1])), format(Sys.time(), '--%Y-%m-%d-%H%M%S'))
results_plots_dir = paste0(results_dir, "/plots/")
results_tables_dir = paste0(results_dir, "/tables/")
dir.create(results_dir)
dir.create(results_plots_dir)
dir.create(results_tables_dir)

# Execution Code ----------------------------------------------------------

# Import copy number data
cn = read.table(file = cn_file, header = T)

# Create copynumber object using the adjusted log2 ratio, chromosome, and region start position from the input
CNA.object = CNA(genomdat = cn[,6], 
                 chrom = cn[,1], 
                 maploc = cn[,2], 
                 data.type = 'logratio')
# # From [Koboldt, et. al, 2013]
# CNA.object = CNA(genomdat = cn$adjusted_log_ratio, 
#                  chrom = cn$chrom, 
#                  maploc = cn$chr_start, 
#                  data.type = 'logratio')

# Perform smoothing on the CNA object to remove outliers
CNA.smoothed = smooth.CNA(CNA.object)

# Run CBS on the smoothed data
segs = segment(CNA.smoothed, verbose = 0, min.width = 2, undo.SD = 3) 
# Note that the undo.SD parameter follows recommendations from the authors of the DNAcopy package, who suggest that change-points of less than two or three standard deviations should be removed. The user is encouraged to experiment with different undo.SD values to achieve a level of segmentation that is best fitted to the data. Empirical experience suggests that undo.SD values between 1 and 4 are usually appropriate.

# Calculate change-point p-values for segments
p_segment = segments.p(segs)

# Plot the results
plot(segs, type = "w")
# Note that the DNAcopy library overrides the plot() function in R, to provide appropriate plotting of raw data and segmented results.

# Output the segmented regions to a tab-delimited file
# From web workflow:
segs2 = segs$output
write.table(segs2[,2:6], 
            file = paste0(results_tables_dir, out_basename, "-DNAcopy-segs2.txt"), 
            row.names = F, col.names = T, quote = F, sep = "\t")
# From [Koboldt, et. al, 2013]:
write.table(p_segment, 
            file = paste0(results_tables_dir, out_basename, "-DNAcopy-p_segment.txt"), sep = "\t", quote = F, col.names = T, row.names = F)
