# #!/usr/bin/env Rscript  

# # Introduction ##############################

# # 1) About --------------------------------------------------------------

# # This script take as input a single VCF file and uses the R package MutationalPatterns to analyze their mutation spectra.

# # TODO: Implement taking VCF file arg from command line (after script function is done)
# # TODO: indicate POLE-exo-mut label somewhere

# # 2) Usage --------------------------------------------------------------

# # Rscript mutSigs.r [VCF_FILE] [SAMPLE_OR_COMPARISON_NAME] [TISSUE_TYPES] [MATCHED_NORMALS]
#     # Where:
#     # VCF_FILE = path to single VCF file
#     # SAMPLE_OR_COMPARISON_NAME = string
#     # TSSUE_TYPES = string
#     # MATCHED_NORMALs = string

# # Run from within the sample/comparison folder in the "expt_06/" folder.
# # e.g. "project_03/expt_06/sample_01/mutSigs.r"

# # FIXME: replace with running script from the "expt_06/{sample|comparison}" dir. Should be able to delete this.
# test_wd = "/Users/leo/Documents/documents-Leo/03--work--Pursell-lab/03--projects--work--Pursell-lab/project_002--DNA-seq-analysis-of-mouse-POLE-exo-mutant-tumor-data/04--analysis/expt_06--SNP-SNV-mutation-signature-analysis/m079-test"
# setwd(test_wd)

# # Setup ----------------------------------------------------------------

# # 1) Load packages ----

# library("BiocManager")
# library("BSgenome")
# library("MutationalPatterns")
# library("gridExtra")
# library("NMF")
# library("ggplot2")
# # Download and load reference genome (mm10)
# # Make sure the genome is installed first (use "install.packages()" to install)
# mm10RefGenome = "BSgenome.Mmusculus.UCSC.mm10"
# library(mm10RefGenome,
#         character.only = T)

# # 2) Arguments ----------------------------------------------------------

# # Take args from the command line
# # args = commandArgs(trailingOnly = T) 
# # TODO: Replace with this ^ when taking commands from CLI
# args = ""
# args[1] = "00--INPUT/m079-tn-WES.snp.Somatic.hc.filter.mm10_multianno.vcf"
# args[2] = "m079--tumor"
# args[3] = "spleen"
# args[4] = "tail"

# # Set variables according to command line args
# # FIXME: Change to accept any VCF file passed from the command line
# vcf_files = args[1]
# sample_names = args[2]
# tissue_types = args[3]
# matched_normals = args[4]

# # 3) Make results directories -------------------------------------------

# # Make results dir for VCF file being analyzed
# results_dir = paste0("results--", as.character(basename(vcf_files)), format(Sys.time(), '--%Y-%m-%d-%H%M%S'))
# results_plots_dir = paste0(results_dir, "/plots/")
# results_tables_dir = paste0(results_dir, "/tables/")
# dir.create(results_dir)
# dir.create(results_plots_dir)
# dir.create(results_tables_dir)

# # Analysis  --------------------------------------------------------------------

# # Adapted from "Introduction to MutationalPatters" guide (Blokzijl, 2019)

# # 1) Load reference genome ----------------------------------------------
# # [in section "Setup: 2) Arguments"]

# # 2) Load sample data ---------------------------------------------------

# # Load VCF files into a GRangesList
# vcfs = read_vcfs_as_granges(vcf_files = vcf_files, 
#                             sample_names = sample_names, 
#                             mm10RefGenome)
# summary(vcfs)

# # 3) Mutation Characteristics ---------------------------------------------

# # 3.1) Base substitution types ----------------------------------------------

# # Retrieve base substitutions from VCF GRanges object as "REF>ALT"
# muts = mutations_from_vcf(vcfs[[1]])

# # Convert mutations to the 6 conventional substitution types
# types = mut_type(vcfs[[1]])

# # Retrieve sequence context of the base substitutions in the VCF object from the reference genome
# contexts = mut_context(vcfs[[1]], mm10RefGenome)

# # Retrieve the types and contexts for all positions in the VCF GRanges object. For base substitutions that are converted to conventional 6 base subs, the reverse complement is returned
# type_contexts = type_context(vcfs[[1]], mm10RefGenome)
# lapply(type_contexts, head, 12)

# # Count mutation type occurrences for all VCF objects in the GRanges list. For C>T mutations, a distinction is made between C>T at CpG sites and other sites.
# type_occurrences = mut_type_occurrences(vcfs, mm10RefGenome)

# write.table(type_occurrences, file = paste0(results_tables_dir, "01--mutation-counts.txt"), 
#             sep = '\t', col.names = NA, quote = F)

# # 3.2) Mutation Spectrum ----------------------------------------------------

# # A mutation spectrum shows the relative contribution of each mutation type in the base substitution catalogs. The plot_spectrum function plots the mean relative contribution of each of the 6 base substitution types over all samples. Error bars indicate standard deviation over all samples. The total number of mutations is indicated.

# # Plot mut spectrum 
# p1 = plot_spectrum(type_occurrences = type_occurrences)

# # plot mut spectrum with distinction between C>T at CpG sites and other sites
# p2 = plot_spectrum(type_occurrences = type_occurrences, CT = T)

# # Combine multiple plots
# pdf(file = paste0(results_plots_dir, "01--mut-type-relative-contrib.pdf"),
#     width = 8, height = 4)
# grid.arrange(p1, p2, ncol = 2, widths=c(3,3))
# dev.off()

# # 3.3) 96 Trinucleotide mutation profile ------------------------------------

# # Make a 96 trinucleotide mutation count matrix
# v_mut_matrix = mut_matrix(vcf_list = vcfs, ref_genome = mm10RefGenome)
# write.table(v_mut_matrix, file = paste0(results_tables_dir, "02--mut-trint-contexts.txt"), sep = '\t', col.names = NA, quote = F)
# # TODO: Add this style of trinuc to the VCF info field!

# # Plot 96 profile
# pdf(file = paste0(results_plots_dir, "03--96-mut-spectrum.pdf"),
#     width = 6, height = 2)
# plot_96_profile(v_mut_matrix, condensed = F)
# dev.off()

# # Plot 96 profile - condensed
# pdf(file = paste0(results_plots_dir, "04--96-mut-spectrum-condensed.pdf"),
#     width = 6, height = 2)
# plot_96_profile(v_mut_matrix, condensed = T)
# dev.off()

# 4) Mutational signatures -------------------------------------

# 4.1) De novo mutational signature extraction using NMF ------------------

# Add a small pseudocount to the mutation count matrix
v_mut_matrix_nmf = v_mut_matrix + 0.0001

# NOTE: Only works when comparing multiple VCFs. Refer to the multi-sample script

# 4.2) Find optimal contribution of COSMIC signatures to reconstruct 96 mut profiles --------

# 4.2.1) COSMIC mutational signatures -------------------------------------
# TODO: Update to new versions of COSMIC signatures

# download mutational signatures from the COSMIC website
spUrl = paste("https://cancer.sanger.ac.uk/cancergenome/assets/", "signatures_probabilities.txt", sep = "")
cosmic_signatures = read.table(file = spUrl, sep = '\t', header = T)
# match the order of the mut types to MutationalPatters standard
new_order = match(row.names(v_mut_matrix_nmf), cosmic_signatures$Somatic.Mutation.Type)
# Reorder cancer signatures dataframe
cosmic_signatures = cosmic_signatures[as.vector(new_order),]
# add trinucleotide changes names as row.names
row.names(cosmic_signatures) = cosmic_signatures$Somatic.Mutation.Type
# Keep only 96 contributions of the signatures in matrix
cosmic_signatures = as.matrix(cosmic_signatures[,4:33])

# plot mut profile of POLE, APOBEC, and MMR COSMIC signatures
pdf(file = paste0(results_plots_dir, "25--POLE-APOBEC-MMR-COSMIC-sigs.pdf"), width = 12, height = 8)
plot_96_profile(cosmic_signatures[,c(10,28,2,13,6,15,20,26)], condensed = T)
dev.off()

# Hierarchically cluster the COSMIC signatures based on their similarity with average linkage
hclust_cosmic_sigs = cluster_signatures(cosmic_signatures, method = 'average')

# store signatures in new order
cosmic_order = colnames(cosmic_signatures)[hclust_cosmic_sigs$order]
pdf(file = paste0(results_plots_dir, "26--COSMIC-sigs-cluster.pdf"), width = 8, height = 6)
plot(hclust_cosmic_sigs)
dev.off()

# 4.2.2) Similarity between mut profiles and COSMIC sigs --------

# FIXME: plot this in a single row heatmap
# calculate pairwise cosine similarity between mut profile and COSMIC sigs
cos_sim_to_cosmic_sigs = cos_sim_matrix(v_mut_matrix_nmf, cosmic_signatures)
barplot(cos_sim_to_cosmic_sigs)
image(cos_sim_to_cosmic_sigs)

# Plot heatmap with specified signature order

pdf(file = paste0(results_plots_dir, "26--cos-sim-to-COSMIC-sigs.pdf"), width = 8, height = 6)
plot_cosine_heatmap(cos_sim_to_cosmic_sigs, 
                    col_order = cosmic_order, 
                    cluster_rows = T)
dev.off()

# # fix mut matrix to the COSMIC mut sigs
# fitRes = fit_to_signatures(v_mut_matrix_nmf, cosmic_signatures)
# 
# # plot the optimal contribution of the COSMIC signatures in each sample as a stacked barplot
# # select sigs with some contribution
# selectSigs = which(rowSums(fitRes$contribution) > 10)
# # plot contribution barplot
# pdf(file = paste0(results_plots_dir, "28--COSMIC-sigs-contrib-barplot.pdf"), width = 8, height = 4)
# plot_contribution(fitRes$contribution[selectSigs,],
#                   cosmic_signatures[selectSigs,],
#                   coord_flip = F,
#                   mode = 'absolute')
# dev.off()
# 
# # plot relative contribution of the cancer signatures in each sample as a heatmap with sample clustering
# pdf(file = "results/plots/29--COSMIC-sigs-contrib-heatmap.pdf", width = 10, height = 4)
# plot_contribution_heatmap(fitRes$contribution,
#                           cluster_samples = T,
#                           method = 'complete')
# dev.off()
# 
# # compare the COSMIC-signature-reconstructed mut profile of each sample with its original mut profile
# 
# pdf(file = "results/plots/30--m079-COSMIC-sigs-cosSim.pdf", width = 8, height = 6)
# plot_compare_profiles(v_mut_matrix_nmf[,1],
#                       fitRes$reconstructed[,1],
#                       profile_names = c('Original','Reconstructed'),
#                       condensed = T)
# dev.off()
# 
# pdf(file = "results/plots/31--m084-COSMIC-sigs-cosSim.pdf", width = 8, height = 6)
# plot_compare_profiles(v_mut_matrix_nmf[,2],
#                       fitRes$reconstructed[,2],
#                       profile_names = c('Original','Reconstructed'),
#                       condensed = T)
# dev.off()
# 
# pdf(file = "results/plots/32--m122-COSMIC-sigs-cosSim.pdf", width = 8, height = 6)
# plot_compare_profiles(v_mut_matrix_nmf[,3],
#                       fitRes$reconstructed[,3],
#                       profile_names = c('Original','Reconstructed'),
#                       condensed = T)
# dev.off()
# 
# pdf(file = "results/plots/33--m124-COSMIC-sigs-cosSim.pdf", width = 8, height = 6)
# plot_compare_profiles(v_mut_matrix_nmf[,4],
#                       fitRes$reconstructed[,4],
#                       profile_names = c('Original','Reconstructed'),
#                       condensed = T)
# dev.off()
# 
# pdf(file = "results/plots/34--m157-COSMIC-sigs-cosSim.pdf", width = 8, height = 6)
# plot_compare_profiles(v_mut_matrix_nmf[,5],
#                       fitRes$reconstructed[,5],
#                       profile_names = c('Original','Reconstructed'),
#                       condensed = T)
# dev.off()
# 
# pdf(file = "results/plots/35--m1098-COSMIC-sigs-cosSim.pdf", width = 8, height = 6)
# plot_compare_profiles(v_mut_matrix_nmf[,6],
#                       fitRes$reconstructed[,6],
#                       profile_names = c('Original','Reconstructed'),
#                       condensed = T)
# dev.off()
# 
# 
# # calculate the cosine sim between all original and reconstructed mut profiles
# # calc all pairwise cos similarities
# cosSimOriRec = cos_sim_matrix(v_mut_matrix_nmf, fitRes$reconstructed)
# # extract cosine sim per sample between original and reconstructed
# cosSimOriRec = as.data.frame(diag(cosSimOriRec))
# write.table(x = cosSimOriRec, file = "results/tables/01--COSMIC-cosSim-all-samples.txt", sep = '\t', quote = F)
# 
# # make barplot of cosine sims between originals and reconstructed mut profiles of each sample
# # adjust data frame for plotting w. ggplot
# colnames(cosSimOriRec) = "cos_sim"
# cosSimOriRec$sample = row.names(cosSimOriRec)
# # make barplot
# ggplot(cosSimOriRec, aes(y=cos_sim, x=sample)) +
#     geom_bar(stat='identity', fill = 'skyblue4') +
#     coord_cartesian(ylim=c(0.8, 1)) +
#     # coord _ flip(ylim=c(0.8,1)) +
#     ylab("Cosine similarity\n original VS reconstructed") +
#     xlab("") +
#     # Reverse order of the samples such that first is up
#     # xlim(rev(levels(factor(cos _ sim _ ori _ rec$sample)))) +
#     theme_bw() +
#     theme(panel.grid.minor.y=element_blank(), panel.grid.major.y=element_blank()) +
#     # Add cut.off line
#     geom_hline(aes(yintercept=.95))
# ggsave(filename = "results/plots/36--COSMIC-cosSim-allSamples-barplot.pdf", width = 4, height = 3, dpi = 320)
# 
