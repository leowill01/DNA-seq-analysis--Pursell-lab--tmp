# About ##############################

# This script is called from a bash script, takes an input of a single VCF file, and uses the R package MutationalPatterns to analyze their mutation spectra.
# FIXME: input = single VCF??

# Arguments ##############################

# Get args from those passed to the sbatch script
# args = commandArgs(trailingOnly = TRUE)

# Code ##############################

# 1. Setup ------------------------------

# Adapted from "Introduction to MutationalPatters" guide (Blokzijl, 2019)
# TODO: label samples with tumor tissue for plotting
# TODO: indicate POLE-exo-mut label somewhere

dir.create("results")
dir.create("results/plots")
dir.create("results/tables")

library("BSgenome")
library("MutationalPatterns")
library("gridExtra")
library("NMF")
library("ggplot2")


# # 2. Data --------------------------------------------------------------------

# # 2.1. Load reference genome 

# # List available genomes
# available.genomes()

# # Download and load reference genome (mm10)
# mm10RefGenome = "BSgenome.Mmusculus.UCSC.mm10"
# library(mm10RefGenome, character.only = T)


# # 2.2. Load sample data ---------------------------------------------------

# # Load MutationalPatterns
# #   under intro

# # TODO--Make a vector of VCF files to input
# vcfFiles = c("/Users/leo/Documents/documents-Leo/work/lab--Pursell/lab-work/projects/project_002--DNA-seq-analysis-of-mouse-POLE-exo-mutant-tumor-data/04--analysis-and-results/expt_04_REMOTE--VarScan2-variant-annotation-ANNOVAR/m079-annnovar/m079-tn-WES.clean.snp.Somatic.hc.somaticFilter.mm10_multianno.vcf", 
#               "/Users/leo/Documents/documents-Leo/work/lab--Pursell/lab-work/projects/project_002--DNA-seq-analysis-of-mouse-POLE-exo-mutant-tumor-data/04--analysis-and-results/expt_04_REMOTE--VarScan2-variant-annotation-ANNOVAR/m084-annnovar/m084-tn-WES.clean.snp.Somatic.hc.somaticFilter.mm10_multianno.vcf", 
#               "/Users/leo/Documents/documents-Leo/work/lab--Pursell/lab-work/projects/project_002--DNA-seq-analysis-of-mouse-POLE-exo-mutant-tumor-data/04--analysis-and-results/expt_04_REMOTE--VarScan2-variant-annotation-ANNOVAR/m122-annnovar/m122-tn-WES.clean.snp.Somatic.hc.somaticFilter.mm10_multianno.vcf", 
#               "/Users/leo/Documents/documents-Leo/work/lab--Pursell/lab-work/projects/project_002--DNA-seq-analysis-of-mouse-POLE-exo-mutant-tumor-data/04--analysis-and-results/expt_04_REMOTE--VarScan2-variant-annotation-ANNOVAR/m124-annnovar/m124-tn-WES.clean.snp.Somatic.hc.somaticFilter.mm10_multianno.vcf", 
#               "/Users/leo/Documents/documents-Leo/work/lab--Pursell/lab-work/projects/project_002--DNA-seq-analysis-of-mouse-POLE-exo-mutant-tumor-data/04--analysis-and-results/expt_04_REMOTE--VarScan2-variant-annotation-ANNOVAR/m157-annnovar/m157-tn-WES.clean.snp.Somatic.hc.somaticFilter.mm10_multianno.vcf", 
#               "/Users/leo/Documents/documents-Leo/work/lab--Pursell/lab-work/projects/project_002--DNA-seq-analysis-of-mouse-POLE-exo-mutant-tumor-data/04--analysis-and-results/expt_04_REMOTE--VarScan2-variant-annotation-ANNOVAR/m1098-annnovar/m1098-tn-WES.clean.snp.Somatic.hc.somaticFilter.mm10_multianno.vcf")

# # Define sample names for VCF files
# sampleNames = c("m079", "m084", "m122", "m124", "m157", "m1098")

# # Load VCF files into a GRangesList
# vcfs = read_vcfs_as_granges(vcf_files = vcfFiles, sample_names = sampleNames, mm10RefGenome)
# summary(vcfs)

# # Define relevant metadata on the samples, such as tissue type
# tissue = c("spleen", "thymus", "thymus", "spleen", "thymus", "thymus")
# matchedNormal = c(rep("tail", 6))



# # 3. Mutation Characteristics ---------------------------------------------

# # 3.1. Base substitution types --------------------------------------------

# # Retrieve base substitutions from VCF GRanges object as "REF>ALT"
# muts = mutations_from_vcf(vcfs[[1]])
# head(muts, 12)

# # Convert mutations to the 6 conventional substitution types
# types = mut_type(vcfs[[1]])
# head(types, 12)

# # Retrieve sequence context of the base substitutions in the VCF object from the reference genome
# context = mut_context(vcfs[[1]], mm10RefGenome)
# head(context, 12)

# # Retrieve the types and contexts for all positions in the VCF GRanges object. For base substitutions that are converted to conventional 6 base subs, the reverse complement is returned
# typeContext = type_context(vcfs[[1]], mm10RefGenome)
# lapply(typeContext, head, 12)

# # Count mutation type occurrences for all VCF objects in the GRanges list. For C>T mutations, a distinction is made between C>T at CpG sites and other sites. 
# typeOccurrences = mut_type_occurrences(vcfs, mm10RefGenome)
# typeOccurrences


# # 3.2. Mutation Spectrum --------------------------------------------------
# # A mutation spectrum shows the relative contribution of each mutation type in the base substitution catalogs. The plot _ spectrum function plots the mean relative contribution of each of the 6 base substitution types over all samples. Error bars indicate standard deviation over all samples. The total number of mutations is indicated.

# p1 = plot_spectrum(type_occurrences = typeOccurrences)

# # plot mut spectrum with distinction between C>T at CpG sites and other sites
# p2 = plot_spectrum(type_occurrences = typeOccurrences, CT = T)

# # Plot spectrum without legend
# # p3 = plot_spectrum(type_occurrences = typeOccurrences, CT = T, legend = F)

# # Combine multiple plots
# pdf(file = "results/plots/01--total-muts-rel-contribution.pdf", 
#     width = 8, height = 4)
# grid.arrange(p1, p2, ncol = 2, widths=c(3,3))
# dev.off()

# # Facet the per sample group, e.g. plot spectrum for each tissue separately
# pdf(file = "results/plots/02--total-muts-rel-contribution-by-tissue.pdf", 
#     width = 8, height = 4)
# plot_spectrum(typeOccurrences, by = tissue, CT = T, legend = T)
# dev.off()

# # Define your own 7 colors for spectrum plotting
# # palette = c("pink", "orange", "blue", 'lightblue', 'green', 'red','purple')
# # p5 = plot_spectrum(type_occurrences = typeOccurrences, CT = T, legend = T, colors = palette)
# # grid.arrange(p4, p5, ncol=2, widths=c(4, 2.3))



# # 3.3. 96 Trinucleotide mutation profile ----------------------------------

# # Make a 96 trinucleotide mutation count matrix
# mutMat = mut_matrix(vcf_list = vcfs, ref_genome = mm10RefGenome)
# head(mutMat)

# # Plot the 96 profile of all samples
# # plot_96_profile(mutMat)

# # Plot 96 profile of all samples in a more condensed plotting format
# pdf(file = "results/plots/03--96-spectrum-all-samples.pdf", 
#     width = 8, height = 6)
# plot_96_profile(mutMat, condensed = T)
# dev.off()


# # 4. Mutational Signatures ------------------------------------------------

# # Add a small pseudocount to your mutation count matrix
# mutMatNmf = mutMat + 0.001

# # Use NMF package to generate an estimate rank plot
# estimate = nmf(x = mutMatNmf, rank = 1:6, method = "brunet", nrun=30, seed=123456)
# pdf(file = "results/plots/04--NMF-rank-survey-rank6-nrun30.pdf", width = 8, height = 6)
# plot(estimate)
# dev.off()
# # NOTE: based on this, I choose rank of 3 or 4


# # 4.1 NMF Rank 3 ----------------------------------------------------------
# # Extract 3 mutational signatures from the mutation count matrix with extract_signatures(). (For larger datasets it is wise to perform more iterations by changing the nrun parameter to achieve stability and avoid local minima)
# nmfRes3 = extract_signatures(mut_matrix = mutMatNmf, rank = 3, nrun = 30)

# # assign signature names
# colnames(nmfRes3$signatures) = c('Sig A', 
#                                 'Sig B', 
#                                 'Sig C')
#                                 # 'Sig D')
# rownames(nmfRes3$contribution) = c('Sig A', 
#                                   'Sig B', 
#                                   'Sig C')
#                                   # 'Sig D')

# # Plot the 96 trinucleotide profile of the signatures
# pdf(file = "results/plots/05--nmfr3-mut-sigs.pdf")
# plot_96_profile(nmfRes3$signatures, condensed = T)
# dev.off()

# # visualize the contribution of the signatures in a barplot
# pc1 = plot_contribution(nmfRes3$contribution, nmfRes3$signature, mode = 'relative')

# # visualize the controbution of the signatures in absolute number of mutations
# pc2 = plot_contribution(nmfRes3$contribution, nmfRes3$signature, mode = 'absolute')

# # combine the two plots
# pdf(file = "results/plots/06--nmfr3-sig-contributions.pdf")
# grid.arrange(pc1, pc2)
# dev.off()

# # flip x and y coordinates
# pdf(file = "results/plots/07--nmfr3-sig-contributions-flipxy.pdf")
# plot_contribution(nmfRes3$contribution, nmfRes3$signature, mode = 'absolute', coord_flip = T)
# dev.off()

# # plot signature contribution as a heatmap with sample clustering dendrogram and a specified signature order
# pch1 = plot_contribution_heatmap(nmfRes3$contribution, sig_order = c('Sig A', 'Sig B', 'Sig C'))

# # plot signature contribution as a heatmap without sample clustering
# pch2 = plot_contribution_heatmap(nmfRes3$contribution, cluster_samples = F)

# # combine the two plots
# pdf(file = "results/plots/08--nmfr3-sig-contributions-heatmap.pdf")
# grid.arrange(pch1, pch2, ncol = 2, widths = c(2, 1.6))
# dev.off()

# # compare the reconstructed mutational profile with the original mutational profile for each sample
# pdf(file = "results/plots/09--nmfr3-m079-cosine-reconstruct.pdf")
# plot_compare_profiles(mutMatNmf[,1], 
#                       nmfRes3$reconstructed[,1], 
#                       profile_names = c('Original', 'Reconstructed'), 
#                       condensed = T)
# dev.off()
# pdf(file = "results/plots/10--nmfr3-m084-cosine-reconstruct.pdf")
# plot_compare_profiles(mutMatNmf[,2], 
#                       nmfRes3$reconstructed[,2], 
#                       profile_names = c('Original', 'Reconstructed'), 
#                       condensed = T)
# dev.off()
# pdf(file = "results/plots/11--nmfr3-m122-cosine-reconstruct.pdf")
# plot_compare_profiles(mutMatNmf[,3], 
#                       nmfRes3$reconstructed[,3], 
#                       profile_names = c('Original', 'Reconstructed'), 
#                       condensed = T)
# dev.off()
# pdf(file = "results/plots/12--nmfr3-m124-cosine-reconstruct.pdf")
# plot_compare_profiles(mutMatNmf[,4], 
#                       nmfRes3$reconstructed[,4], 
#                       profile_names = c('Original', 'Reconstructed'), 
#                       condensed = T)
# dev.off()
# pdf(file = "results/plots/13--nmfr3-m157-cosine-reconstruct.pdf")
# plot_compare_profiles(mutMatNmf[,5], 
#                       nmfRes3$reconstructed[,5], 
#                       profile_names = c('Original', 'Reconstructed'), 
#                       condensed = T)
# dev.off()
# pdf(file = "results/plots/14--nmfr3-m1098-cosine-reconstruct.pdf")
# plot_compare_profiles(mutMatNmf[,6], 
#                       nmfRes3$reconstructed[,6], 
#                       profile_names = c('Original', 'Reconstructed'), 
#                       condensed = T)
# dev.off()



# # 4.2 NMF Rank 4 ----------------------------------------------------------
# # Do the same thing for NMF rank of 4
# nmfRes4 = extract_signatures(mut_matrix = mutMatNmf, rank = 4, nrun = 30)
# colnames(nmfRes4$signatures) = c('Sig A','Sig B','Sig C','Sig D')
# rownames(nmfRes4$contribution) = c('Sig A','Sig B','Sig C','Sig D')

# pdf(file = "results/plots/15--nmfr4-96tnt-mut-sigs.pdf")
# plot_96_profile(nmfRes4$signatures, condensed = T)
# dev.off()

# pc1 = plot_contribution(nmfRes4$contribution, nmfRes4$signature, mode = 'relative')
# pc2 = plot_contribution(nmfRes4$contribution, nmfRes4$signature, mode = 'absolute')

# pdf(file = "results/plots/16--nmfr4-sig-contributions.pdf")
# grid.arrange(pc1, pc2)
# dev.off()

# pdf(file = "results/plots/17--nmfr4-sig-contributions-flipxy.pdf")
# plot_contribution(nmfRes4$contribution, nmfRes4$signature, mode = 'absolute', coord_flip = T)
# dev.off()

# pch1 = plot_contribution_heatmap(nmfRes4$contribution, sig_order = c('Sig A', 'Sig B', 'Sig C', 'Sig D'))
# pch2 = plot_contribution_heatmap(nmfRes4$contribution, cluster_samples = F)

# pdf(file = "results/plots/18--nmfr4-sig-contributions-heatmap.pdf")
# grid.arrange(pch1, pch2, ncol = 2, widths = c(2, 1.6))
# dev.off()
# # cosine similarity reconstuction
# pdf(file = "results/plots/19--nmfr4-m079-cosine-reconstruct.pdf")
# plot_compare_profiles(mutMatNmf[,1], 
#                       nmfRes4$reconstructed[,1], 
#                       profile_names = c('Original', 'Reconstructed'), 
#                       condensed = T)
# dev.off()

# pdf(file = "results/plots/20--nmfr4-m084-cosine-reconstruct.pdf")
# plot_compare_profiles(mutMatNmf[,2], 
#                       nmfRes4$reconstructed[,2], 
#                       profile_names = c('Original', 'Reconstructed'), 
#                       condensed = T)
# dev.off()

# pdf(file = "results/plots/21--nmfr4-m122-cosine-reconstruct.pdf")
# plot_compare_profiles(mutMatNmf[,3], 
#                       nmfRes4$reconstructed[,3], 
#                       profile_names = c('Original', 'Reconstructed'), 
#                       condensed = T)
# dev.off()

# pdf(file = "results/plots/22--nmfr4-m124-cosine-reconstruct.pdf")
# plot_compare_profiles(mutMatNmf[,4], 
#                       nmfRes4$reconstructed[,4], 
#                       profile_names = c('Original', 'Reconstructed'), 
#                       condensed = T)
# dev.off()

# pdf(file = "results/plots/23--nmfr4-m157-cosine-reconstruct.pdf")
# plot_compare_profiles(mutMatNmf[,5], 
#                       nmfRes4$reconstructed[,5], 
#                       profile_names = c('Original', 'Reconstructed'), 
#                       condensed = T)
# dev.off()

# pdf(file = "results/plots/24--nmfr4-m1098-cosine-reconstruct.pdf")
# plot_compare_profiles(mutMatNmf[,6], 
#                       nmfRes4$reconstructed[,6], 
#                       profile_names = c('Original', 'Reconstructed'), 
#                       condensed = T)
# dev.off()



# # 4.2. Find optimal contribution of known signatures ------------------------
# # 4.2.1. COSMIC mutational signatures -------------------------------------

# # download mutational signatures from the COSMIC website
# spUrl = paste("https://cancer.sanger.ac.uk/cancergenome/assets/", "signatures_probabilities.txt", sep = "")
# cancerSignatures = read.table(file = spUrl, sep = '\t', header = T)
# # match the order of the mut types to MutationalPatters standard
# newOrder = match(row.names(mutMatNmf), cancerSignatures$Somatic.Mutation.Type)
# # Reorder cancer signatures dataframe
# cancerSignatures = cancerSignatures[as.vector(newOrder),]
# # add trinucleotide changes names as row.names
# row.names(cancerSignatures) = cancerSignatures$Somatic.Mutation.Type
# # Keep only 96 contributions of the signatures in matrix
# cancerSignatures = as.matrix(cancerSignatures[,4:33])

# # plot mut profile of POLE, APOBEC, and MMR COSMIC signatures
# pdf(file = "results/plots/25--relevant-COSMIC-sigs.pdf", width = 12, height = 8)
# plot_96_profile(cancerSignatures[,c(10,28,2,13,6,15,20,26)], condensed = T)
# dev.off()

# # Hierarchically cluster the COSMIC signatures based on their similarity with average linkage
# hclustCosmic = cluster_signatures(cancerSignatures, method = 'average')
# # store signatures in new order
# cosmicOrder = colnames(cancerSignatures)[hclustCosmic$order]
# pdf(file = "results/plots/26--COSMIC-sigs-cluster.pdf", width = 8, height = 6)
# plot(hclustCosmic)
# dev.off()


# # 4.2.2. Similarity between mut profiles and COSMIC sigs --------

# # calculate cosine similarity between two mut profiles/sigs 
# cos_sim(mutMatNmf[,1], cancerSignatures[,1])
# # TODO: do these specific comparisons for each sample vs desired sigs? or is this redundant bc of next steps?

# # calculate pairwise cosine similarity between mut profiles and COSMIC sigs
# cosSimSamplesSigs = cos_sim_matrix(mutMatNmf, cancerSignatures)
# # plot heatmap with specified signature order
# pdf(file = "results/plots/27--cosSim-heatmap-COSMIC-sigs.pdf", width = 8, height = 4)
# plot_cosine_heatmap(cosSimSamplesSigs, col_order = cosmicOrder, cluster_rows = T)
# dev.off()

# # 4.2.3. Find optimal contribution of COSMIC signatures to reconstruct 96 mut profiles --------

#  # fix mut matrix to the COSMIC mut sigs
# fitRes = fit_to_signatures(mutMatNmf, cancerSignatures)

# plot the optimal contribution of the COSMIC signatures in each sample as a stacked barplot
# select sigs with some contribution
selectSigs = which(rowSums(fitRes$contribution) > 10)
# plot contribution barplot
pdf(file = "results/plots/28--COSMIC-sigs-contrib-barplot.pdf", width = 8, height = 4)
plot_contribution(fitRes$contribution[selectSigs,],
                  cancerSignatures[selectSigs,],
                  coord_flip = F,
                  mode = 'absolute')
dev.off()

# # plot relative contribution of the cancer signatures in each sample as a heatmap with sample clustering
# pdf(file = "results/plots/29--COSMIC-sigs-contrib-heatmap.pdf", width = 10, height = 4)
# plot_contribution_heatmap(fitRes$contribution, 
#                           cluster_samples = T, 
#                           method = 'complete')
# dev.off()

# # compare the COSMIC-signature-reconstructed mut profile of each sample with its original mut profile

# pdf(file = "results/plots/30--m079-COSMIC-sigs-cosSim.pdf", width = 8, height = 6)
# plot_compare_profiles(mutMatNmf[,1], 
#                       fitRes$reconstructed[,1], 
#                       profile_names = c('Original','Reconstructed'), 
#                       condensed = T)
# dev.off()

# pdf(file = "results/plots/31--m084-COSMIC-sigs-cosSim.pdf", width = 8, height = 6)
# plot_compare_profiles(mutMatNmf[,2], 
#                       fitRes$reconstructed[,2], 
#                       profile_names = c('Original','Reconstructed'), 
#                       condensed = T)
# dev.off()

# pdf(file = "results/plots/32--m122-COSMIC-sigs-cosSim.pdf", width = 8, height = 6)
# plot_compare_profiles(mutMatNmf[,3], 
#                       fitRes$reconstructed[,3], 
#                       profile_names = c('Original','Reconstructed'), 
#                       condensed = T)
# dev.off()

# pdf(file = "results/plots/33--m124-COSMIC-sigs-cosSim.pdf", width = 8, height = 6)
# plot_compare_profiles(mutMatNmf[,4], 
#                       fitRes$reconstructed[,4], 
#                       profile_names = c('Original','Reconstructed'), 
#                       condensed = T)
# dev.off()

# pdf(file = "results/plots/34--m157-COSMIC-sigs-cosSim.pdf", width = 8, height = 6)
# plot_compare_profiles(mutMatNmf[,5], 
#                       fitRes$reconstructed[,5], 
#                       profile_names = c('Original','Reconstructed'), 
#                       condensed = T)
# dev.off()

# pdf(file = "results/plots/35--m1098-COSMIC-sigs-cosSim.pdf", width = 8, height = 6)
# plot_compare_profiles(mutMatNmf[,6], 
#                       fitRes$reconstructed[,6], 
#                       profile_names = c('Original','Reconstructed'), 
#                       condensed = T)
# dev.off()


# calculate the cosine sim between all original and reconstructed mut profiles
# calc all pairwise cos similarities
cosSimOriRec = cos_sim_matrix(mutMatNmf, fitRes$reconstructed)
# extract cosine sim per sample between original and reconstructed
cosSimOriRec = as.data.frame(diag(cosSimOriRec))
write.table(x = cosSimOriRec, file = "results/tables/01--COSMIC-cosSim-all-samples.txt", sep = '\t', quote = F)

# make barplot of cosine sims between originals and reconstructed mut profiles of each sample
# adjust data frame for plotting w. ggplot
colnames(cosSimOriRec) = "cos_sim"
cosSimOriRec$sample = row.names(cosSimOriRec)
# make barplot
ggplot(cosSimOriRec, aes(y=cos_sim, x=sample)) +
    geom_bar(stat='identity', fill = 'skyblue4') +
    coord_cartesian(ylim=c(0.8, 1)) +
    # coord _ flip(ylim=c(0.8,1)) +
    ylab("Cosine similarity\n original VS reconstructed") +
    xlab("") +
    # Reverse order of the samples such that first is up
    # xlim(rev(levels(factor(cos _ sim _ ori _ rec$sample)))) +
    theme_bw() +
    theme(panel.grid.minor.y=element_blank(), panel.grid.major.y=element_blank()) +
    # Add cut.off line
    geom_hline(aes(yintercept=.95))
ggsave(filename = "results/plots/36--COSMIC-cosSim-allSamples-barplot.pdf", width = 4, height = 3, dpi = 320)
