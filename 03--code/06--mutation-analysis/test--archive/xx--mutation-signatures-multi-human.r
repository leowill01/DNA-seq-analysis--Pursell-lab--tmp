#!/usr/bin/env Rscript

# title: mutational signature analysis - multi
# author: Leo Williams
# date: 2019-08-08

# About -------------------------------------------------------------------

# This script takes a group of VCF files as input and analyzes their mutation catalogues with the R package MutationalPatterns and NMF to extract de novo mutational signatures and to compare against known COSMIC mutation signatures.

# Usage: ====

# Rscript mutsigs.r [INFO_TABLE]

# Where INFO_TABLE is a csv listing the vcf filepaths, their sample names, and other information like tissue of origin, matched normal, etc. This is dependent on the experimental design.

# Setup -------------------------------------------------------------------

# 1) Set working directory ====
# Run script within sample-specific dir within the experiment dir
# e.g. "project_02/expt_06/sample_01/mut-sigs.r"

# 2) Enable testing or production mode ====
    # 2.1) Enable for production:
    # args = commandArgs(trailingOnly = T)

    # 2.2) Enable for testing:
    setwd("/Users/leo/Documents/documents-Leo/02--work--career/lab--Pursell/03--projects--Pursell-lab/project_002--DNA-seq-analysis-of-mouse-POLE-exo-mutant-tumor-data/04.02--analysis--LOCAL/expt_06--SNP-SNV-mutation-signature-analysis/multi/GDC-human-POLE-P286R-samples")
    args = "" # init test command line rgs
    args[1] = "input-files/02--multi-sample-info-table.csv" # multisample info table

# 3) Load info table  ====  
input_info_table = read.csv(file = args[1], header = T)

# 4) Load requirements ====

library('BSgenome')
library('MutationalPatterns')
library('gridExtra')
library('NMF')
library('ggplot2')
library('RColorBrewer')

# 5) Make results dir ====

results_dir = paste0("results--", as.character(basename(args[1])), format(Sys.time(), '--%Y-%m-%d-%H%M%S'))
results_plots_dir = paste0(results_dir, "/plots/")
results_tables_dir = paste0(results_dir, "/tables/")
dir.create(results_dir)
dir.create(results_plots_dir)
dir.create(results_tables_dir)
    
# Analysis ----------------------------------------------------------------

# 1) Introduction ---------------------------------------------------------

# Adapted from the guide "Introduction to MutationalPatterns" (2019, Blokzijl)
# This analysis is for a single VCF file. Another script is available for analyzing pooled multiple VCFs in the same run.

# 2) Data -----------------------------------------------------------------

    # 2.1) List reference genome ====
    
    # List available genomes
    # available.genomes()

    # NOTE: MUST INSTALL GENOME AS ITS OWN PACKAGE BEFORE USAGE USING:
    # if (!requireNamespace("BiocManager", quietly = TRUE))
    #   install.packages("BiocManager")
    # 
    # BiocManager::install("GENOME.NAME.FROM.PREVIOUS.COMMAND")
    
    # Download and load reference genome of interest
    ref_genome = "BSgenome.Hsapiens.UCSC.hg38"
    library(ref_genome, 
            character.only = T)
    
    # 2.2) Load example data ====
    
    vcf_files = as.character(input_info_table$vcf_filepath)
    sample_names = as.character(input_info_table$sample_name)
    vcfs = read_vcfs_as_granges(vcf_files, sample_names, ref_genome)
    
    summary(vcfs)
    
    # Define relevant sample metadata
    tissue = as.character(input_info_table$tissue_of_origin)
    
    # 2.3) Add reference genome contig lengths to vcfs (originally my VCFs didn't have this. Solved by Francis Blokzijl) ====

    # Get and store contig lengths
    chr_lengths = lengths(get(ref_genome))

    # Add chromosome names
    names(chr_lengths) <- names(get(ref_genome))

    # Check which chromosomes are in your vcfs object
    chr_select <- names(seqlengths(vcfs))

    # Get lengths of these chromosomes from reference genome seq lengths and add them to your vcfs object. This fixes the problem with rainfall plot in section 6.1 not showing up
    seqlengths(vcfs) <- chr_lengths[chr_select]

# 3) Mutation characteristics ---------------------------------------------

    # 3.1) Base substitution types ====

    # Retrieve base substitutions from the VCF GRanges object as "REF>ALT"
    muts = mutations_from_vcf(vcfs[[1]])
    head(muts)

    # Retrieve the base substitutions from the VCF GRanges object and convert them to the 6 types of base substitution types that are distinguished by convention: C>A, C>G, C>T, T>A, T>C, T>G.
    types = mut_type(vcfs[[1]])
    head(types)

    # Retrieve the trinucleotide mutation sequence
    context = mut_context(vcfs[[1]], ref_genome)
    head(context)

    # Retrieve the types and contexts for all positions in the VCF GRanges object. For the base substitutions that are converted to the conventional base substitution types, the reverse complement of the sequence context is returned.
    type_context = type_context(vcfs[[1]],
                                ref_genome = ref_genome)
    lapply(type_context, head)
    head(type_context)
    # TODO: EXPORT THIS AS TRINUCLEOTIDE CONTEXT! ("$type" + "$context")
    # > paste0(type_context$types[[1]], ":", type_context$context[[1]])
    # [1] "T>C:GTT"

    # Count mutation type occurrences for all VCF objects in the GRangesList
    type_occurrences = mut_type_occurrences(vcf_list = vcfs, ref_genome)
    head(type_occurrences)
    # Save mutation type occurrences as a table
    write.table(x = type_occurrences,
                file = paste0(results_tables_dir,
                              "01--mutation-type-occurrences.txt"),
                sep = "\t", col.names = NA, quote = F)

    # 3.2) Mutation spectrum ====

    # Plot the mutation spectrum with distinction between C>T at CpG sites and other sites
    p1 = plot_spectrum(type_occurrences = type_occurrences, CT = T)

    # Save combined mutation spectra
    pdf(file = paste0(results_plots_dir, "01--mut-spectrum-all.pdf"), width = 8, height = 4)
    p1
    dev.off()
    
    # Plot the mutation spectra faceted by metadata e.g. tissue
    p2 = plot_spectrum(type_occurrences, by = tissue, CT = T)
  
    # Save mutation spectra by tissue
    pdf(file = paste0(results_plots_dir, "02--mut-spectrum-by-tissue.pdf"), width = 8, height = 4)
    p2
    dev.off()
    
    # 3.3) 96 mutational profile ====

    # Make a 96 trinucleodide mutation count matrix:
    mut_mat = mut_matrix(vcf_list = vcfs, ref_genome = ref_genome)
    head(mut_mat)
    write.table(mut_mat,
                file = paste0(results_tables_dir, "02--96-mut-matrix-all.txt"),
                sep = '\t', col.names = NA, quote = F)

    # Plot the 96 profile of all samples
    pdf(file = paste0(results_plots_dir, "03--96-mut-profiles-all.pdf"),
        width = 12, height = 6)
    plot_96_profile(mut_mat)
    dev.off()
    
    # Plot the 96 profile of each sample
    # Separate mut_mat into each sample
    for (i in 1:length(colnames(mut_mat))) {
      col_sample_name = colnames(mut_mat)[i]
      sample_mat = matrix(mut_mat[,i])
      rownames(sample_mat) = rownames(mut_mat)
      colnames(sample_mat) = col_sample_name
      pdf(file = paste0(results_plots_dir, "03--96-mut-profile-", col_sample_name, ".pdf"), 
          width = 8, height = 2)
      print(plot_96_profile(sample_mat))
      dev.off()
    }
    
# 4) Mutational signatures ------------------------------------------------

    # 4.1) De novo mutational signature extraction using NMF ====
    # Mutational signatures can be extracted from your mutation count matrix 
    
    # First add small pseudocount to mutation count matrix
    mut_mat = mut_mat + 0.0001
    
    # Use the NMF package to generate an estimate rank plot:
    # From "An introduction to NMF package" (Gaujoux, 2018):
    
    # A critical parameter in NMF is the factorization rank r. It deﬁnes the number of metagenes used to approximate the target matrix. Given a NMF method and the target matrix, a common way of deciding on r is to try diﬀerent values, compute some quality measure of the results, and choose the best value according to this quality criteria.
    # 
    # Several approaches have then been proposed to choose the optimal value of r. For example, (Brunet2004) proposed to take the ﬁrst value of r for which the cophenetic coeﬃcient starts decreasing, (Hutchins2008) suggested to choose the ﬁrst value where the RSS curve presents an inﬂection point, and (Frigyesi2008) considered the smallest value at which the decrease in the RSS is lower than the decrease of the RSS obtained from random data.
    
    # When the seeding method is stochastic, multiple runs are usually required to achieve stability or a resonable result. This can be done by setting argument nrun to the desired value. For performance reason we use nrun=5 here, but a typical choice would lies between 100 and 200:
    
    estimate = nmf(mut_mat,
                   rank = 1:5,
                   method = "brunet",
                   nrun = 200,
                   seed = 123456)
    
    # Plot NMF rank estimate plots
    pdf(file = paste0(results_plots_dir, "04--nmf-rank-estimate.pdf"),
        width = 12, height = 8)
    plot(estimate)
    dev.off()
    
    # Extract _ signatures from mutation count matrix based on NMF factorization rank 
    # Rank = 3
    nmf_res_3 = extract_signatures(mut_mat, rank = 3, nrun = 200)
    de_novo_sigs_3 = c('A','B','C')
    colnames(nmf_res_3$signatures) = de_novo_sigs_3
    rownames(nmf_res_3$contribution) = de_novo_sigs_3
    # Plot 96 profile of de novo signatures
    pdf(file = paste0(results_plots_dir, "05--3-de-novo-signatures.pdf"),
        width = 12, height = 8)
    plot_96_profile(nmf_res_3$signatures)
    dev.off()
    
    # Rank = 4
    nmf_res_4 = extract_signatures(mut_mat, rank = 4, nrun = 200)
    de_novo_sigs_4 = c('A','B','C','D')
    colnames(nmf_res_4$signatures) = de_novo_sigs_4
    rownames(nmf_res_4$contribution) = de_novo_sigs_4
    # Plot 96 profile of de novo signatures
    pdf(file = paste0(results_plots_dir, "05--4-de-novo-signatures.pdf"),
        width = 12, height = 8)
    plot_96_profile(nmf_res_4$signatures)
    dev.off()
    
    # Visualize contributions of the 3 signatures in a barplot, heatmap, and reconstruction
    # Barplots
    pc1 = plot_contribution(nmf_res_3$contribution, nmf_res_3$signatures, mode = 'relative')
    pc2 = plot_contribution(nmf_res_3$contribution, nmf_res_3$signatures, mode = 'absolute')
    # Save as PDF
    pdf(file = paste0(results_plots_dir, "06--3-sigs-contribution-all.pdf"),
        width = 6, height = 4)
    grid.arrange(pc1, pc2)
    dev.off()
    # Heatmaps
    pch1 = plot_contribution_heatmap(nmf_res_3$contribution)
    pch2 = plot_contribution_heatmap(nmf_res_3$contribution, cluster_samples = F)
    # Save as PDF
    pdf(file = paste0(results_plots_dir, "06--3-sigs-contribution-heatmap-all.pdf"),
        width = 6, height = 4)
    grid.arrange(pch1, pch2)
    dev.off()
    # Reconstructions
    for (i in 1:length(colnames(mut_mat))) {
      pdf(file = paste0(results_plots_dir, "06--3-sigs-reconstructed-", colnames(mut_mat)[i], ".pdf"), width = 12, height = 6)
      print(plot_compare_profiles(mut_mat[,i], 
                            nmf_res_3$reconstructed[,i], 
                            profile_names = c('Original','Reconstructed')))
      dev.off()
    }
    
    # Visualize contributions of the 4 signatures in a barplot and heatmap
    pc3 = plot_contribution(nmf_res_4$contribution, nmf_res_4$signatures, mode = 'relative')
    pc4 = plot_contribution(nmf_res_4$contribution, nmf_res_4$signatures, mode = 'absolute')
    # Save as PDF
    pdf(file = paste0(results_plots_dir, "07--4-sigs-contribution-all.pdf"),
        width = 6, height = 4)
    grid.arrange(pc3, pc4)
    dev.off()
    # Heatmaps
    pch3 = plot_contribution_heatmap(nmf_res_4$contribution)
    pch4 = plot_contribution_heatmap(nmf_res_4$contribution, cluster_samples = F)
    # Save as PDF
    pdf(file = paste0(results_plots_dir, "07--4-sigs-contribution-heatmap-all.pdf"),
        width = 3, height = 4)
    grid.arrange(pch3, pch4)
    dev.off()
    # Reconstructions
    for (i in 1:length(colnames(mut_mat))) {
      pdf(file = paste0(results_plots_dir, "07--4-sigs-reconstructed-", colnames(mut_mat)[i], ".pdf"), width = 12, height = 6)
      print(plot_compare_profiles(mut_mat[,i], 
                                  nmf_res_4$reconstructed[,i], 
                                  profile_names = c('Original','Reconstructed')))
      dev.off()
    }
    

    # 4.2) Find optimal contribution of known signatures ====

      # 4.2.1) COSMIC mutational signatures ####
  
      # Load the COSMIC Mutational Signatures v3 - single-base substitution (SBS) signatures:
    cosmicv3_sigs_no_arts = read.csv(file = "/Users/leo/Documents/documents-Leo/02--work--career/lab--Pursell/03--projects--Pursell-lab/project_002--DNA-seq-analysis-of-mouse-POLE-exo-mutant-tumor-data/02--data/02--reference-data/COMIC-mutation-signatures/COSMIC-mutation-signatures-v3-no-artif-formatted.csv", header = T)
      # Old - for COSMICv2:
      # Download mutational signatures from COSMIC website:
      # sp_url = paste("https://cancer.sanger.ac.uk/cancergenome/assets/",
      #                "signatures_probabilities.txt",
      #                sep = "")
      # cosmicv3_sigs_no_arts = read.table(sp_url, sep = "\t", header = T)
  
      # Match the order of the mutation types to MutationalPatterns standard
      new_order = match(row.names(mut_mat), cosmicv3_sigs_no_arts$Somatic.Mutation.Type)
      # Reorder cancer signatures dataframe
      cosmicv3_sigs_no_arts = cosmicv3_sigs_no_arts[as.vector(new_order), ]
      # Add trinucleotide changes names as row.names
      row.names(cosmicv3_sigs_no_arts) = cosmicv3_sigs_no_arts$Somatic.Mutation.Type
      # Keep only 96 contributions of the signatures in matrix
      cosmicv3_sigs_no_arts = as.matrix(cosmicv3_sigs_no_arts[ ,4:52])

      # Plot mutational profile of some specified COSMIC signatures:
      # POLE and DNA damage/repair-related signatures from COSMICv3:
      relevant_signatures = c('SBS28', 'SBS10a', 'SBS10b', 
                              'SBS14', 'SBS20', 
                              'SBS6', 'SBS15', 'SBS26', 'SBS44', 
                              'SBS1', 'SBS2', 'SBS13', 'SBS84', 'SBS85', 
                              'SBS3', 'SBS9', 'SBS30', 'SBS36')
      pdf(file = paste0(results_plots_dir, "relevant-COSMICv3-signatures.pdf"), 
          width = 12, height = 16)
      plot_96_profile(cosmicv3_sigs_no_arts[ , relevant_signatures])
      dev.off()
      
      # Hierarchically cluster the COSMIC signatures based on their similarity with average linkage:
      hclust_cosmic = cluster_signatures(cosmicv3_sigs_no_arts, method = "average")
      # Store signatures in new order:
      cosmic_order = colnames(cosmicv3_sigs_no_arts)[hclust_cosmic$order]
      pdf(file = paste0(results_plots_dir, "COSMICv3-signatures-hcluster.pdf"),
          width = 12, height = 8)
      plot(hclust_cosmic)
      dev.off()

      # 4.2.2) Similarity between mutational proﬁles and COSMIC signatures ####

      # Calculate cosine similarity of sample mutation matrix with one COSMIC signature. In this case, Signature 10:
      cos_sim_sig_10 = cos_sim(mut_mat[ ,1], cosmicv3_sigs_no_arts[ ,14])
      cos_sim_sig_10
  
      # Calculate pairwise cosine similarity between the sample mutational proﬁle and COSMIC signatures:
      cos_sim_samples_signatures = cos_sim_matrix(mut_mat, cosmicv3_sigs_no_arts)
  
      # Plot heatmap with specified signature order
      pdf(file = paste0(results_plots_dir, "08--cos-sim-sample-COSMICv3-sigs.pdf"), 
          width = 20, height = 10)
      plot_cosine_heatmap(cos_sim_samples_signatures, 
                          col_order = cosmic_order, 
                          cluster_rows = T)
      dev.off()

    # 4.2.3) Find optimal contribution of COSMIC signatures to reconstruct 96 mutational proﬁles ####

    # Fit mutation matrix to the COSMIC mutational signatures:
    fit_res = fit_to_signatures(mut_mat, cosmicv3_sigs_no_arts)

    # Plot the optimal contribution of the COSMIC signatures in each sample as a stacked barplot.
    # Select signatures with some contribution
    select = which(rowSums(fit_res$contribution) > 10)
    # Consolidate to matrix with sample name
    fit_res_contrib_mat = as.matrix(fit_res$contribution[select, ])
    # Order by average contribution
    fit_res_contrib_mat = fit_res_contrib_mat[order(rowSums(fit_res_contrib_mat), decreasing = F), ]
    # Rearrange to put signatures of interest at top (POLE, MMR, APOBEC, combos)
    # Which relevant sigs are present in the contrib matrix:
    relevant_sigs_in_contrib = fit_res_contrib_mat[relevant_signatures[relevant_signatures %in% rownames(fit_res_contrib_mat)],]
    # Make vector of names of relvant sigs present in matrix
    relevant_sigs_in_contrib_names = rownames(relevant_sigs_in_contrib)
    # Move choice sigs to top of contrib matrix
    fit_res_contrib_mat = fit_res_contrib_mat[c(relevant_sigs_in_contrib_names, setdiff(rownames(fit_res_contrib_mat), relevant_sigs_in_contrib_names)),] # mat[c(select, setdiff(rownames(mat), select)),]
    # Select contrib mat rownames and reverse order
    rev_order_contrib_rownames = rev(rownames(fit_res_contrib_mat))
    fit_res_contrib_mat_rev = fit_res_contrib_mat[rev_order_contrib_rownames,]
    # Define custom color palette for main sigs present (some will be missing bc not in matrix)
        # POLE sigs (10a, 10b, 28)
        col_pole = brewer.pal(3, 'Blues')
        # Combo sigs (14)
        col_combo = c('#b457b4') # purple
        # MMR sigs (6, 15, 26)
        col_mmr = brewer.pal(3, 'Reds')
        # APOBEC sigs (1, 84, 85)
        col_apobec = c('orange')
        # Other relevant sigs (3, 9, 30, 36)
        col_other_relevant = brewer.pal(4, 'Greens')
        # Colors for other present signatures
        col_others = colorRampPalette(brewer.pal(9, 'Greys'))(25)
        # Total custom color palette of 30 colors
        color_palette_37 = c(col_others, col_other_relevant, col_apobec, col_mmr, col_combo, col_pole)
        
    # Plot contribution barplot
    # Absolute
    pdf(file = paste0(results_plots_dir, "09--COSMICv3-sig-abs-contrib.pdf"), width = 6, height = 8)
    plot_contribution(fit_res_contrib_mat_rev,
                      cosmicv3_sigs_no_arts[ ,select],
                      coord_flip = F, 
                      mode = "absolute", 
                      palette = color_palette_37)
    dev.off()
    # Relative
    pdf(file = paste0(results_plots_dir, "09--COSMICv3-sig-rel-contrib.pdf"), width = 6, height = 8)
    plot_contribution(fit_res_contrib_mat_rev,
                      cosmicv3_sigs_no_arts[ ,select],
                      coord_flip = F,
                      mode = "relative", 
                      palette = color_palette_37)
    dev.off()
 
    # Plot relative contribution of each COSMICv3 signature as a heatmap with sample clustering
    pdf(file = paste0(results_plots_dir, "09--COSMICv3-sig-rel-contrib-heatmap.pdf"), width = 20, height = 8)
    plot_contribution_heatmap(fit_res$contribution, 
                              cluster_samples = T ,
                              method = "complete")
    dev.off()

    # Compare the reconstructed mutational proﬁle of sample 1 with its original mutational proﬁle:
    for (i in 1:length(colnames(mut_mat))) {
      pdf(file = paste0(results_plots_dir, "10--COSMICv3-sigs-reconstructed-", colnames(mut_mat)[i], ".pdf"), width = 12, height = 6)
      print(plot_compare_profiles(mut_mat[,i], 
                                  fit_res$reconstructed[,i], 
                                  profile_names = c('Original','Reconstructed')))
      dev.off()
    }
    
    # Calculate the cosine similarity between all original and reconstructed mutational profiles:
    # Calculate all pairwise cosine similarities
    cos_sim_ori_rec = cos_sim_matrix(mut_mat, fit_res$reconstructed)
    # Extract cosine similarities per sample between original and reconstructed
    cos_sim_ori_rec = as.data.frame(diag(cos_sim_ori_rec))

    # Use ggplot2 to make a barplot of the cosine similarities between the original and reconstructed mutational profile of each sample.
    # Adjust data frame for plotting with gpplot
    colnames(cos_sim_ori_rec) = "cos_sim"
    cos_sim_ori_rec$sample = row.names(cos_sim_ori_rec)
    # Make barplot
    ggplot(cos_sim_ori_rec, aes(y = cos_sim, x = sample)) +
        geom_bar(stat = "identity", fill = "cornflowerblue") +
        coord_cartesian(ylim = c(0.8, 1)) +
        # coord_flip(ylim = c(0.8, 1)) +
        ylab("Sample-COSMIC Cosine Similarity\nOriginal v. Reconstructed") +
        xlab("") +
        # Reverse order of the samples such that first is up
        # xlim(rev(levels(factor(cos_sim_ori_rec$sample)))) +
        theme_classic() +
        theme(panel.grid.minor.y = element_blank(),
              panel.grid.major.y = element_blank()) +
        # Add cut off line
        geom_hline(aes(yintercept = 0.95))
    ggsave(filename = paste0(results_plots_dir, "11--cos-sim-sample-COSMIC-sigs.pdf"), width = 6, height = 4, dpi = "retina")

# 5) Strand bias analyses -------------------------------------------------

    # 5.1) Transcriptional strand bias analysis ====

    # Get gene deﬁnitions for your reference genome:
    # BiocManager::install("TxDb.Mmusculus.UCSC.mm10.knownGene")
    library("TxDb.Mmusculus.UCSC.mm10.knownGene")
    genes_mm10 = genes(TxDb.Mmusculus.UCSC.mm10.knownGene)
    head(genes_mm10)

    # Get transcriptional strand information for all positions in the ﬁrst VCF object
    strand = mut_strand(vcfs[[1]], genes_mm10)
    head(strand)

    # Make mutation count matrix with transcriptional strand information
    mut_mat_s = mut_matrix_stranded(vcfs, ref_genome, genes_mm10)
    head(mut_mat_s)

    # Count the number of mutations on each strand, per tissue, per mutation type:
    strand_counts = strand_occurrences(mut_mat_s, by = tissue)
    head(strand_counts)

    # Perform Poisson test for strand asymmetry signiﬁcance testing:
    strand_bias = strand_bias_test(strand_counts)
    strand_bias

    # Plot the mutation spectrum with strand distinction:
    ps1 = plot_strand(strand_counts, mode = "relative")
    # Plot the eﬀect size (log2(untranscribed/transcribed) of the strand bias. Asteriks indicate signiﬁcant strand bias.
    ps2 = plot_strand_bias(strand_bias)
    # Combine both plots
    pdf(file = paste0(results_plots_dir, "12--strand-biases.pdf"), width = 8, height = 6)
    grid.arrange(ps1, ps2)
    dev.off()

    # 5.2) Replicative strand bias analysis ====
    # Can't robustly do without doing a Repli-seq experiment
    # "replication timing is dynamic and cell-type speciﬁc, which makes replication strand determination less straightforward than transcriptional strand bias analysis. Replication timing proﬁles can be generated with Repli-Seq experiments. Once the replication direction is deﬁned, a strand asymmetry analysis can be performed similarly as the transcription strand bias analysis." - Intro to Mut. Pat.
    
    # 5.3) Extract signatures with strand bias
    
    # TODO: re-run NMF rank estimate on this matrix to find best rank?
    # Extract 2 signatures from mutation count matrix with strand features
    nmf_res_strand = extract_signatures(mut_mat_s, rank = 2)
    # Provide signature names
    colnames(nmf_res_strand$signatures) = c('A','B')
    # Plot signatures with 192 features
    a = plot_192_profile(nmf_res_strand$signatures)
    # Plot strand bias per mutation type for each signature with significance test
    # FIXME: not plotting anything
    b = plot_signature_strand_bias(nmf_res_strand$signatures)
    # Combine plots into one figure
    pdf(file = paste0(results_plots_dir, "13--strand-bias-sigs.pdf"), width = 16, height = 6)
    grid.arrange(a, b, ncol = 2, widths = c(5, 1.8))
    dev.off()
    
# 6) Genomic distribution -------------------------------------------------

    # 6.1) Rainfall plot ====

    # Make rainfall plot of sample 1 over all autosomal chromosomes
    # Define autosomal chromosomes
    chromosomes = seqnames(get(ref_genome))[1:25] # Goes from chr1-chr22,X,Y,M. A number >25 starts to include e.g. "chr1_GL383518v1_alt" chromosomes. The total number of named chromosomes chromosome categories is 455
    chromosomes
    # Make rainfall plot
    for (i in 1:length(colnames(mut_mat))) {
        assign(paste0("rainfall_", colnames(mut_mat)[i]), plot_rainfall(vcfs[[i]], title = names(vcfs[i]), chromosomes = chromosomes, cex = 1.5, ylim = 1e+09))
      pdf(file = paste0(results_plots_dir, "14--genomic-rainfall-", colnames(mut_mat)[i], ".pdf"), width = 12, height = 4)
        print(plot_rainfall(vcfs[[i]], title = names(vcfs[i]), chromosomes = chromosomes, cex = 1.5, ylim = 1e+09))
      dev.off()
    }
    # Combine rainfall plots of all samples
    pdf(file = paste0(results_plots_dir, "14--genomic-rainfall-all.pdf"), width = 24, height = 9)
    grid.arrange(`rainfall_AN-A046`, `rainfall_AX-A05Z`, `rainfall_B5-A0JY`, `rainfall_D1-A17Q`, `rainfall_F5-6814`, `rainfall_IB-7651`, `rainfall_QF-A5YS`)
    dev.off()

      
#     # 6.2) Enrichment or depletion of mutations in genomic regions ====
#     # TODO: Don't have the background to do this yet
#     
#     # 6.3) Test for signiﬁcant depletion or enrichment in genomic regions ====
#     # TODO: Dont have the background to do this yet
#     
#     