#!/usr/bin/env Rscript

# title: mutational signature analysis - single - COSMIC signatures v3
# author: Leo Williams
# date: 2019-07-02

# About -------------------------------------------------------------------

# This script takes a single VCF file as input and analyzes its mutation spectra with the R package MutationalPatterns

# Usage: ====

# Rscript mut-sigs-single.r [VCF_FILE] [SAMPLE_NAME]

# Setup -------------------------------------------------------------------

# ENABLE FOR GENERAL USE:
# args = commandArgs(trailingOnly = T)

# ENABLE FOR TESTING:
setwd("/Users/leo/Documents/documents-Leo/02--work--career/lab--Pursell/03--projects--Pursell-lab/project_002--DNA-seq-analysis-of-mouse-POLE-exo-mutant-tumor-data/04.02--analysis--LOCAL/expt_06--SNP-SNV-mutation-signature-analysis/single/test") # set testing working dir
args = "" # init test command line rgs
args[1] = "input/m079-tn-WES.snp.Somatic.hc.filter.mm10_multianno.vcf" # test VCF file
args[2] = "m079-test" # test sample name

    # 1) Set working directory ====
    # Run script within sample-specific dir within the experiment dir
    # e.g. "project_02/expt_06/sample_01/mut-sigs.r"
    
    # 2) Load requirements ====
    
    library('BSgenome')
    library('MutationalPatterns')
    library('gridExtra')
    library('NMF')
    library('ggplot2')
    
    # 3) Set variables ====
    
    print(paste0("VCF file: ", args[1]))
    print(paste0("Sample name: ", args[2]))
    
    # Set CLI args to variables
    in_sample_vcf = args[1]
    in_sample_name = args[2]
    
    # 4) Make dirs for results ====
    
    results_dir = paste0("results--", as.character(basename(in_sample_vcf)), format(Sys.time(), '--%Y-%m-%d-%H%M%S'))
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
        
        # Download and load reference genome of interest
        ref_genome = "BSgenome.Mmusculus.UCSC.mm10"
        library(ref_genome, 
                character.only = T)
        
        # 2.2) Load example data ====
        
        # Load the VCF file into a GRanges list
        vcfs = read_vcfs_as_granges(vcf_files = in_sample_vcf, 
                                    sample_names = in_sample_name, 
                                    genome = ref_genome)
        summary(vcfs)
        
        # 2.3) Add reference genome contig lengths to vcfs (originally my VCFs didn't have this. Solved by Francis Blokzijl)
        
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
        context = mut_context(vcfs[[1]], 
                              ref_genome = ref_genome)
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
        type_occurrences = mut_type_occurrences(vcf_list = vcfs, 
                                                ref_genome = ref_genome)
        head(type_occurrences)
        # Save mutation type occurrences as a table
        write.table(x = type_occurrences, 
                    file = paste0(results_tables_dir, 
                                  "01--mutation-type-occurrences.txt"), 
                    sep = "\t", col.names = NA, quote = F)
        
        # 3.2) Mutation spectrum ====
        
        # Plot the mutation spectrum with distinction between C>T at CpG sites and other sites
        p1 = plot_spectrum(type_occurrences = type_occurrences, CT = T)
        
        # Combine multiple plots with the gridExtra package
        pdf(file = paste0(results_plots_dir, "01--mut-spectrum.pdf"), width = 8, height = 4)
        p1
        dev.off()
        
        # 3.3) 96 mutational profile ====
        
        # Make a 96 trinucleodide mutation count matrix:
        mut_mat = mut_matrix(vcf_list = vcfs, ref_genome = ref_genome)
        head(mut_mat)
        write.table(mut_mat, 
                    file = paste0(results_tables_dir, "02--96-trint-mut-matrix.txt"), 
                    sep = '\t', col.names = NA, quote = F)
        
        # Plot the 96 profile of the sample
        pdf(file = paste0(results_plots_dir, "02--96-trint-mut-profile.pdf"),
            width = 6, height = 2)
        plot_96_profile(mut_mat)
        dev.off()
        
    # 4) Mutational signatures ------------------------------------------------
    
        # 4.2) Find optimal contribution of known signatures ====
        
        # 4.2.1) COSMIC mutational signatures ####
        
        # Load the COSMIC Mutational Signatures v3 - single-base substitution (SBS) signatures:
        cancer_signatures = read.csv(file = "/Users/leo/Documents/documents-Leo/02--work--career/lab--Pursell/03--projects--Pursell-lab/project_002--DNA-seq-analysis-of-mouse-POLE-exo-mutant-tumor-data/02--data/02--reference-data/COMIC-mutation-signatures/COSMIC-mutation-signatures-v3--formatted.csv", header = T)
            # # TODO: Incorporate new COSMIC signatures release!
            # # Download mutational signatures from COSMIC website:
            # sp_url = paste("https://cancer.sanger.ac.uk/cancergenome/assets/",
            #                "signatures_probabilities.txt",
            #                sep = "")
            # cancer_signatures = read.table(sp_url, sep = "\t", header = T)
        
        # Match the order of the mutation types to MutationalPatterns standard
        new_order = match(row.names(mut_mat), cancer_signatures$Somatic.Mutation.Type)
        # Reorder cancer signatures dataframe
        cancer_signatures = cancer_signatures[as.vector(new_order), ]
        # Add trinucleotide changes names as row.names
        row.names(cancer_signatures) = cancer_signatures$Somatic.Mutation.Type
        # Keep only 96 contributions of the signatures in matrix
        cancer_signatures = as.matrix(cancer_signatures[ ,4:70])
        
        # Plot mutational profile of some specified COSMIC signatures:
        pdf(file = paste0(results_plots_dir, "relevant-COSMIC-signatures.pdf"), 
            width = 12, height = 8)
        plot_96_profile(cancer_signatures[ ,c(2,6,10,13,15,20,26,28)], 
                        condensed = T)
        dev.off()
        
        # Hierarchically cluster the COSMIC signatures based on their similarity with average linkage:
        hclust_cosmic = cluster_signatures(cancer_signatures, method = "average")
        # Store signatures in new order:
        cosmic_order = colnames(cancer_signatures)[hclust_cosmic$order]
        pdf(file = paste0(results_plots_dir, "COSMIC-signatures-hcluster.pdf"), 
            width = 12, height = 8)
        plot(hclust_cosmic)
        dev.off()
        
        # 4.2.2) Similarity between mutational proﬁles and COSMIC signatures ####
        
        # Calculate cosine similarity of sample mutation matrix with one COSMIC signature. In this case, Signature 10:
        cos_sim_sig_10 = cos_sim(mut_mat[ ,1], cancer_signatures[ ,10])
        cos_sim_sig_10
        
        # Calculate pairwise cosine similarity between the sample mutational proﬁle and COSMIC signatures:
        cos_sim_samples_signatures = cos_sim_matrix(mut_mat, cancer_signatures)
        
        # Plot with simple barplot for a single sample
        # FIXME: This isn't that useful. Not every signature is even extracted from the sample.
        pdf(file = paste0(results_plots_dir, "03--cos-sim-sample-COSMIC-sigs.pdf"), width = 8, height = 12)
        par(mar = c(5, 6, 4, 1)+.1)
        barplot(cos_sim_samples_signatures, 
                main = "Cosine Similarity with COSMIC Signatures", 
                horiz = T, xlim = c(0,1),
                names.arg = colnames(cancer_signatures), 
                las=1)
        dev.off()
        
        # 4.2.3) Find optimal contribution of COSMIC signatures to reconstruct 96 mutational proﬁles ####
        
        # Fit mutation matrix to the COSMIC mutational signatures:
        fit_res = fit_to_signatures(mut_mat, cancer_signatures)
        
        # Plot the optimal contribution of the COSMIC signatures in each sample as a stacked barplot.
        # Select signatures with some contribution
        select = which(rowSums(fit_res$contribution) > 10)
        # Consolidate to matrix with sample name
        fit_res_contrib_mat = as.matrix(fit_res$contribution[select, ])
        colnames(fit_res_contrib_mat) = in_sample_name
        # Plot contribution barplot
        pdf(file = paste0(results_plots_dir, "04--COSMIC-sig-abs-contrib.pdf"), width = 4, height = 8)
        plot_contribution(fit_res_contrib_mat,
                          cancer_signatures[ ,select],
                          coord_flip = F,
                          mode = "absolute")
        dev.off()
        
        # Compare the reconstructed mutational proﬁle of sample 1 with its original mutational proﬁle:
        pdf(file = paste0(results_plots_dir, "05--COSMIC-reconstruct-v-orig.pdf"), width = 12, height = 8)
        plot_compare_profiles(mut_mat[ ,1], 
                              fit_res$reconstructed[ ,1], 
                              profile_names = c("Original", "Reconstructed"))
        dev.off()
        
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
            geom_bar(stat = "identity", fill = "skyblue4") + 
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
        ggsave(filename = paste0(results_plots_dir, "06--cos-sim-sample-COSMIC-sigs.pdf"), width = 2, height = 4, dpi = "retina")
         
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
        mut_mat_s[1:5, ]
        
        # Count the number of mutations on each strand, per tissue, per mutation type:
        strand_counts = strand_occurrences(mut_mat_s)
        # TODO: MULTI: add `by = tissue`
        head(strand_counts)
        
        # Perform Poisson test for strand asymmetry signiﬁcance testing:
        strand_bias = strand_bias_test(strand_counts)
        strand_bias
        
        # Plot the mutation spectrum with strand distinction:
        ps1 = plot_strand(strand_counts, mode = "relative")
        # Plot the eﬀect size (log2(untranscribed/transcribed) of the strand bias. Asteriks indicate signiﬁcant strand bias.
        ps2 = plot_strand_bias(strand_bias)
        # Combine both plots
        pdf(file = paste0(results_plots_dir, "07--strand-biases.pdf"), width = 4, height = 6)
        grid.arrange(ps1, ps2)
        dev.off()
        
        # 5.2) Replicative strand bias analysis ====
        # TODO: Need "replicationDirectionRegions.bed" file?
        
    # 6) Genomic distribution -------------------------------------------------
    
        # 6.1) Rainfall plot ====
        
        # Make rainfall plot of sample 1 over all autosomal chromosomes
        # Define autosomal chromosomes
        chromosomes = seqnames(get(ref_genome))[1:22] # Goes from chr1-chrM. A number >22 starts to include e.g. "chr1_GL456210_random" chromosomes. The total number of named chromosomes chromosome categories is 66.
        chromosomes
        # Make rainfall plot
        # FIXME: Chromosome delimiters not showing up
        pdf(file = paste0(results_plots_dir, "08--genomic-dist-rainfall.pdf"), width = 12, height = 4)
        plot_rainfall(vcfs[[1]], title = names(vcfs[1]), chromosomes = chromosomes, cex = 1.5, ylim = 1e+09)
        dev.off()
        
        # 6.2) Enrichment or depletion of mutations in genomic regions ====
        # TODO: Don't have the background to do this yet
        
        # 6.3) Test for signiﬁcant depletion or enrichment in genomic regions ====
        # TODO: Dont have the background to do this yet

