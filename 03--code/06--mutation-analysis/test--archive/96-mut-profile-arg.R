#!/usr/local/bin/Rscript

args = commandArgs(trailingOnly = T)
# args[1] = vcf filepath
# args[2] = sample name

setwd("~/Desktop/test-work/")

library("BSgenome")
library("MutationalPatterns")
library("gridExtra")

# Download and load reference genome (mm10)
mm10RefGenome = "BSgenome.Mmusculus.UCSC.mm10"
library(mm10RefGenome, character.only = T)

vcfFiles = args[1]

sampleName = args[2]

# Load VCF files into a GRangesList
vcfs = read_vcfs_as_granges(vcf_files = vcfFiles, sample_names = sampleName, mm10RefGenome)

mutMat = mut_matrix(vcf_list = vcfs, ref_genome = mm10RefGenome)
head(mutMat)

pdf(file = paste(args[1], "96-mut-profile.pdf", sep = '-'), width = 8, height = 2)
plot_96_profile(mutMat, condensed = T)
dev.off()