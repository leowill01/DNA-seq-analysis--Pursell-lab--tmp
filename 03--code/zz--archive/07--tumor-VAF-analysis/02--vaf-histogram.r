#!/usr/local/bin/Rscript

# Make a histogram of tumor VAF from a list of the VAF percentages (without the "%")

# Usage: Rscript <vaf-hist.r> <VAF-list>.txt <sample-name>

args = commandArgs(trailingOnly = T)

# import the VAF list as a vector of numbers
vafVec = read.csv(args[1])
head(vafVec)
str(vafVec)

pdf(file = paste(args[2],"Tumor-VAF.pdf", sep = '-'))
hist(vafVec[[1]], 
     xlim = c(0,100), 
     main = paste(args[2], "Tumor VAF", sep = ' ') )
dev.off()
