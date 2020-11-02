# Navigate to directory with '\n'-separated list of VAFs
getwd()
setwd("/Users/leo/Desktop/temp/isec/m33_m22--isec--m33_m17")

# import CSV file as a variable
csv_file = read.csv(file = "0002_AF.csv", header = FALSE, sep = "\n")

# make histogram
hist(csv_file[,1], breaks = seq(0,100,10), main = "m22_T_het -isec- m17_T_het Uniques", xlab = "Allele Frequency")
