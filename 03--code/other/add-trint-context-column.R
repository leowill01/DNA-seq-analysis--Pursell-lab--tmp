library(tibble)

sigv2 <- read.csv(file = "COSMIC-signatures-v2.csv", header = T)
sigv3 <- read.csv(file ="sigProfiler_SBS_signatures_2019_05_22.csv", header = T)

sigv2 = sigv2[order(sigv2$Trinucleotide),]
sigv3 = sigv3[order(sigv3$SubType),]

trint = as.data.frame(sigv2$Somatic.Mutation.Type)

sigv3 <- add_column(.data = sigv3, trint, .after = "SubType")

write.csv(x = sigv3, file = "COSMIC-signatures-v3-SBS-formatted.csv", row.names = F, quote = F)
