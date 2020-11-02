library(GenomicRanges)
gr = GRanges(seqnames = Rle(c('chr1','chr2','chr1','chr3'), 
                            c(1, 3, 2, 4)), 
             ranges = IRanges(101:110, end = 111:120, names = head(letters, 10)), 
             strand = Rle(strand(c('-','+','*','+','-')), 
                          c(1,2,2,3,2)), 
             score = 1:10, 
             GC = seq(1,0,length=10))
gr

seqnames(gr)
ranges(gr)
strand(gr)
strand(gr)
granges(gr)
mcols(gr)
mcols(gr)$score

seqlengths(gr) <- c(249250621, 243199373, 198022430) # assigned to each seqname
seqlengths(gr)

names(gr)
length(gr)

sp = split(gr, rep(1:2, each=5))
sp
c(sp[[1]], sp[[2]])

gr
gr[2:3]

gr[2:3, 'GC']

singles = split(gr, names(gr))
grMod = gr
grMod[2] = singles[[1]]
head(grMod, n=3)

