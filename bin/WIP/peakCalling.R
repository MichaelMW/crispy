## replicate peak calling from Rongxin


insPacs = installed.packages()[,"Package"]
if(!("GenomicRanges" %in% insPacs)){
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("GenomicRanges", version = "3.8")}

library(GenomicRanges)


## examples 
gr <- GRanges(
  seqnames=Rle(paste("chr", c(1, 2, 1, 3), sep=""), c(1, 3, 2, 4)),
  ranges=IRanges(1:10, width=10:1, names=letters[1:10]),
  strand=Rle(strand(c("-", "+", "*", "+", "-")), c(1, 2, 2, 3, 2)),
  score=1:10,
  GC=seq(1, 0, length=10)
)
gr

gr1 <- GRanges(seqnames="chr2", ranges=IRanges(3, 6),
               strand="+", score=5L, GC=0.45)
gr2 <- GRanges(seqnames="chr1",
               ranges=IRanges(c(10, 7, 19), width=5),
               strand=c("+", "-", "+"), score=3:5, GC=c(0.3, 0.5, 0.66))
gr3 <- GRanges(seqnames=c("chr1", "chr2"),
               ranges=IRanges(c(1, 4), c(3, 9)),
               strand=c("-", "-"), score=c(6L, 2L), GC=c(0.4, 0.1))
grl <- GRangesList(gr1=gr1, gr2=gr2, gr3=gr3)
grl


reduce(gr)
gr2 <- reduce(gr, with.revmap=TRUE)


