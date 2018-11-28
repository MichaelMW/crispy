#!/usr/bin/env Rscript

### this script input two tsv files with the same sgRNA entries
### output scatter plot of sgRNA pvalues. 
### eg. input sox2.gfp.tsv and sox2.mcherry.tsv

# set current directory as wd
#setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

### load data
#inFile1 = "../results.test/tiling_gfp_q4.sgRNA.tsv"
#inFile2 = "../results.test/tiling_mcherry_q4.sgRNA.tsv"
inFile1 = "../results.test/tiling_gfp.sgRNA.tsv"
inFile2 = "../results.test/tiling_mcherry.sgRNA.tsv"
label1 = "GFP"
label2 = "mCherry"

### get df
df1 = data.frame(read.table(inFile1)[,c(4,5)])
df2 = data.frame(read.table(inFile2)[,c(4,5)])
df = merge(df1, df2, by=0)

### plot
library(ggplot2)
ggplot(df, aes(x = -log(PValue.x), 
               y = -log(PValue.y), 
               color = status.x)) + 
  geom_point() +
  xlab(label1) +
  ylab(label2) +
  theme_classic()
