#!/usr/bin/env Rscript

######### checking dependencies, installing libraries. #######
insPacs = installed.packages()[,"Package"]
if(!("statmod" %in% insPacs)) install.packages("statmod", repos = "http://cran.us.r-project.org")
if(!("ggplot2" %in% insPacs)) install.packages("ggplot2", repos = "http://cran.us.r-project.org")
if(!("gridExtra" %in% insPacs)) install.packages("gridExtra", repos = "http://cran.us.r-project.org")
if(!("preprocessCore" %in% insPacs)) install.packages("preprocessCore", repos = "http://cran.us.r-project.org")
if(!("ggfortify" %in% insPacs)) install.packages("ggfortify", repos = "http://cran.us.r-project.org")
if(!("edgeR" %in% insPacs)){
  source("https://bioconductor.org/biocLite.R")
  biocLite("edgeR")}

library(preprocessCore)
library(statmod)
library(edgeR)
library(gridExtra)
library(ggfortify)
library(ggplot2)

################## Input args and input. #####################
args <- commandArgs(TRUE)
if(length(args) < 1) args <- c("--help")
if("--help" %in% args) {
  cat("
	Arguments:
	--inFile=[inFile.tsv, please include a header]
	--fg=[colomn names as foreground/test group, eg. 'exp1,exp2']
  --bg=[colomn names as background/control group, eg. 'exp3,exp4']
  --min_cpm=[minimal read count cutoff for sgRNA to be deemed as express. Use 0 to disable this filter. Default:5]
  --min_cpm_ratio=[minimal ratio of FGs or BGs that have sgRNA read count higher than MINCPM. Use 0 to disable this filter. Default:0.5]
  --qnorm=[1 for quantile normalization of reads within fgs and within bgs. 2 for all experiments (warning: this is a strong hypothesis). Or you can specify by providing the column names, eg. '-q cis1,cis2;ctr1,ctr2;high1,high2' . Default is 0 for no qnrom.]
	--outDir=[name of output dir]
  --prefix=[optional prefix for each file name. eg.'GPF+mCherry-']
  --plotFormat=[pdf or png. default=pdf.]
	--help

	Example:
	call.gRNA.R --inFile='reads.tsv' --fg=exp1,exp2 --bg=exp3,exp4 --qnorm=1 --outDir='sox2.screen'\n")
  q(save="no")
}
parseArgs <- function(x) strsplit(sub("^--", "", x), "=")
argsDF <- as.data.frame(do.call("rbind", parseArgs(args)))
argsL <- as.list(as.character(argsDF$V2))
names(argsL) <- argsDF$V1

## check args and assign defaults
if(is.null(argsL$inFile) || is.null(argsL$bg) || is.null(argsL$fg)) {
  cat("empty input file or fg/bg group! use --help for more info\n")
  q(save="no")
}
if(is.null(argsL$outDir)) argsL$outDir="crispyOut"
if(is.null(argsL$prefix)) argsL$prefix="gRNA"
if(is.null(argsL$plotFormat)) argsL$plotFormat="pdf"
if(is.null(argsL$qnorm)) argsL$qnorm=0
if(is.null(argsL$min_cpm)) argsL$min_cpm=5
if(is.null(argsL$min_cpm_ratio)) argsL$min_cpm_ratio=0.5

inFile = argsL$inFile
fgs = unlist(strsplit(argsL$fg,","))
bgs = unlist(strsplit(argsL$bg,","))
prefix = as.character(argsL$prefix)
qnorm = as.character(argsL$qnorm)
min_cpm = as.numeric(argsL$min_cpm)
min_cpm_ratio = as.numeric(argsL$min_cpm_ratio)
plotFormat = as.character(argsL$plotFormat)
outDir=argsL$outDir
dir.create(file.path(outDir), showWarnings = FALSE)
message = paste0("using infile = ", inFile, "\n",
               "using as fg columns = ", list(fgs), "\n",
               "using as bg columns = ", list(bgs), "\n",
               "output directory to ", outDir, "\n", 
               "prefix = ", prefix, "\n",
               "qnorm = ", qnorm,"\n",
               "min_cpm = ",min_cpm,"\n",
               "min_cpm_ratio = ",min_cpm_ratio,"\n")
cat(message)

####################### debug input start from here. #######################
# #debug
# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
# inFile = "../demos/d1.Yarui/data.tsv"
# fgs=c("cis1")
# bgs=c("high1","high2","high3")
# qnorm = 0

### read data. 
dat = read.table(inFile, header = T)

### derived values and more warnings
nCol = dim(dat)[2]
allExps = c(2:(nCol-1))

## check qnorm parameter. 
# convert qnorm string input to list 
str2lists <- function(inStr){
  a = strsplit(inStr,";")[[1]]
  b = strsplit(a,",")
  return(b)
}
if(grepl(",", qnorm, fixed=TRUE)){
  qnormLists = str2lists(qnorm)
  qnormFlat = unlist(qnormLists)
  if(all(qnormFlat %in% colnames(dat)[allExps]) == FALSE) stop("Specified qnorm experiments not in read count header!")
  for (qnormList in qnormLists){
    qnormArray = unlist(qnormList)
    if(length(qnormArray)<=1) stop("Needs 2 or more specified experiments for each qnorm group!")
  }
}
if(qnorm==2) cat("qnorm=2. Quantile normalization on reads from all experiments (warning: strong hypothesis!) ...\n")

## check number of rep
if(length(fgs)==1 || length(bgs)==1){
  cat("#### Warning: ####\n
      no biological replicates provided, using pooled fg/bg for dispersion estimation\n
      For details, checkout estimateGLMCommonDisp in edgeR.\n\n")
  hasRep = 0
}else{
  hasRep = 1
}

################## preprocessing ##################### 
### quantile normalization of reads
qnormFun <- function(df){
  mat = as.matrix(df)
  df.qnorm = as.data.frame(normalize.quantiles(mat))
  rownames(df.qnorm) = rownames(df)
  colnames(df.qnorm) = colnames(df)
  return(df.qnorm)
  }

## qnorm at given named columns, in place.
qnormCols <- function(df, cols){
  df[,cols] = qnormFun(df[,cols])
  return(df)
  }

## run qnorm 
if(qnorm=="1"){
  cat("qnorm within fgs and within bgs\n")
  if(length(fgs)>1) dat = qnormCols(dat, fgs)
  if(length(bgs)>1) dat = qnormCols(dat, bgs)
}else if(qnorm=="2"){
  cat("qnorm in all experiments\n")
  dat = qnormCols(dat, allExps)
}else if(grepl(",", qnorm, fixed=TRUE)){
  cat("qnorm in the specified columns\n")
  for(qnormList in str2lists(qnorm)){
    qnormArray = unlist(qnormList)
    dat = qnormCols(dat, qnormArray)
  }
}

### boxplot on all reads. ###
cat("Plotting distribution of sgRNA reads ...\n")
X = dat[, allExps]
X_stack = stack(X)
colnames(X_stack) <- c("Read","Exp")

distr1 <- ggplot(X_stack, aes(x=Exp, y=log10(Read+1))) + 
  geom_violin() +
  geom_boxplot(width = 0.1, outlier.shape = NA) +
  xlab("Experiments") +
  ylab("Reads counts (log10(Read+1))") +
  theme_classic()

### PCA on all reads. ###
cat("Running PCA ...\n")
## pca1 # view by sgRNAs
ppca1 <- autoplot(prcomp(X), data=dat, colour = colnames(dat)[dim(dat)[2]],
                  loadings = T,
                  loadings.colour = 'blue', 
                  loadings.label = T) + theme_classic()
## pca2 # view by experiments
tX <- as.data.frame(t(X))
ppca2 <- autoplot(prcomp(tX), 
                  label = TRUE, 
                  shape = FALSE) + theme_classic()

### keep sgRNA that "expressed in at least one conditions" ### following guideline from edgeR
# # the original crest-seq cpm filter
# min_cpm = 3
# min_cpm_ratio = 1/3
# row.use = (rowMeans(data.frame(dat[,c(fgs,bgs)] > min_cpm)) >= min_cpm_ratio)

# the CRISPY cpm filter, scalable and fixes unbalanced issure
# min_cpm = 5
# min_cpm_ratio = 0.5
row.use.fgs = (rowMeans(data.frame(dat[,c(fgs)] > min_cpm)) >= min_cpm_ratio)
row.use.bgs = (rowMeans(data.frame(dat[,c(bgs)] > min_cpm)) >= min_cpm_ratio)
row.use = (row.use.fgs | row.use.bgs)

#
dat.filtered = dat[row.use, ]
nsgRNA = dim(dat)[1]
nsgRNA.filtered = dim(dat.filtered)[1]
msg = paste0("min_cpm: ", round(min_cpm,2), "\n",
             "min_cpm_ratio: ", round(min_cpm_ratio,2), "\n",
             nsgRNA.filtered, "/", nsgRNA," (", round(nsgRNA.filtered/nsgRNA*100,2), "%)")

### cumulative percentile ### only use fg and bg ###
## average
X = data.frame(FGs = rowMeans(data.frame(dat.filtered[,c(fgs)])),
               BGs = rowMeans(data.frame(dat.filtered[,c(bgs)])))
X_stack = stack(X)
colnames(X_stack) <- c("Read","Exp")
distr2 <- ggplot(X_stack, aes(x=log10(Read+1), color=Exp)) +
  stat_ecdf(geom = "line") +
  xlab("Filtered normalized reads counts -- Log10(N+1)") +
  ylab("Cumulative frequency") +
  annotate("text",x=max(log10(X_stack$Read+1)),y=0,label=msg, hjust=1, vjust = 0) +
  theme_classic()

## new way to see distro. 
X = data.frame(FGs = rowMeans(log10(data.frame(dat.filtered[,c(fgs)])+1)),
               BGs = rowMeans(log10(data.frame(dat.filtered[,c(bgs)])+1)),
               Group = dat.filtered[,nCol])
distr3 <- ggplot(X, aes(x=BGs, y=FGs, color=Group)) +
  geom_point(size = 0.5) +
  xlab("Filtered normalized reads counts in BGs -- Log10(N+1)") +
  ylab("Filtered normalized reads counts in FGs") +
  theme_classic()

# ### reads percentile table. ###
# qtab <- round(apply(X,2,quantile, probs = c(0,0.25,0.5, 0.75,1)), 2)
# #libray(ggpubr)
# #distr2_tab <- ggtexttable(qtab, theme = ttheme("classic"))
# outTsv=paste0(outDir, "/", paste0(prefix, ".readsQt.tsv"))
# #outTsv = "test.tsv"
# write.table(qtab, file=outTsv, quote=FALSE, sep='\t')

### negative binomial test ###
cat("Negative binomial test on sgRNAs ...\n")
status = as.factor(dat.filtered[,nCol])
group = c(rep("fg", length(fgs)), rep("bg", length(bgs)))
design = model.matrix(~group)
reads = dat.filtered[,c(fgs,bgs)]
dlist = DGEList(as.matrix(reads))
if(hasRep==1){
  d = estimateDisp(calcNormFactors(dlist), design)
  fit = glmQLFit(d, design, robust=TRUE)
  results = glmQLFTest(fit)
}else{
  d = estimateGLMCommonDisp(dlist, method="deviance", robust=TRUE, subset=NULL)
  fit <- glmFit(d, design)
  results <- glmLRT(fit)
}

tab = cbind(results$table, status)
rownames(tab) = dat.filtered[,1]

## preview and output
cat("preview top results:\n")
topTags(results)
outTsv=paste0(outDir, "/", paste0(prefix, ".sgRNA.tsv"))
write.table(format(tab,digits =4), file=outTsv, quote=FALSE, sep='\t')

## plotting QC
# fc.vs.cpm
p1 <- ggplot(tab, aes(x=logCPM,y=logFC)) +
  geom_point(aes(colour = status), size = 1) + theme_classic()

# # pval.vs.fc
# p2 <- ggplot(tab, aes(x=logFC,y=PValue)) +
#   geom_point(aes(colour = status), size = 1) + theme_classic()

## output all figures to file
outQcPlot.width = 15
outQcPlot.height = 9
png.res = 200 # ppi
if(plotFormat=="png"){
  ## png, fast but lower quality
  outQcPlot = paste0(outDir, "/", paste0(prefix, ".qc.png"))
  cat(paste0("plot QC file = ", outQcPlot,"\n"))
  png(outQcPlot, width = outQcPlot.width, height = outQcPlot.height, units = "in", res=png.res)
}else{
  ## defalut, pdf (high quality but slow rendering)
  outQcPlot = paste0(outDir, "/", paste0(prefix, ".qc.pdf"))
  cat(paste0("plot QC file = ", outQcPlot,"\n"))
  pdf(outQcPlot, width = outQcPlot.width, height = outQcPlot.height)
}
grid.arrange(distr1, ppca1, ppca2, distr2, distr3, p1, nrow = 2)
dev.off()


### TODO: Clean up duplicated X, X_stack, dat, reads situation. 