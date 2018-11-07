#!/usr/bin/env Rscript

######### checking dependencies, installing libraries. #######
insPacs = installed.packages()[,"Package"]
if(!("statmod" %in% insPacs)){
  install.packages("statmod", repos = "http://cran.us.r-project.org")
  }
if(!("ggplot2" %in% insPacs)){
  install.packages("ggplot2", repos = "http://cran.us.r-project.org")
}
if(!("gridExtra" %in% insPacs)){
  install.packages("gridExtra", repos = "http://cran.us.r-project.org")
}
if(!("preprocessCore" %in% insPacs)){
  install.packages("preprocessCore", repos = "http://cran.us.r-project.org")
}
if(!("ggfortify" %in% insPacs)){
  install.packages("ggfortify", repos = "http://cran.us.r-project.org")
}
if(!("edgeR" %in% insPacs)){
  source("https://bioconductor.org/biocLite.R")
  biocLite("edgeR")
}

library(preprocessCore)
library(statmod)
library(edgeR)
library(gridExtra)
library(ggfortify)
library(ggplot2)

################## Input args and input. #####################
args <- commandArgs(TRUE)
if(length(args) < 1) {
  args <- c("--help")
}
if("--help" %in% args) {
  cat("

	Arguments:
	--inFile=[inFile.tsv, please include a header]
	--fg=[colomn names as foreground/test group, eg. 'exp1,exp2']
  --bg=[colomn names as background/control group, eg. 'exp3,exp4']
  --qnorm=[1 for quantile normalization of reads within fgs and within bgs. 2 for all experiments (warning: this is a strong hypothesis). Or you can specify by providing the column names, eg. '-q ctr1,ctr2,high1,high2' . Default is 0 for no qnrom.]
	--outDir=[name of output dir]
  --prefix=[optional prefix for each file name. eg.'GPF+mCherry-']
	--help

	Example:
	call.gRNA.R --inFile='reads.tsv' --fg=exp1,exp2 --bg=exp3,exp4 --qnorm=1 --outDir='sox2.screen'\n")
  q(save="no")
}
parseArgs <- function(x) strsplit(sub("^--", "", x), "=")
argsDF <- as.data.frame(do.call("rbind", parseArgs(args)))
argsL <- as.list(as.character(argsDF$V2))
names(argsL) <- argsDF$V1

## check args and defaults
if(is.null(argsL$inFile) || is.null(argsL$bg) || is.null(argsL$fg)) {
  cat("empty input file or fg/bg group! use --help for more info\n")
  q(save="no")
  }
if(is.null(argsL$outDir)) {
  argsL$outDir="crispyOut"
  }
if(is.null(argsL$prefix)) {
  argsL$prefix="gRNA"
}

inFile=argsL$inFile
fgs = strsplit(argsL$fg,",")
bgs = strsplit(argsL$bg,",")
prefix = as.character(argsL$prefix)
qnorm = as.character(argsL$qnorm)
outDir=argsL$outDir
dir.create(file.path(outDir), showWarnings = FALSE)
message = paste("using infile = ", inFile, "\n",
               "using as fg columns = ", fgs, "\n",
               "using as bg columns = ", bgs, "\n",
               "output directory to ", outDir, "\n", 
               "prefix:", prefix, "\n",
               "qnorm:", qnorm,"\n",
               sep ="")
cat(message)

################## Building gRNA model ##################### 
# #debug
# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
# inFile = "../demos/d1.Yarui/data.tsv"
# fgs=list("cis1", "cis2", "cis3", "cis4", "cis5")
# bgs=list("ctr1", "ctr2")
# qnorm="ctr1,ctr2,high1,high2"
# qnorm=1

# check on rep
hasRep = 1
if(length(unlist(fgs))==1 || length(unlist(bgs))==1){
  hasRep = 0
  cat("#### Warning: ####\nno biological replicates provided, using pooled fg/bg for dispersion estimation\nFor details, checkout estimateGLMCommonDisp in edgeR.\n\n")
}

### input preprocessing
dat = read.table(inFile, header = T)
status = as.factor(dat[,dim(dat)[2]])
ncol = length(c(unlist(fgs), unlist(bgs)))
group = c(rep("fg", length(unlist(fgs))), rep("bg", length(unlist(bgs))))
design = model.matrix(~group)
reads.fgs = dat[,unlist(fgs)]
reads.bgs = dat[,unlist(bgs)]
nSamp = c(2:(dim(dat)[2]-1))

# quantile normalization of reads
qnormFun <- function(df){
  mat = as.matrix(df)
  df.qnorm = normalize.quantiles(mat)
  rownames(df.qnorm) = rownames(df)
  colnames(df.qnorm) = colnames(df)
  return(df.qnorm)
}
# qnorm condition
if(qnorm=="1"){
  if(length(unlist(fgs))>1){
    reads.fgs = qnormFun(reads.fgs)
    cat("using quantile normalization on foreground ...\n")
  }
  if(length(unlist(bgs))>1){
    reads.bgs = qnormFun(reads.bgs)
    cat("using quantile normalization on background ...\n")
  }
  # patch dat with qnormed reads in fg and bg
  reads = cbind(reads.fgs,reads.bgs)
  idxReplace = match(colnames(reads),colnames(dat))
  dat[,idxReplace] <- reads
  X = dat[,nSamp]
}else if(qnorm=="2"){
  cat("using quantile normalization on reads from all experiments (warning: strong hypothesis!) ...\n")
  reads = cbind(reads.fgs,reads.bgs)
  nSamp = c(2:(dim(dat)[2]-1))
  X = dat[,nSamp]
  X = qnormFun(X)
}else if(grepl(",",qnorm,fixed=TRUE)){
  # qnorm the specified columns from input string "qnorm"
  qnormArray = unlist(strsplit(qnorm,","))
  cat(paste0("using specified experiments: ", qnorm, " to perform qnorm ...\n"))
  # check arguments. 
  if(all(qnormArray %in% colnames(dat)[nSamp]) == FALSE){
    stop("Error: Provided qnorm experiments doesn't agree with provided read count header!")
  }
  if(length(qnormArray)<=1){
    stop("Error: too few specified experiments for qnorm!")
  }
  dat.qnormed = qnormFun(dat[,qnormArray])
  dat[,qnormArray] <- dat.qnormed
  reads.fgs = dat[,unlist(fgs)]
  reads.bgs = dat[,unlist(bgs)]
  reads = cbind(reads.fgs,reads.bgs)
  X = dat[,nSamp]
}else{
  cat("no quantile normalization is performed\n")
  reads = cbind(reads.fgs,reads.bgs)
  nSamp = c(2:(dim(dat)[2]-1))
  X = dat[,nSamp]
}

### plot PCA with all reads. 
# pca1
cat("Running PCA on sgRNAs ...\n")
ppca1 <- autoplot(prcomp(X), data=dat, colour = colnames(dat)[dim(dat)[2]],
                  loadings = T,
                  loadings.colour = 'blue', 
                  loadings.label = T) + theme_classic()
# pca2
cat("Running PCA on experiments ...\n")
tX <- as.data.frame(t(dat[,nSamp]))
ppca2 <- autoplot(prcomp(tX), label = TRUE, shape = FALSE) + theme_classic()

### negative binomial test
cat("Negative binomial test on sgRNAs ...\n")
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
rownames(tab) = dat[,1]

## preview and output
cat("preview top results:\n")
topTags(results)
outTsv=paste0(outDir, "/", paste0(prefix, ".sgRNA.tsv"))
write.table(format(tab,digits =4), file=outTsv, quote=FALSE, sep='\t')

## plotting QC
# fc.vs.cpm
p1 <- ggplot(tab, aes(x=logCPM,y=logFC)) +
  geom_point(aes(colour = status), size = 1) + theme_classic()

# pval.vs.fc
p2 <- ggplot(tab, aes(x=logFC,y=PValue)) +
  geom_point(aes(colour = status), size = 1) + theme_classic()

## output all figures to file
## pdf (high quality but slow rendering)
outPdf = paste0(outDir, "/", paste0(prefix, ".qc.pdf"))
cat(paste0("plot QC file = ", outPdf,"\n"))
pdf(outPdf, width = 12, height = 10)
grid.arrange(ppca1, ppca2, p1, p2, nrow = 2)
dev.off()

## png
#outPng = paste0(outDir, "/", paste0(prefix, ".qc.png"))
#cat(paste0("plot QC file = ", outPng,"\n"))
#png(outPng, width = 900, height = 700)
#grid.arrange(ppca1, ppca2, p1, p2, nrow = 2)
#dev.off()
