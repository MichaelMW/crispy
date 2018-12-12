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
  --qnorm=[0 for no quantile normalization and 1 for separate qnorm within fgs and within bgs. Or you can specify by providing the column names, eg. '-q cis1,cis2;ctr1,ctr2;high1,high2' . Default:0.]
  --pvalCut =[Use this cutoff to highlight positive sgRNAs passing the filter. No highlight if FALSE. Default:F]
  --direction = [1 or -1. If pvalCut exists, use this cutoff to highlight positive sgRNAs passing the filter. No highlight if FALSE. Default:1 ]
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
if(is.null(argsL$pvalCut)) argsL$pvalCut=FALSE
if(is.null(argsL$direction)) argsL$direction=1

inFile = argsL$inFile
fgs = unlist(strsplit(argsL$fg,","))
bgs = unlist(strsplit(argsL$bg,","))
prefix = as.character(argsL$prefix)
qnorm = as.character(argsL$qnorm)
min_cpm = as.numeric(argsL$min_cpm)
min_cpm_ratio = as.numeric(argsL$min_cpm_ratio)
pvalCut = as.numeric(argsL$pvalCut)
postPvalMode = as.character(argsL$postPvalMode)
direction = as.numeric(argsL$direction)
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
#debug1
# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
# inFile = "../demos/d1.Yarui/data.tsv"
# fgs=c("cis1","cis2","cis3","cis4","cis5")
# bgs=c("ctr1","ctr2")
# qnorm = 0
# min_cpm = 5
# min_cpm_ratio = 0.5
# 
#debug2
# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
# inFile = "../demos/d2.SOX2/reads.tsv"
# bgs=c("S1.Unsorted","S2.Unsorted")
# fgs=c("S1.GpMn","S2.GpMn")
# qnorm = 0
# min_cpm = 5
# min_cpm_ratio = 0.5
# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

### hardcoded parameter
negLabel = "negative"
testLabel = "test"
qcutoffs = c(0, 0.001, 0.005, 0.01, 0.05, 0.1) # pvalue guidelines help to aim at these qcutoffs.

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
  mat = data.matrix(df)
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
}else if(grepl(",", qnorm, fixed=TRUE)){
  cat("qnorm in the specified columns\n")
  for(qnormList in str2lists(qnorm)){
    qnormArray = unlist(qnormList)
    dat = qnormCols(dat, qnormArray)
  }
}

### boxplot on all reads. ###
cat("Plotting distribution of sgRNA reads ...\n")
X = data.frame(dat[, allExps])
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
                  loadings.label = T,
                  size = 0.5) + theme_classic()
## pca2 # view by experiments
tX <- as.data.frame(t(X))
ppca2 <- autoplot(prcomp(tX), 
                  label = TRUE, 
                  loadings.label.repel=T,
                  shape = FALSE) + theme_classic()

### keep sgRNA that "expressed in at least one conditions" ### following guideline from edgeR
## the CRISPY cpm filter, scalable and fixes unbalanced issure
row.use.fgs = (rowMeans(data.frame(dat[,c(fgs)] > min_cpm)) >= min_cpm_ratio)
row.use.bgs = (rowMeans(data.frame(dat[,c(bgs)] > min_cpm)) >= min_cpm_ratio)
row.use = (row.use.fgs | row.use.bgs)

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
  geom_point(alpha = 0.3, size = 0.5) +
  geom_point(size = 0.5, data = subset(X, Group!=testLabel)) +
  xlab("Filtered normalized reads counts in BGs -- Log10(N+1)") +
  ylab("Filtered normalized reads counts in FGs") +
  theme_classic()

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
# highlight postive sgRNA with pvalue < cutoff. 
if(pvalCut > 0){
  tab$status <- as.character(tab$status)
  # default. only positive FC
  if(direction == 1){
    filteredIdx = (tab$PValue < pvalCut & tab$status==testLabel & tab$logFC>0)
    posSgrnaMsg = paste0("Enriched sgRNA p < ", pvalCut)
  }
  # only negative FC
  if(direction == -1){
    filteredIdx = (tab$PValue < pvalCut & tab$status==testLabel & tab$logFC<0)
    posSgrnaMsg = paste0("Depleted sgRNA p < ", pvalCut)
  }
  tab[filteredIdx, "status"] = posSgrnaMsg
}

# fc.vs.cpm. The MA plot
p1 <- ggplot(tab, aes(x=logCPM,y=logFC, colour = status)) +
  geom_point(size = 0.5) + 
  geom_point(size = 0.5, data = subset(tab, status!=testLabel)) +
  theme_classic()

# pval.vs.fc
p2 <- ggplot(tab, aes(x=logFC,y=PValue, colour = status)) +
  geom_point(size = 0.5) +
  geom_point(size = 0.5, data = subset(tab, status!=testLabel)) +
  theme_classic()

## plot pval distribution to guide pvalue cutoff choice. 
tab_pvals = tab[,c("PValue","logFC","status")]
pvalDistro <- ggplot(tab_pvals, aes(x=status, y=-log10(PValue))) + 
  geom_violin() +
  geom_boxplot(width = 0.1) +
  xlab("sgRNA group") +
  ylab("-log10(pval)") +
  scale_y_continuous(breaks=seq(0,max(-log10(tab_pvals$PValue)),1)) +
  theme(panel.grid.major.y = element_line(colour = "black", linetype = "dashed"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.background = element_rect(fill = "white",colour = "black"))

## add pvalue cutoffs guided by FDR in qcutoffs ##
## convert qvalue cutoff to the highest pvalue that satifies.
qval2pval <- function(qcutoffs, pvals, status, testLabel, negLabel){
  df.pval = data.frame(pvals= pvals,status = status)
  df.pval = subset(df.pval, status %in% c(testLabel, negLabel))
  df.pval$status = as.factor(as.character(df.pval$status))
  df.pval.sorted = df.pval[order(df.pval$pvals),]
  N_recNeg = cumsum(df.pval.sorted$status == negLabel)
  N_allNeg = sum(status==negLabel)
  dfm = data.frame(df.pval.sorted, 
                   FDR = N_recNeg/N_allNeg)
  idx = 1:nrow(dfm)
  df.q2p = data.frame()
  for(qcutoff in qcutoffs){
    pcutoff = dfm[max(idx[dfm$FDR <= qcutoff]),]$pval
    pcutoff = round(pcutoff,5)
    df.q2p = rbind(df.q2p, data.frame(qcutoff = qcutoff,
                                      pcutoff = pcutoff))
  }
  return(df.q2p)
}
msg = qval2pval(qcutoffs, tab$PValue, tab$status, testLabel, negLabel)
q2pGuideTab = tableGrob(msg, rows = NULL)  

## output all figures to file
outQcPlot.width = 15
outQcPlot.height = 12
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
grid.arrange(distr1, ppca1, ppca2, distr2, distr3, p1, p2, pvalDistro, q2pGuideTab, nrow = 3)
dev.off()
