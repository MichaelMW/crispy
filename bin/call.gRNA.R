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
if(!("edgeR" %in% insPacs)){
  source("https://bioconductor.org/biocLite.R")
  biocLite("edgeR")
}

library(edgeR)

################## Input args and input. #####################
args <- commandArgs(TRUE)
if(length(args) < 1) {
  args <- c("--help")
}
if("--help" %in% args) {
  cat("
	crispr.r using edgeR, statmod and ggplot2 libraries.
  needs multiple biological reps in both treatment and control group to work.
	
	Arguments:
	--inFile=[inFile.tsv, please include a header]
	--fg=[colomn names as foreground/test group, eg. 'exp1,exp2']
  --bg=[colomn names as background/control group, eg. 'exp3,exp4']
	--outDir=[name of output dir]
  --prefix=[optional prefix for each file name. eg.'GPF+mCherry-']
	--help

	Example:
	call.gRNA.R --inFile='reads.tsv' --fg=exp1,exp2 --bg=exp3,exp4 --outDir='sox2.screen'\n")
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
outDir=argsL$outDir
dir.create(file.path(outDir), showWarnings = FALSE)
message = paste("using infile = ", inFile, "\n",
               "using as fg columns = ", fgs, "\n",
               "using as bg columns = ", bgs, "\n",
               "output directory to ", outDir, "\n", 
               "prefix:", prefix, "\n",
               sep ="")
cat(message)

# check on rep
hasRep = 1
if(length(unlist(fgs))==1 || length(unlist(bgs))==1){
  hasRep = 0
  cat("#### Warning: ####\nno biological replicates provided, using pooled fg/bg for dispersion estimation\nFor details, checkout estimateGLMCommonDisp in edgeR.\n\n")
}


################## Building gRNA model ##################### 
dat = read.table(inFile, header = T)
reads = dat[,c(unlist(fgs), unlist(bgs))]
status = as.factor(dat[,dim(dat)[2]])
ncol = dim(reads)[2]
group = c(rep("fg", length(unlist(fgs))), rep("bg", length(unlist(bgs))))
design = model.matrix(~group)
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
outTsv=paste0(outDir, "/", paste0(prefix, ".pvalues.tsv"))
write.table(format(tab,digits =4), file=outTsv, quote=FALSE, sep='\t')

## plotting QC
library(ggplot2)
library(gridExtra)
outPdf = paste0(outDir, "/", paste0(prefix, ".qc.pdf"))
cat(paste0("plot QC file = ", outPdf,"\n"))
pdf(outPdf, width = 7, height = 10)
# fc.vs.cpm
p1 <- ggplot(tab, aes(x=logCPM,y=logFC)) +
  geom_point(aes(colour = status), size = 1)
# pval.vs.fc
p2 <- ggplot(tab, aes(x=logFC,y=PValue)) +
  geom_point(aes(colour = status), size = 1)
grid.arrange(p1, p2, nrow = 2)
dev.off()
