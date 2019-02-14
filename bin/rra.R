#!/usr/bin/env Rscript
################ checking dependencies ###############
insPacs = installed.packages()[,"Package"]
if(!("RobustRankAggreg" %in% insPacs)){
  install.packages("RobustRankAggreg", repos = "http://cran.us.r-project.org")
}

library(RobustRankAggreg)

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
      --inFile=[inFile of regions with comma separated ranked values (0-1). ]
      --outFile=[outFile]
	  --method=[RRA,min,geom.mean,median,stuart], default:RRA
      --help

      Example:
      RRA.R --inFile=sox2.binRank.tsv --outFile=ranked.tsv --method=RRA\n")
  q(save="no")
}

parseArgs <- function(x) strsplit(sub("^--", "", x), "=")
argsDF <- as.data.frame(do.call("rbind", parseArgs(args)))
argsL <- as.list(as.character(argsDF$V2))
names(argsL) <- argsDF$V1
inFile = argsL$inFile
outFile = argsL$outFile
method = argsL$method

if(is.null(argsL$method)) {
  method = "RRA"
  }


############### parameters for debug #############
# library(RobustRankAggreg)
# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
# inFile = "../tmp.4"
# method = "RRA"

################ read and convert data #############
dat = read.table(inFile)
loci = paste(dat[,1], dat[,2], dat[,3], sep="_")
signalLists = dat[,4] ## -log(pval) * sign(log(fg/bg))  ## so the higher is stronger enrichment in fg.  
signalLists = lapply(strsplit(as.character(signalLists), ","), unlist)
df <- plyr::ldply(signalLists, rbind)
mat <- as.matrix(df)
mat2 <- apply(mat, 2, as.numeric)
rownames(mat2) = loci
################ run rra ##############
rraResult <- aggregateRanks(rmat = mat2, method = method, full=T)
CountSgRNA = data.frame(rowSums(!is.na(mat2)))
#results$Count = CountSgRNA
result = merge(rraResult, CountSgRNA, by = 0)[,-1]
colnames(result)=c("Name","Score","Count")
result <- subset(result, Score<1)
result$Score = p.adjust(result$Score, method = "fdr") # change the pval score to FDR; But this reduces resolution of region peaks. 
result <- subset(result, Score<1)
write.table(format(result,digits =4), file=outFile, quote=FALSE, sep='\t', row.names = F)

