########## PCA ##############
### load data
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

### data visual
## pca
library(ggfortify)
inFile = "peaks.table.txt"
dat <- as.data.frame(read.table(inFile, header = T, check.names = F))
X <- as.matrix(dat[,3:dim(dat)[2]])

#plot(X.pca)
autoplot(prcomp(X), data=dat, colour='T1', size=1, 
         loadings = TRUE, 
         loadings.colour = 'blue', 
         loadings.label = TRUE)


