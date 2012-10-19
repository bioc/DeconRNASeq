library(DeconRNASeq)
## multi_tissue: expression profiles for 10 mixing samples from multiple tissues 
data(multi_tissue)  
  
datasets <- x.data[,2:11]
signatures <- x.signature.filtered.optimal[,2:6]
proportions <- fraction

DeconRNASeq(datasets, signatures, proportions, checksig=FALSE, known.prop = TRUE, use.scale = TRUE)


##########Example 2: For microarray data, GSE19830#########################

data(rat_liver_brain)

datasets <- all.datasets

proportions <- array.proportions

### From two tissue types, we take the mean to generate the signature matrix
liver_means <- rowMeans(array.signatures[,1:3])
brain_means <- rowMeans(array.signatures[,4:6])
mean_signature <- cbind(liver_means, brain_means)
signatures <- as.data.frame(mean_signature)

## For microarray data, sometimes it is better to use the expression data directly. 
## Here, we set use.scale = F to advoid center and/or scale the columns of our data.
DeconRNASeq(datasets, signatures, proportions, checksig=FALSE, known.prop = TRUE, use.scale = FALSE)
