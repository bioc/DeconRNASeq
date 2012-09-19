library(DeconRNASeq)
## multi_tissue: expression profiles for 10 mixing samples from multiple tissues 
data(multi_tissue)  
  
datasets <- x.data[,2:11]
signatures <- x.signature.filtered.optimal[,2:6]
proportions <- fraction

DeconRNASeq(datasets, signatures, proportions, checksig=FALSE, known.prop = TRUE)