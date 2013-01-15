library(DeconRNASeq)
## multi_tissue: expression profiles for 10 mixing samples from multiple tissues 
data(multi_tissue)  
  
datasets <- x.data[,2:11]
signatures <- x.signature.filtered.optimal[,2:6]
proportions <- fraction

DeconRNASeq(datasets, signatures, proportions, checksig=FALSE, known.prop = TRUE, use.scale = TRUE, fig = TRUE)

#if the true proportions are known, we can estimate the confidence interval for the fraction estimation.
# We used the bootstrapping 100s by selecting the same number of differentially expressed genes as for the 
#optimal signature matrix (i.e. 1570) and recorded the deconvolution results. 
# We computed the mean proportions for all the tissues we estimated and calculated the confidence interval of the estimated proportions using a t-test.


n.iter <- 100

# create a 3 dimensional array to store the simulation outcomes
sim.proportions <- array(dim=c(ncol(datasets), ncol(signatures), n.iter))

# run the bootstrap; store each matrix in sim.proportions
for (i in 1:n.iter){
    tmp <- x.signature.filtered[x.signature.filtered$p_val.filtered.diff==0,]
    signatures.alt <- tmp[sample(1:dim(tmp)[[1]],dim(signatures)[[1]]),2:6]
    sim.proportions[,,i] <- DeconRNASeq (datasets, signatures.alt, proportions, checksig=FALSE, known.prop = TRUE, fig = FALSE)$out.all
}

# create another 3 dimensional array to store means and 95% confidence intervals
sim.estimates <- array(dim=c(ncol(datasets), ncol(signatures), 3), dimnames=list(rownames(proportions),colnames(proportions), c("mean","left CI","rigth CI")))

# loop over the 1st two dimensions of the simulated data
# compute mean and confidence intervals across the iterations
for (i in 1:ncol(datasets)) {
    for (j in 1:ncol(signatures)) {
        this.t <- t.test(sim.proportions[i,j,])
        sim.estimates[i,j,1] <- this.t$estimate
        sim.estimates[i,j,2] <- this.t$conf.int[1]
        sim.estimates[i,j,3] <- this.t$conf.int[2]
    }
}


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
