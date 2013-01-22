decon.bootstrap <- function(data.set, possible.signatures, n.sig, n.iter) {

    if (is.null(n.iter)) { n.iter <- 100 }

    if (n.sig >= nrow(possible.signatures)) { stop ("attempt to sample too many transcripts from candidate transcripts") }

    ## create an array to store the outcome from each bootstrap interval
    sim.proportions <- array(dim=c(ncol(data.set), ncol(possible.signatures), n.iter))

    ## run the bootstrap; store each matrix in sim.proportions
    for (i in 1:n.iter){
        signatures.alt       <- possible.signatures[sample(nrow(possible.signatures), n.sig),]
        sim.proportions[,,i] <- DeconRNASeq (data.set, signatures.alt, checksig=FALSE, fig = FALSE)$out.all
    }

    ## create another 3 dimensional array to store means and 95% confidence intervals
    sim.estimates <- array(dim=c(ncol(data.set), ncol(possible.signatures), 3),
                           dimnames=list(colnames(data.set),colnames(possible.signatures), c("mean","lower CI","upper CI")))

    ## loop over the 1st two dimensions of the simulated data
    ## compute mean and confidence intervals across the iterations
    for (i in 1:dim(sim.estimates)[1]) {
        for (j in 1:dim(sim.estimates)[2]) {
            this.t <- t.test(sim.proportions[i,j,])
            sim.estimates[i,j,1] <- this.t$estimate
            sim.estimates[i,j,2] <- this.t$conf.int[1]
            sim.estimates[i,j,3] <- this.t$conf.int[2]
        }
    }

    return(sim.estimates)

}
