\name{condplot}
\alias{condplot}
\title{
      Draw the plot of the condition numbers vs. the number of genes in the signature matrix
}
\description{
      A function is used to draw the plot of the condition number of signature matrices of all sizes, from a handful of genes in one extreme to the whole signature in the other.}
\usage{
      condplot(step, cond)
}
\arguments{
     \item{step}{an array with the number of genes used to calculate the condition numbers of signature matrices, default stepwise = 20}
     \item{cond}{an array with the condition numbers of signature matrices}
}
\value{
     a plot for the condition numbers of signature matrices
}
\references{
      Gong, T., et al. (2011) 
      Optimal Deconvolution of Transcriptional Profiling Data Using Quadratic Programming with Application to Complex Clinical Blood Samples, 
      PLoS One, 6, e27156.
}
\author{
        Ting Gong \email{tinggong@gmail.com}
        Joseph D. Szustakowski \email{joseph.szustakowski@novartis.com}
}
\examples{
library(DeconRNASeq)
####################################################################
## toy data example:


      step <- seq(20,1000, by=20) #every 20 genes
      ## cell type-specific gene expression matrix:
      x.signature <- matrix(rexp(2000),ncol=2)
      sig.cond <- sapply(step, function(x) kappa(scale(x.signature[1:x,]))) 
      function (step, cond)
}
