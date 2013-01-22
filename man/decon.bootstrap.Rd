\name{decon.bootstrap}
\alias{decon.bootstrap}
\title{
       Estimate the confidence interval for the proportions predicted by deconvolution
}
\description{
       A function is used to estimate the the confidence interval for the proportions predicted by deconvolution through bootstrapping.
}
\usage{
      decon.bootstrap(data.set, possible.signatures, n.sig, n.iter)
}
\arguments{
  \item{data.set}{the data object for mixing samples}
  \item{possible.signatures}{a data frame providing the expression values from pure tissue samples}
  \item{n.sig}{the number of genes/transcripts used for estimation of proportions from our deconvolution}
  \item{n.iter}{the number of bootstraps for our deconvolution}
}
\value{
    A three dimentional array to store means and 95\% confidence interval
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
