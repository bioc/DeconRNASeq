\name{rmse}
\alias{rmse}
\title{
       Calculate the differences between proportions predicted by deconvolution and the values actually measured
}
\description{
       A function is used to calculate the root-mean-square error (RMSE) for the accurracy of estimated proportions.
}
\usage{
      rmse(x, y)
}
\arguments{
  \item{x}{proportions from the actual measurement}
  \item{y}{estimated proportions from our deconvolution}
}
\value{
    A number for RMSE
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
