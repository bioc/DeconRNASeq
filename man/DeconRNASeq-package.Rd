\name{DeconRNASeq-package}
\alias{DeconRNASeq_package}
\alias{DeconRNASeq.package}
\docType{package}
\title{
package DeconRNASeq contains function "DeconRNASeq", implementing the decomposition
of RNA-Seq expression profilings of heterogeneous tissues into cell/tissue type specific expression and cell type concentration
based on cell-type-specific reference measurements
}
\description{
Main function "DeconRNASeq" implements an nonnegative decomposition by quadratic programming as datasets = signature*A,
where "datasets"  are the originally measured data matrix (e.g. genes by samples),
"signature" is the signature matrix (genes by cell types) and "A" the cell type concentration matrix (cell types by samples) 
}
\details{
\tabular{ll}{
Package: \tab DeconRNASeq\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2012-05-25\cr
License: \tab GPL version 2 or later\cr
}
DeconRNASeq(datasets, signature)
}
\author{
Ting Gong \email{tinggong@gmail.com}
Joseph D. Szustakowski \email{joseph.szustakowski@novartis.com} 
}
\references{
Gong, T., et al. (2011) 
Optimal Deconvolution of Transcriptional Profiling Data Using Quadratic Programming with Application to Complex Clinical Blood Samples, 
PLoS One, 6, e27156.
}
\keyword{methods}
\keyword{DeconRNASeq}

