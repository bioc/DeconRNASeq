\name{x.signature.filtered.optimal}
\docType{data}
\alias{x.signature.filtered.optimal}
\title{selected signatures from multi-tissues pure samples}
\description{
 A data frame providing the expression values from pure brain, muscle, lung, liver and heart samples used for the example in DeconRNA-Seq.
}
\usage{x.signature.filtered.optimal}
\format{
 A martix whose rows are gene symbols and columns are RPKM expressions from pure tissues: based on the filtered signature matrix, the optimal number of genes have been selected for the deconvolution according to the condition numbers
}
\author{Ting Gong \email{tinggong@gmail.com}
        Joseph D. Szustakowski \email{joseph.szustakowski@novartis.com}
}
\keyword{datasets}
\examples{
data(multi_tissue)
}


