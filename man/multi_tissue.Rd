\name{multi_tissue}
\docType{data}
\alias{multi_tissue}
\title{data objects for multi-tissues mixing samples}
\description{
 a list containing: 
 
  1) x.data:a data frame providing the RPKM of nine mixing samples.   
  
  2) x.signatures: a data frame providing the expression values from pure brain, muscle, lung, liver and heart samples.
  
  3) x.signatures.filtered: a data frame providing the expression values from pure brain, muscle, lung, liver and heart samples after filtering.
  
  4) x.signatures.filtered.optimal: a data frame providing the expression values from pure brain, muscle, lung, liver and heart samples used for the example in DeconRNA-Seq.
  
  5)fraction: a data frame providing the fractions from 5 tissues in the mixing samples
}
\usage{multi_tissue}
\format{
 A list
  1) a matrix with all the genes' expression in the mixing samples: the first two columns are corresponding to the RefSeq accession numbers and gene symbols
  
  2) a martix whose rows are gene symbols and columns are RPKM expressions from pure tissues.
  
  3) a martix whose rows are gene symbols and columns are RPKM expressions from pure tissues: the genes with RPKM less than 200 within any of the five tissues have been filtered. 
  
  4) a martix whose rows are gene symbols and columns are RPKM expressions from pure tissues: based on the filtered signature matrix, the optimal number of genes have been selected for the deconvolution according to the condition numbers

  5) a martix whose rows are mixing samples' name and columns are fractions from pure tissues including brain, muscle, lung, liver and heart
}
\author{Ting Gong \email{tinggong@gmail.com}
        Joseph D. Szustakowski \email{joseph.szustakowski@novartis.com}
}
\keyword{datasets}
\examples{
data(multi_tissue)
}


