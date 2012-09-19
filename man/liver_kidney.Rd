\name{liver_kidney}
\docType{data}
\alias{liver_kidney}
\title{data objects for liver and kidney mixing samples}
\description{
 a list containing: 
 
  1) datasets:a data drame providing the RPKM of seven mixing samples.   
  
  2) proportions: a data frame providing the fractions for liver and kidney in the mixing samples

  3) signatures: a data frame providing the expression values from pure liver and kidney samples
}
\usage{liver_kidney}
\format{
 A list
  1) a data frame with 31979 genes' expression on the 7 mixing samples: reads.1, reads.2, reads.3, reads.4, reads.5, reads.6, reads.7
  
  2) a martix whose rows are mixing samples' name and columns are fractions from pure live and kidney tissues
  
  3) a data matrix with 630 expressions from pure liver and kidney tissues
}
\author{Ting Gong \email{tinggong@gmail.com}
        Joseph D. Szustakowski \email{joseph.szustakowski@novartis.com}
}
\keyword{datasets}
\examples{
data(liver_kidney)
}


