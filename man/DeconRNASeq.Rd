\name{DeconRNASeq}
\alias{DeconRNASeq}
\title{Function for Deconvolution of Complex Samples from RNA-Seq.}
\description{
		This function predicts proportions of constituting cell types from gene expression data generated from RNA-Seq data. 
		Perform nonnegative quadratic programming to get per-sample based globally optimized solutions for constituting cell types .
}
\usage{
DeconRNASeq(datasets, signatures, proportions = NULL, checksig = FALSE, known.prop = FALSE)
}
\arguments{
  \item{datasets}{measured mixture data matrix, genes (transcripts) e.g. gene counts by samples, . The user can choose the appropriate counts, RPKM, FPKM etc.. }
  \item{signatures}{signature matrix from different tissue/cell types, genes (transcripts) by cell types. For gene counts, the user can choose the appropriate counts, RPKM, FPKM etc.. }
  \item{proportions}{proportion matrix from different tissue/cell types. }
  \item{checksig}{whether the condition number of signature matrix should be checked}
  \item{known.prop}{whether the proportions of cell types have been known in advanced for proof of concept}
}
\details{
      Data in the originally measured mixuture sample matrix: datasets and reference matrix: signatures, need to be non-negative.
  	We recommend to deconvolute without log-scale.
}
\value{
      Function DeconRNA-Seq returns a list of results
	\item{out.all}{estimated cell type fraction matrix for all the mixture samples}
      \item{out.pca}{svd calculated PCA on the mixture samples to estimate the number of pure sources according to the cumulative R2}
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
\keyword{methods}
\keyword{DeconSeq}

\examples{
##source("DeconRNASeq.R")
### multi_tissue: expression profiles for 10 mixing samples from multiple tissues 
#data(multi_tissue.rda)  
   
#datasets <- x.data[,2:11]
#signatures <- x.signature.filtered.optimal[,2:6]
#proportions <- fraction

#DeconRNASeq(datasets, signatures, proportions, checksig=FALSE, known.prop = TRUE)
#
}

