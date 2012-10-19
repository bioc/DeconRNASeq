###############################################################
#DeconRNASeq
#last update: 04/05/2012
#
#ARGUMENTS:
#datasets: expression profiling data in matrix format.
#Arrange each data set such that rows are genes/transcripts and columns are mixtures.
#signatures: baseline signatures from pure cell/tissue types. Arrange the same way as the datasets
#proportions: known proportions for pure tissue/cell types
#checksig: boolean variable to check the condition number of the matrix consisting of baseline signatures
#known.prop: boolean variable to check whether the fractions are available or not

#VALUE:
#out.all: mix matrix
################################################################

#library(limSolve)
#library(pcaMethods)
#library(ggplot2)

DeconRNASeq = function(datasets, signatures, proportions=NULL, checksig=FALSE, known.prop = FALSE, use.scale = TRUE){

  if (is.null(datasets)) 
      stop(" Missing the mixture dataset, please provide a tab-delimited text file for mixture samples.")
  if (is.null(signatures)) 
      stop(" Missing the signature dataset, please provide a tab-delimited text file for pure tissue/cell types.")
  if (is.null(proportions)&&known.prop) 
      stop(" Missing the known proprotions, please provide a tab-delimited text file containing known fractions for pure tissue/cell types.")


  ## load data 

  #x.signature <- read.table(signatures, header=T, row.names=1, sep="\t")
  #x.data <- read.table(datasets, sep="\t", header=T, row.names=1)
  x.signature <- signatures
  x.data <- datasets

  #error checks
  if (is.data.frame(x.signature)==F) 
      stop("signature datasets must be a dataframe")
  if (sum(is.na(x.signature))>0)
      stop("signature data cannot have NAs. please exclude or impute missing values.")
  if (is.data.frame(x.data)==F)
      stop("mixture datasets must be a dataframe")
  if (sum(is.na(x.data))>0)
      stop("mixture data cannot have NAs. please exclude or impute missing values.")


  numofg <- nrow(x.signature)
  Numofx <- ncol(x.signature)

  #Underdetermined systems check
  if (numofg < Numofx)
      stop("The number of genes is less than the number of cell types, which means less independent equations than unknowns.")


  #PCA analysis to estimate the number of pure cell types in the mixture
  x.data.temp <- prep(x.data, scale = "none", center = TRUE)
  x.data.pca <- pca(x.data.temp, method = "svd", center = FALSE, nPcs = Numofx)
  #summary of PCA
  out.pca <- summary(x.data.pca)
  Var <- R2cum(x.data.pca)
  numofmix <- order(Var>0.99,decreasing=T)[1]
  
  if (Numofx != numofmix) {
     cat("\n Attention: the number of pure cell types =", Numofx, " defined in the signature matrix;\n" )
     cat("\n PCA results indicate that the number of cell types in the mixtures =", numofmix, "\n" )
  }
  
  if (checksig && numofg >=40){
      #generate the steps to check the condition number of the signature matrix 
      step <- seq(20,numofg, by=20) #every 20 genes
      #step-wise compute the condition number for the subsets of sigature
      sig.cond <- sapply(step, function(x) kappa(scale(x.signature[1:x,]))) 
      condplot(step, sig.cond)
  }
  #x.data.pca <- prcomp(x.data, center=T, scale.=T, retx=T)
  

  ## reorder the rows to match the signature file
  common.signature <- rownames(x.signature) %in% rownames(x.data)
  common.data <- rownames(x.data) %in% rownames(x.signature)
  x.data <- x.data[common.data,]
  x.signature <- x.signature[common.signature,]
  
  x.subdata <- x.data[rownames(x.signature),]

  ## number of cell/tissue types
  Numofx <- ncol(x.signature)

  ## quadratic programming preparation
  if (use.scale) {
     AA <- scale(x.signature)
  }
  else {
     AA <- x.signature
  }

  EE <- rep(1, Numofx)
  FF <- 1
  GG <- diag(nrow=Numofx)
  HH <- rep(0, Numofx)

  out.all <- c()

  for (i in colnames(x.subdata)) {

    BB <- x.subdata[,i]
    if (use.scale) {
       BB <- scale(BB)
    }
  
    out <- lsei(AA, BB, EE, FF, GG, HH)

    out.all <- rbind(out.all, out$X)

  }
  
  if (known.prop)
  {
   x.proportions <- proportions

   x.proportions <- x.proportions[colnames(x.data),]

   parray <- ggplot()
   length(parray) <- ncol(out.all)
   
   for (i in 1:ncol(out.all)){
       A.error <- rmse(x.proportions[,i], out.all[,i])
       x <- out.all[,i]
       y <- x.proportions[,i]
  
       xlabel <- paste('estimated ',colnames(x.proportions)[i])
       ylabel <- paste('actual ', colnames(x.proportions)[i])
       main <- sprintf("Scatter plot of proportions,\n RMSE = %.3f", A.error)
       
       parray[[i]] <- ggplot(data.frame(x,y), aes(x,y)) + geom_point(alpha=.3)+opts(title=main)+geom_abline(intercept=0, slope=1, colour = "red", size = 1)+ xlab(xlabel) + ylab(ylabel)
        
   }
   
   multiplot(plotlist = parray, cols=2)
   
  }
  return(list(out.all = out.all, out.pca = out.pca))

}

# Create Line Chart for condition number

condplot = function(step,cond){
  yrange <- range(cond) 
  xrange <- c(1, length(step))

  #png('condnum.png', width=500, height=500)
  # set up the plot
  plot(xrange, yrange, type="n", xlab="Number of genes in the signature", ylab="Condition number", xaxt="n")
  colors <- rainbow(1)
  linetype <- 1
 
  # add lines
  lines(1:length(step), cond, type="b", lwd=1.5,lty=linetype, col=colors)
  
  # draw an axis on the bottom
  step.ind <- seq(1, length(step), by = 3)
  axis(1, at=step.ind,labels=step[step.ind])

  # add a title and subtitle
  title("Condition number of the signature matrix")
  #dev.off()
  
}


