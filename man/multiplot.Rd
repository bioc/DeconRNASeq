\name{multiplot}
\alias{multiplot}
\title{
     Draw the plots of proportions of cells determined from deconvolution vs. proportions of the cells actually mixed  (when available) with RMSE.
}
\description{
     A function is used to draw the multiple plots of proportions of cells determined from deconvolution vs. proportions of the cells actually mixed. 
     Each plot corresponds to one tissue/cell.  
}
\usage{
multiplot(..., plotlist = NULL, cols)
}
\arguments{
  \item{\dots}{any number of the plot objects that store the scatter plots for all the cells/tissue types }
  \item{plotlist}{any other plot objects}
  \item{cols}{columns of the plots, default = 1}
}
\value{
  A pdf file with the plots of proportions of cells determined from deconvolution vs. proportions of the cells actually mixed with RMSE
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
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (..., plotlist = NULL, cols) 
{
    pdf("scatterplots.pdf")
    require(grid)
    plots <- c(list(...), plotlist)
    numPlots = length(plots)
    plotCols = cols
    plotRows = ceiling(numPlots/plotCols)
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(plotRows, plotCols)))
    vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
    for (i in 1:numPlots) {
        curRow = ceiling(i/plotCols)
        curCol = (i - 1)\%\%plotCols + 1
        print(plots[[i]], vp = vplayout(curRow, curCol))
    }
    dev.off()
  }
}
