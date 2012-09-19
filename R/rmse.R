##########################################################
##RMSE with the knowing proportions
#last update: 05/21/2012
#
#ARGUMENTS:
#x: actual fractions
#y: estimated fractions
#
#VALUE:
#RMSW between actual fractions and estimated fractions
################################################################

rmse <- function(x, y) {
  if (length(x) != length(y)) {
    stop(sprintf("x and y are different lengths : %d, %d",
                 length(x), length(y)))
  }

  residuals <- x-y
  res2      <- residuals^2
  res2.mean <- sum(res2) / length(x)
  rms       <- sqrt(res2.mean)

  return(rms)
  
}
