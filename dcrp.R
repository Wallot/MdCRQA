dcrp <- function(CRP,lags,ang) {
# This function computes a diagonal cross-recurrence profile
# 
#
# Inputs:
#
#  CRP is a sparse binary (cross-)recurrence matrix
#
#  lags is the number of lags around the diagonal of that matrix for which
#  percent recurrence should be computed
#  
#  ang is the angle with which the cross-recurrence matrix is rotated.
#  If the input matrix is rotated so that the main diagonal (cross-recurrence at lag0) runs
#  from lower-left to upper-right, the matrix needs to be rotated by 90 degrees (i.e., ang = 1).
#  Otherwise, the input matrix does not need to be rotated (i.e., ang = 0).
#  The plots from the mdcrqa-function are rotated by 0 degrees. Hence, ang = 0.
#
#
# Outputs:
#
#  REC = percent recurrence at a respective diagonal around the central
#  diagonal of the recurrence matrix
#
#  lag = lag index for REC
#
#
# Reference:
#
#  Wallot, S. (2017). Multidimensional Cross-Recurrence Quantification
#  Analysis (MdCRQA) - a method for quantifying correlation between
#  multivariate time-series. ???
#
# Version:
#
# v1.0, 04. October 2017
# by Sebastian Wallot, Max Planck Insitute for Empirical Aesthetics, Frankfurt, Germany

  # Load Matrix form Matrix
  library(Matrix)
  
  # check input variables
  if (exists("CRP")) {
  } else {
    print("No input data has been specified")
  }
  
  if (exists("lags")) {
  } else {
    print("lasgs has not been specified")
  }
  
  if (exists("ang")) {
  	if (ang == 1) {
  } else {
    ang <- 0
  }
  }
  
if (ang == 1) {
CRP <- t(apply(CRP, 2, rev))
}
REC <- 0 CRP <- split(CRP, row(CRP) - col(CRP)) for (i in seq(lags+round(length(CRP)/2)+1, -lags+round(length(CRP)/2)+1)) {
	print(i)tempDiagLine <- unlist(CRP[i])REC <- append(REC, sum(tempDiagLine)/length(tempDiagLine))} 
REC <- REC[2:length(REC)]
lag <- seq(-lags, lags)
  
output <- list(REC = REC, lag = lag)
  
  return(output)
}