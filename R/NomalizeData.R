
#' @title Quantile Normalization using a specified target distribution vector
#'
#' @description This function normalizes the columns of a matrix based upon a specified normalization distribution.
#' @param x1 input training matrix [n1xp]
#' @param x2 input testing matrix [n2xp]
#'
#' @return output normalized matrices

NormalizeData <- function(x1, x2){
  
  library(preprocessCore)
  # Find target
  target <- normalize.quantiles.determine.target(x1,target.length=NULL,subset=NULL)
  # Normalized matrices
  x1.norm <- normalize.quantiles.use.target(x1,targ,copy=TRUE,subset=NULL)
  x2.norm <- normalize.quantiles.use.target(x2,targ,copy=TRUE,subset=NULL)
  
  output <- list(t(x1.norm),t(x2.norm))
  
  return(output)
}
