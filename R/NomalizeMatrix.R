
#' @title Quantile Normalization
#'
#' @description This function normalizes the columns of a matrix based upon a specified normalization distribution
#' @param x input training matrix [n1xp]
#'
#' @return output normalized matrix

NormalizeMatrix<- function(x){
  
  library(preprocessCore)
  # Find target
  target <- normalize.quantiles.determine.target(x,target.length=NULL,subset=NULL)
  # Normalized matrices
  output <- normalize.quantiles.use.target(x,target,copy=TRUE,subset=NULL)
  
  return(output)
}
