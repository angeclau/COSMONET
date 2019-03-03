
#' @title Model Validation and Prediction
#'
#' @description This function performes the model validation and prediction using the testing set
#' @param x2 input testing matrix [n2xp]
#' @param y2 response variable, y2 should be a two-column matrix with columns named time and status. The latter is a binary variable, with 1 indicating event, and 0 indicating right censored.
#' @param screenVars screened variables obtained from BMD- or DAD-, or BMD+DAD-screening
#' @param beta regression cofficients estimated from variable selection methods on the training set
#' @param opt.cutoff optimal cutoff selected on the training set
#'
#' @return PI.test, p.value, high-and low-risk groups and survival curves

CosmonetTesting <- function(x2,y2,screenVars,beta,opt.cutoff){

  indexScreen <- match(screenVars,colnames(x2))
  beta <- as.matrix(beta)
  testingValue <- validationTest(x2[,indexScreen],y2,beta,opt.cutoff)
  p.value <- testingValue$p.value
  HL.group <- testingValue$group
  PI.test <- testingValue$PI.test

  return(list(PI.test=PI.test,p.value=p.value, HL.group=HL.group))
}
