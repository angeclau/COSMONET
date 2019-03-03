
#' @title Fit network Cox regression penalized methods
#'
#' @description This function fits penalized Cox regression methods in order to incorporate gene regulatory relationships and to select a subset
#' of potential biomarkers by using the training set
#' @param x1 input training matrix [n1xp]
#' @param y1 response variable, y1 should be a two-column matrix with columns named time and status. The latter is a binary variable, with 1 indicating event, and 0 indicating right censored.
#' @param screenVars screened variables obtained from BMD- or DAD-, or BMD+DAD-screening
#' @param family Cox proportional hazards regression model
#' @param penalty penalty function
#' @param alpha regularization parameter
#' @param Omega network functional matrix
#' @param ncv number of iterations for the cross validation (number of folds by default is 5)
#'
#' @return fitTrain, PI.train, cutoff.opt, p.value

CosmonetTraining <- function(x1,y1,screenVars,family="Cox",penalty=c("AdaLnet","ADMMnet"),alpha=NULL,Omega=NULL,ncv=NULL){

  ## Fit network Cox methods on I, i.e., the set of screened variables
  indexScreen <- match(screenVars,colnames(x1))
  penalty <- match.arg(penalty)

    if(family=="Cox") {
      fitTrain <- switch(penalty,
                         "AdaLnet"=NetworkCox(x1[,indexScreen],y1,penalty=c("AdaLnet"),alpha,Omega[indexScreen,indexScreen],ncv),
                         "ADMMnet"=NetworkCox(x1[,indexScreen],y1,penalty=c("ADMMnet"),alpha,Omega[indexScreen,indexScreen],ncv))

      fitTrain$family <- "cox"
    }

  ## Compute the optimal cutoff on $T$
  select.cutoff <- SelectOptimalCutoff(x1[,indexScreen],y1,fitTrain$beta)
  opt.cutoff <- select.cutoff$opt.cutoff
  p.value <- select.cutoff$p.value
  PI.train <- select.cutoff$PI.train

  return(list(fitTrain=fitTrain,PI.train=PI.train,opt.cutoff=opt.cutoff,p.value=p.value))
}
