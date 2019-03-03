#' @title Survival prediction via screening-network Cox methods
#'
#' @description This function fits penalized Cox regression methods in order to incorporate gene regulatory relationships and to select a subset
#' of potential biomarkers by using the training set. Then, the validation and prediction are performed using the testing set.
#' @param x1 input training matrix [n1xp]
#' @param y1 response variable, y1 should be a two-column matrix with columns named time and status. The latter is a binary variable, with 1 indicating event, and 0 indicating right censored.
#' @param x2 input testing matrix [n2xp]
#' @param y2 response variable, y2 should be a two-column matrix with columns named time and status. The latter is a binary variable, with 1 indicating event, and 0 indicating right censored.
#' @param screenVars screened variables obtained from BMD- or DAD-, or BMD+DAD-screening
#' @param family Cox proportional hazards regression model
#' @param penalty penalty function
#' @param alpha regularization parameter
#' @param Omega network functional matrix
#' @param ncv number of iterations for the cross validation (number of folds by default is 5)
#'
#' @return fitTrain, PI.train, cutoff.opt, p.value
#'
#' @examples Cosmonet('Breast')
#'

Cosmonet <-
  function(x1,y1,x2,y2,screenVars,family="Cox",penalty=c("Fastcox","Coxnet","ADMMnet"),alpha=NULL,Omega=NULL,ncv=NULL){

    fit <- CosmonetTraining(x1,y1,screenVars,family="Cox",penalty=c("AdaLnet","ADMMnet"),alpha,Omega,ncv)
    beta <- as.matrix(fit$fit$beta)
    opt.cutoff <- fit$opt.cutoff

    validation <- CosmonetTesting(x2,y2,screenVars,beta,opt.cutoff)

    class(fit)="Cosmonet"
    return(list(fit=fit, validation=validation))
}

# devtools::load_all()
