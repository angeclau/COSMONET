
#' @title Network Cox regression models
#'
#' @description This function fits penalized Cox regression models with net penalty based on AdaLnet method (Sun, H., et al., 2014).
#' @param x1 input training matrix [n1xp]
#' @param y1 response variable, y1 should be a two-column matrix with columns named time and status. The latter is a binary variable, with 1 indicating event, and 0 indicating right censored.
#' @param penalty penalty type (AdaLnet or ADMMnet)
#' @param alpha ratio between L1 and Laplacian for "Net". The penalty "Net" is defined as
# \lambda*{\alpha*||\beta||_1+(1-\lapha)/2*(\beta^{T}L\beta)}, where L is a Laplacian matrix calculated from Omega.
#' @param ncv number of iterations for the cross validation
#'
#' @return beta, opt.tuningPars

NetworkCox<- function(x1,y1,penalty=c("AdaLnet", "ADMMnet"),alpha=NULL,Omega=NULL,ncv=NULL){
  
  set.seed(1234)
  count = 0;
  
  colnames(y1) <- c("time","status")
  lambdavec = NULL
  bestlambda <- NULL
  for(j in 1:length(alpha)){
    
    for(i in 1:ncv){
      # print(i)
      
      if(penalty=="AdaLnet"){
      cv <- Coxnet::Coxnet(x = x1, y = y1, Omega=Omega, penalty="Net",
                            alpha = alpha[j], nlambda=100, nfolds = 5)
      }
      
      if(penalty=="ADMMnet"){
      cv <- ADMMnet::ADMMnet(x = x1,y = y1,family="cox",penalty="Net",Omega=Omega,
                              alpha = alpha[j], nlambda=100,nfolds = 5)
      }
      
      #lambda.min <- cv$lambda.min
      # print(lambda.min)
      lambdavec <- rbind(lambdavec,cv$lambda.opt)
      
      x <- as.matrix(log(cv$fit[1]))
      y <- as.matrix(cv$fit[2])
      plot(x,y,"p",xlab = "Log(lambda)",ylab ="cvm",col = "red", pch = 20)
      abline(v=log(cv$lambda.opt), col="black",lwd=3, lty=2)
      title(paste("alpha:",alpha[j]," ", "cv:",i, sep=""))
      Sys.sleep(0.1)
    }
    
    lambdaValue <- data.frame(alpha,lambdavec)
    colnames(lambdaValue) <- c("alpha", paste("ncv",i,sep = ""))
    
    lambda.mean <- mean(lambdaValue[j,2:dim(lambdaValue)[2]])
    bestlambda <- rbind(bestlambda, lambda.mean)
    
  }
  
  opt.pars <- data.frame(alpha,bestlambda[,1])
  colnames(opt.pars)[2] <- c("lambda.mean")
  rownames(opt.pars) <- c(1:length(alpha))
  
  ind.lambda.min <- which.min(opt.pars[,2])
  opt.tuningPars <- opt.pars[ind.lambda.min,]
  
  fit <- ADMMnet::ADMMnet(x = x1, y = y1, family="cox", penalty = "Net", Omega = Omega,
                          lambda=opt.tuningPars$lambda.mean,alpha=opt.tuningPars$alpha)
  
  beta = as.matrix(fit$Beta);
  indexBeta = which(beta!=0)
  # length(indexBeta)
  
  if(length(indexBeta)==0){count = count + 1;
  if(count > 5) {print(sprintf("More than 5 beta, file are null"));}
  }
  
  return(list(beta=beta,opt.tuningPars=opt.tuningPars))
  
}
