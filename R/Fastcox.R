
#' @title Fastcox
#'
#' @description This function fits elastic-net penalized Cox's regression in high dimensions models using the cocktail algorithm.
#' @param x1 input training matrix [n1xp]
#' @param y1 response variable, y1 should be a two-column matrix with columns named time and status. The latter is a binary variable, with 1 indicating event, and 0 indicating right censored.
#' @param alpha elasticnet mixing parameter, with 0 < alpha <= 1
#' @param ncv number of iterations for the cross validation
#'
#' @return beta, opt.tuningPars

Fastcox <- function(x1,y1,alpha=NULL,ncv=NULL){

    # print("Fastcox")
    set.seed(1234)
    count = 0

    times = y1[,1]
    status = y1[,2]

    # print(dim(x1))
    # par(mfrow=c(1,ncv))

    lambdavec = NULL
    bestlambda <- NULL
    for(j in 1:length(alpha)){

        for(i in 1:ncv){
            # print(i)
            cv <- fastcox::cv.cocktail(x=x1,y=times,d=status,alpha=alpha[j],nlambda=100,nfolds=5)
            # lambda.min <- cv$lambda.min
            # print(lambda.min)
            lambdavec <- rbind(lambdavec,cv$lambda.min)

            plot(cv)
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

    fit <- fastcox::cocktail(x=x1,y=times,d=status,lambda=opt.tuningPars$lambda.mean,alpha=opt.tuningPars$alpha)

    beta = as.matrix(fit$beta);
    indexBeta = which(beta!=0)
    # length(indexBeta)

    if(length(indexBeta)==0){count = count + 1;
        if(count > 5) {print(sprintf("More than 5 beta, file are null"));}
    }

    return(list(beta=beta,opt.tuningPars=opt.tuningPars))

}
