
#' @title Select Optimal Cutoff on the training set
#'
#' @description This function computes the optimal cutoff PI^{*,T} on T, i.e., the value that corresponds to the best separation
#' in high-and-low risk group with respect to the log-rank test using the training set.
#' @param x1 input training matrix [nxp]
#' @param y1 response variable, y1 should be a two-column matrix with columns named time and status. The latter is a binary variable, with 1 indicating event, and 0 indicating right censored.
#' @param beta regression coefficients selected by the network methods after the variable screening (BMD- or DAD- or BMD+DAD-screening)
#'
#' @return PI.train, cutoff.opt, p.value

SelectOptimalCutoff <- function(x1,y1,beta){

  library(survival)

  PI.train=x1 %*% beta # [nxp]*[px1]=[nx1]

  # Compute the quantile q_gamma, with gamma = 0.2, ..., 0.80
  probs=seq(0.2 ,0.8, by=0.05)
  perc=paste(probs*100, "%", sep="")
  q <- quantile(PI.train, probs=probs)
  p <- vector()
  #par(mfrow=c(3,3))

  for (j in 1:length(q)){
    grouprisk=vector()
    for (i in c(1:nrow(PI.train))){
      if (PI.train[i]>=q[j]) {grouprisk[i]=1}
      else if (PI.train[i]<q[j]) {grouprisk[i]=2}
    }

    #print(table(grouprisk))
    grouprisk <- grouprisk
    x1 <- data.frame(x1)
    logranktest <- survdiff(Surv(y1[,1], y1[,2]) ~ grouprisk, data = x1, rho = 0)
    p.quantile <- 1-pchisq(logranktest$chisq, 1)
    p[j] <- as.vector(signif(p.quantile,3))

    fit <- survfit(Surv(y1[,1], y1[,2]) ~ grouprisk)
    cutoff <- signif(q[j],3)
  }

  q.p.values <- rbind(q,p)
  index.non.zero <- which(q.p.values[2,]!=0)

  if(length(index.non.zero)==1){
    opt.cutoff <- q.p.values[1,index.non.zero]
    opt.cutoff <- signif(opt.cutoff,3)
    # print(paste("Optimal cutoff:",opt.cutoff, sep=" "))
    p.value <- q.p.values[2,index.non.zero]
    p.value <- signif(p.value,3)
  }

  if(length(index.non.zero)!=1 && sum(index.non.zero)!=0){
    p.values <- q.p.values[2,index.non.zero]
    index.p <- which(p.values==min(p.values))
    cutoff.value.min <- min(q.p.values[1,index.p])
    opt.cutoff <- signif(cutoff.value.min,3)
    # print(paste("Optimal cutoff:",opt.cutoff, sep=" "))
    p.value <- q.p.values[2,index.p]
    p.value <- signif(p.value,3)
  }

  if(sum(index.non.zero)==0){
    opt.cutoff <- q.p.values[1,1]
    opt.cutoff <- signif(opt.cutoff,3)
    # print(paste("Optimal cutoff:",opt.cutoff, sep=" "))
    p.value <- q.p.values[2,1]
    p.value <- signif(p.value,3)
  }

  return(list(PI.train=PI.train,opt.cutoff=opt.cutoff,p.value=p.value))
}
