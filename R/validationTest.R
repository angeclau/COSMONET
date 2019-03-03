
#' @title Validation test
#'
#' @description This function performes the validation test using the testing set
#' @param x2 input testing matrix [n2xp]
#' @param y2 response variable, y2 should be a two-column matrix with columns named time and status. The latter is a binary variable, with 1 indicating event, and 0 indicating right censored.
#' @param beta regression cofficients estimated from variable selection methods on the training set
#' @param opt.cutoff optimal cutoff selected on the training set
#'
#' @return PI.test, p.value, high-and low-risk groups and survival curves
#'
validationTest <- function(x2,y2,beta,opt.cutoff){

  library(survival)

  PI.test=x2 %*% beta

  groupRisk=matrix(0,nrow = dim(x2)[1],ncol = 1)
  gr1 = 0;
  gr2 = 0;
  for (i in c(1:nrow(PI.test))){
    if (PI.test[i] >= opt.cutoff) {groupRisk[i]=1; gr1 = 1;}
    else if (PI.test[i] < opt.cutoff) {groupRisk[i]=2; gr2 = 1;}
    #print(i)
  }

  ind.high <- which(groupRisk == 1)
  # high <-data.frame(high_risk_group=rownames(x2)[ind_high])
  ind.low <- which(groupRisk == 2)
  # low <- data.frame(low_risk_group=rownames(x2)[ind_low])

  group <- data.frame(sample=rownames(PI.test),group=groupRisk)
  group$risk <- ""
  group$risk[ind.high] <- "High"
  group$risk[ind.low] <- "Low"

  if(gr1 == 0 | gr2 == 0){
    p.value <- 1  # no splitting - only one group
  } else {
  x2 <-data.frame(x2)
  logranktest <- survdiff(Surv(y2[,1], y2[,2]) ~ groupRisk, data = x2, rho = 0)
  p.quantile <- 1-pchisq(logranktest$chisq, 1)
  p.value <-signif(p.quantile,3)
  }

  # Kaplan-Meier survival curves test set
  fitTest <- survfit(Surv(y2[,1], y2[,2]) ~ groupRisk, data = x2)

  filename = sprintf("survPlot_%s.pdf",paste(dim(beta)[1]))
  #ggsave(filename, print(survp))
  pdf(filename)

  library(survminer)
  survp <- ggsurvplot(
    fitTest,                   # survfit object with calculated statistics.
    data = x2,                 # data used to fit survival curves.
    risk.table = TRUE,         # show risk table.
    pval = TRUE,               # show p-value of log-rank test.
    conf.int = TRUE,           # show confidence intervals for
                               # point estimaes of survival curves.
    #xlim = c(0,max(y2[,1])),  # present narrower X axis, but not affect
                               # survival estimates.
    #break.time.by = 50,       # break X axis in time intervals by 500.
    ggtheme = theme_minimal(), # customize plot and risk table with a theme.
    risk.table.y.text.col = T, # colour risk table text annotations.
    risk.table.y.text = FALSE  # show bars instead of names in text annotations
                               # in legend of risk table
  )
  print(survp, newpage = FALSE)
  dev.off()

  # setEPS()
  # filename = sprintf("KM_plot_%s.eps", th)
  # postscript(filename, height = 4, width = 6)
  # Sys.sleep(0.1);
  # plot(fitTest, lwd = 2, lty = c(1,1),col = c("blue","red"), xlab = 'Time (months)', ylab= 'Estimated Survival Function')
  # legend(250, 0.3, legend=c('low-risk','high-risk'), lty = c(1,1), col = c("red","blue"), lwd = 2, box.lwd = 0, box.col = "white",bg = "white") #bty='n'
  # title("Kaplan-Meier Curves")
  # text(1,0.1, paste("p.value=",p.value), pos=4)
  # dev.off()

  # print(paste("p.value:", p.value, sep=" "))

  return(list(PI.test=PI.test,p.value=p.value, group=group))
}
