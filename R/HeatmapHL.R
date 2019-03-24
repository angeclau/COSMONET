#' @title Heatmap for high- and low-risk group
#'
#' @description This function calculate the heatmap of genes selected by the combination of BMD-screening and network methods by dividing the patiens in high- and low-risk group.
#' @param x1 input training matrix [n1xp]
#' @param x2 input testing matrix [n2xp]
#' @param beta regression coefficients
#' @param PI.train prognostic index on training set
#' @param PI.test prognostic index on testing set
#' @param th threshold for the heatmap (z-score)
#'
#' @return heatmap on traing and  testing set
#'

heatmapHL <- function(x1,x2,beta,PI.train,opt.cutoff,PI.test,th=3.5){

  PI.train.ordered <- as.matrix(PI.train[order(PI.train, decreasing = TRUE),])
  x1.ordered <- x1[order(PI.train, decreasing = TRUE),]
  # dim(x1.ordered)

  # Divide patients in high and low risk group in the training set
  group.risk.train=vector()
  gr1 = 0;
  gr2 = 0;
  for (i in c(1:nrow(PI.train.ordered))){
    if (PI.train.ordered[i] >= PI_star) {group.risk.train[i]="High-risk"; gr1 = 1;} # 1=high-low
    else if (PI.train.ordered[i] < PI_star) {group.risk.train[i]="Low-risk"; gr2 = 1;} # 2=low-risk
    #print(i)
  }

  # beta non zero
  ind.non.zero <- which(beta!=0)

  data.train <- data.frame(colnames(x1.ordered[,ind.non.zero]),t(x1.ordered[,ind_non_zero]))
  colnames(data.train) <- c("genes", rownames(x1.ordered[,ind.non.zero]))
  rownames(data.train) <- c(1:dim(beta[ind.non.zero,])[1])
  # dim(data.train)

  HL.group.train <- data.frame(group.risk.rain)
  rownames(HL.group.train) <- rownames(PI.train.ordered)
  colnames(HL.group.train) <- "Risk-groups"

  s.train <- length(which(group.risk.train == "High-risk"))

  zscore.train <- t(scale(t(as.matrix(data.train[,2:dim(data.train)[2]]))))
  zscore_train[zscore.train < -th] <- -th
  zscore_train[zscore.train > th] <- th
  rownames(zscore.train) <- beta[ind.non.zero,1]

  library(pheatmap)

  colors = colorRampPalette(c("green", "black", "red"))(100)

  heatmap.train <- pheatmap(zscore.train, color = colors, cluster_row = TRUE, cluster_cols = FALSE,
                            gaps_col = s.train, annotation_col = HL.group.train, cutree_col = 2,
                            #cutree_rows = 6,
                            show_colnames = F, fontsize = 6.5,fontsize_row=6, filename = "zscore.train.pdf")

  clust_genes <- heatmap.train$tree_row[["order"]]
  coeff.non.zero <- beta[ind.non.zero,]
  genes.test <- as.character(coeff.non.zero[clust_genes,1])
  order.clust <- cutree(heatmap.train$tree_row,6)[heatmap.train$tree_row[["order"]]]

  # Testing set
  PI.test.ordered <- as.matrix(PI.test[order(PI.test, decreasing = TRUE),])
  x2.ordered <- x2[order(PI.test, decreasing = TRUE),]

  # Divide patients in high and low risk group in the training set
  group.risk.test=vector()
  gr1 = 0;
  gr2 = 0;
  for (i in c(1:nrow(PI.test.ordered))){
    if (PI.test.ordered[i] >= opt.cutoff) {group.risk.test[i]="High-risk"; gr1 = 1;}
    else if (PI.test.ordered[i] < opt.cutoff) {group.risk.test[i]="Low-risk"; gr2 = 1;}
    #print(i)
  }

  data.test <- data.frame(colnames(x2.ordered[,ind.non.zero]),t(x2.ordered[,ind.non.zero]))
  colnames(data.test) <- c("genes", rownames(x2.ordered[,ind.non.zero]))
  rownames(data.test) <- c(1:dim(beta[ind.non.zero,])[1])

  HL.group.test <- data.frame(group.risk.test)
  rownames(HL.group.test) <- rownames(PI.test.ordered)
  colnames(HL.group.test) <- "Risk-groups"

  s.test <- length(which(group.risk.test == "High-risk"))

  zscore.test <- t(scale(t(as.matrix(data_test[clust_genes,2:dim(data_test)[2]]))))
  zscore.test[zscore.test < -th] <- -th
  zscore.test[zscore.test > th] <- th
  rownames(zscore.test) <- genes.test

  heatmap.test <- pheatmap(zscore.test+1, color = colors, cluster_row = FALSE, cluster_cols = FALSE,
                           gaps_col = s.test, annotation_col = HL.group.test, cutree_col = 2,
                           show_colnames = F, fontsize = 6.5,fontsize_row=6, filename = "zscore.test.pdf")
}
