#' @title Co-expression network
#'
#' @description This function calculate the co-expression network
#' @param x input matrix [nxp]
#'
#' @return FL and adjm matrices
#'
#
NetworkMatrix <- function(x){

  # n: number of samples
  # p: number of genes

  n <- dim(x)[1]
  p <- dim(x)[2]

  vec <- matrix(1,n)
  mean <- as.matrix(t(apply(x,2,mean)))
  x <- x - vec%*%mean

  W <- matrix(1,p,p)
  for(i in 1:p){
    for(j in 1:p){
      c <- sum(x[,i]^2)/(sqrt(sum(x[,i]^2))*sqrt(sum(x[,j]^2)))
      W[i,j] <- abs(c)
      W[j,i] <- W[i,j]
    }
  }

  IX <- apply(W, 2, function(x)(order(x,decreasing = TRUE)))
  IXI <- apply(IX, 2, function(x)(order(x,decreasing = FALSE)))

  W <- matrix(0,p,p)
  for(i in 1:p){
    for(j in 1:p){
      if(i==j){W[i,j] <- 0} else {W[i,j] <- 1/(IXI[i,j])/(IXI[j,i])}
    }
  }

  Sum_R <- rowSums(W)
  Sum_C <- colSums(W)
  S <- matrix(0,p,p)

  for(i in 1:p){
    for(j in 1:p){
      S[i,j] <- W[i,j]/sqrt(Sum_R[i])/sqrt(Sum_C[j])
    }
  }

  # Functional linkage network
  FL <- S

  # 0-1 matrix
  adjm <- mat.or.vec(p,p)
  for (i in 1:p){
    for (j in 1:p){
      if(i==j){FL[i,j]=0}
      else if(FL[i,j]>=0.5){adjm[i,j]=1}
      else if(FL[i,j]<0.5){adjm[i,j]=0}
    }
  }

  return(list(FL=FL,adjm=adjm))
}
