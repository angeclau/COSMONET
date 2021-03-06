
#' @title Split dataset in training and testing set
#'
#' @description This function splits data in training and testing set
#' @param x input matrix [nxp]
#' @param q q% of the sample size
#'
#' @return training and testing set
#'
splitData <- function(x,q){
  
  ## q% of the sample size
  sample_size <- floor(q * nrow(x))
  
  ## set the seed to make your partition reproducible
  set.seed(123)
  ind_training <- sample(seq_len(nrow(x)), size = sample_size)
  
  training <- x[ind_training, ]
  testing <- x[-ind_training, ]
  
  return(list(training=training, testing=testing))
}


