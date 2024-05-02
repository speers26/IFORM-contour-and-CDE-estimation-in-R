################################################################################
##### Script for validating density fits using predictive tail likelihood ######
################################################################################

# functions ---------------------------------------------------------------

get_group = function(data, inds, j){
  data[inds==j]
}

split_training_test = function(data, k, threshold_q){
  #' Function that splits sample into training and test sets, ready for k-fold
  #' cross validation on the extreme portion of the data
  #' @param data matrix of dataset x and y
  #' @param k number of folds
  #' @param threshold_q extreme threshold quantile

  u = quantile(data[,1], threshold_q, names=F)
  x_lower = data[,1][data[,1]<=u]; y_lower = data[,2][data[,1]<=u]
  x_upper = data[,1][data[,1]>u]; y_upper = data[,2][data[,1]>u]

  group_ind = rep(1:k, ceiling(length(x_upper)/k))[1:length(x_upper)]
  x_splits = sapply(1:k, get_group, data=x_upper, inds=group_ind)
  y_splits = sapply(1:k, get_group, data=y_upper, inds=group_ind)

  folds = list()
  for (i in 1:k){
    training_x = c(x_lower, unlist(x_splits[-i]))
    training_y = c(y_lower, unlist(y_splits[-i]))
    testing_x = unlist(x_splits[i])
    testing_y = unlist(y_splits[i])

    folds[[i]] = list(train=matrix(c(training_x, training_y), ncol=2),
                      test=matrix(c(testing_x, testing_y), ncol=2)
                      )

  }
  return(folds)
}


get_permutation = function(data){
  #' returns a d dimensional permuted sample
  #' @param data description
  d = dim(data)
  boot_ind = sample(1:length(data[,1]), length(data[,1]), replace=F)
  boot_data = matrix(data=NA, nrow=d[1], ncol=d[2])
  for (i in 1:d[2]){
    boot_data[, i] = data[,i][boot_ind]
  }
  return(boot_data)
}

get_bands = function(sample, n_bands){
  #' split Y into given number of bands by X values
  #'@keywords internal
  X = sample[,1] ; Y = sample[,2]

  bands_X = eq_split(X, n_bands)$chunks
  bands_Y = eq_split(Y, n_bands)$chunks

  return(list(bandsX = bands_X, bandsY=bands_Y))

}

eq_split = function(X, n){
  #' splits a vector into n equal chunks (with overlap if necessary)
  #'@keywords internal
  k = length(X) %/% n
  X_overflow = X[(k*n+1): length(X)]

  chunks = list()
  for (i in 1:n){
    chunks[[i]] = X[(1+ (i-1)*k): (i*k)]
  }
  return(list(chunks=chunks, overflow=X_overflow))
}

