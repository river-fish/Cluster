SingleGibbs <- function(beta, alpha, w, theta, Z, X, K=10){
  require(MCMCpack)
  # beta - number
  # alpha - number
  # w - vector of length K (nr of clusters)
  # theta - matrix with K rows and 84 columns
  # Z - a vector of length 1540 (nr of patients); contains integers between 1 and 10
  # X - data
  
  colsums <- colSums(X)
  K = length(w)
   print (K)
  # update of weights of clusters 
  w_new <- rdirichlet(1, beta + sapply(1:K, function(i) sum(Z==i)))
  
  print(w_new)
  # update of Z
  Z_log_probs <- sapply(1:K, function(k){
    log(w_new[k]) + sum(sapply(1:ncol(theta), function(j) log(theta[k,j])*colsums[j]))
  })
  print(Z_log_probs)
  Z_new <- sample(1:K, size = length(Z), prob = exp(Z_log_probs - max(Z_log_probs)+1), replace = TRUE)
  
  # update theta
  #theta_new <- theta
  list_new_thetas <- lapply(1:K, function(k){
    rdirichlet(1, alpha + colSums(X[Z_new==k,]))
  })
  
  theta_new <- t(matrix(unlist(list_new_thetas), nrow =ncol(theta)))
  
  output <- list(w_new, Z_new, theta_new)
  return(output)
}

library(MCMCpack)
m=84
K = 10
X <- genotypesImputed

beta <- rep(1/K, times=K)
alpha <- rep(1/m, times=m)
w = rep(1/10, times =10)
Z <- sample(1:K, nrow(X), prob=w, replace = TRUE)
theta <- rdirichlet(K, alpha)

out <- SingleGibbs(beta, alpha, w, theta, Z, X, K=10)






