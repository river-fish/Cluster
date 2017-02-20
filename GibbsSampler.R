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

  # update of weights of clusters 
  w_new <- rdirichlet(1, beta + sapply(1:K, function(i) sum(Z==i)))
  
  # update of Z
  Z_probs <- sapply(1:K, function(k){
    w_new[k]*prod(sapply(1:ncol(theta), function(j) theta[k,j]^colsums[j]))
  })
  Z_new <- sample(1:K, size = length(Z), prob = Z_probs)
  
  # update theta
  #theta_new <- theta
  list_new_thetas <- lapply(1:K, function(k){
    rdirichlet(1, alpha + colSums(X[Z==k,]))
  })
  
  theta_new <- t(matrix(unlist(list_new_thetas), nrow =ncol(theta)))
  
  output <- list(w_new, Z_new, theta_new)
  return(output)
}