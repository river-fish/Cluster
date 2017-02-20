# single Gibb's sampler step
SingleGibbs <- function(beta, alpha, w, theta, Z, X){
  require(MCMCpack)
  # beta - number
  # alpha - number
  # w - vector of length K (nr of clusters)
  # theta - matrix with K rows and 84 columns
  # Z - a vector of length 1540 (nr of patients); contains integers between 1 and 10
  # X - data
  
  colsums <- colSums(X)
  K = length(w) # nr of clusters
  n = nrow(X) # nr of patients

  # update of weights of clusters 
  w_new <- rdirichlet(1, beta + sapply(1:K, function(i) sum(Z==i)))
  
 
  # update of Z
  
  # Z_log_probs <- sapply(1:K, function(k){
  #   log(w_new[k]) + sum(sapply(1:ncol(theta), function(j) log(theta[k,j])*colsums[j]))
  # })
  # print(Z_log_probs)
  # Z_new <- sample(1:K, size = length(Z), prob = exp(Z_log_probs - max(Z_log_probs)+1), replace = TRUE)
  # 
  
  list_probs <- lapply(1:n, function(i){
    # want to create a vector of length K
    sapply(1:K, function(k) w_new[k]*prod(theta[k,]^X[i,]))
  })
    
  Z_new <- unlist(lapply(1:n, function(i) sample(1:K, size = 1, prob= list_probs[[i]])))
  # update theta
  #theta_new <- theta
  list_new_thetas <- lapply(1:K, function(k){
    rdirichlet(1, alpha + colSums(X[Z_new==k,]))
  })
  
  theta_new <- t(matrix(unlist(list_new_thetas), nrow =ncol(theta)))
  
  output <- list(w=w_new, Z=Z_new, theta=theta_new)
  return(output)
}

# ----------------------------------------------------------------
# the full function (n_iter iterations of the sampler)
FullGibbs <- function(n_iter, X, w_start, theta_start, Z_start, alpha = rep(1/84, 84), beta = rep(1/10, 10)){
  output_list <- list()
  for (i in 1:n_iter){
    output <- SingleGibbs(beta, alpha, w, theta, Z, X)
    w <- out$w
    theta <- out$theta
    Z <- out$Z
    output_list[[i]] <- output
  }
  return(output_list)
}

# ----------------------------------------------------------------
# preparing initial data 
m=84
K = 10
X <- genotypesImputed

beta <- rep(1/K, times=K)
alpha <- rep(1/m, times=m)
w = rep(1/10, times =10)
Z <- sample(1:K, nrow(X), prob=w, replace = TRUE)
theta <- rdirichlet(K, alpha)


result <- FullGibbs(100, X = X, w_start = w, theta_start = theta, Z_start = Z)
length(result)
summary(as.factor(result[[99]]$Z))
