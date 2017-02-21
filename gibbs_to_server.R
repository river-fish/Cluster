library(MCMCpack)

source('GibbsSampler.R')
load('genotypesImputed.Rdata')
# preparing initial data 
m=84
K = 11
X <- genotypesImputed

beta <- rep(1/K, times=K)
alpha <- rep(1/m, times=m)
w = rep(1/K, times =K)
Z <- sample(1:K, nrow(X), prob=w, replace = TRUE)
theta <- rdirichlet(K, alpha)

N <- 10000
system.time(
result <- FullGibbs(N, X = X, w_start = w, theta_start = theta, Z_start = Z, alpha = alpha, beta=beta)
)
save(result, file=paste0('gibbs_result_',N, '.RData'))
