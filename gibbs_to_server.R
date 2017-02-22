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

N <- 50000
system.time(
result_part1 <- FullGibbs(N, X = X, w_start = w, theta_start = theta, Z_start = Z, alpha = alpha, beta=beta)
)
save(result_part1, file=paste0('gibbs_result_part1_',K,'_',N, '.RData'))
print('hello')

system.time(
  result_part2 <- FullGibbs(N, X = X, w_start = result_part1[[N]]$w, theta_start = result_part1[[N]]$theta, Z_start = result_part1[[N]]$Z, alpha = alpha, beta=beta)
)
save(result_part2, file=paste0('gibbs_result_part2_',K,'_',N, '.RData'))
print('hello2')

##################### 10 clusters
m=84
K = 10
X <- genotypesImputed

beta <- rep(1/K, times=K)
alpha <- rep(1/m, times=m)
w = rep(1/K, times =K)
Z <- sample(1:K, nrow(X), prob=w, replace = TRUE)
theta <- rdirichlet(K, alpha)

N <- 50000
system.time(
  result_part1 <- FullGibbs(N, X = X, w_start = w, theta_start = theta, Z_start = Z, alpha = alpha, beta=beta)
)
save(result_part1, file=paste0('gibbs_result_part1_',K,'_',N, '.RData'))
print('hello')

system.time(
  result_part2 <- FullGibbs(N, X = X, w_start = result_part1[[N]]$w, theta_start = result_part1[[N]]$theta, Z_start = result_part1[[N]]$Z, alpha = alpha, beta=beta)
)
save(result_part2, file=paste0('gibbs_result_part2_',K,'_',N, '.RData'))
print('hello2')


################## 9 clusters

m=84
K = 9
X <- genotypesImputed

beta <- rep(1/K, times=K)
alpha <- rep(1/m, times=m)
w = rep(1/K, times =K)
Z <- sample(1:K, nrow(X), prob=w, replace = TRUE)
theta <- rdirichlet(K, alpha)

N <- 50000
system.time(
  result_part1 <- FullGibbs(N, X = X, w_start = w, theta_start = theta, Z_start = Z, alpha = alpha, beta=beta)
)
save(result_part1, file=paste0('gibbs_result_part1_',K,'_',N, '.RData'))
print('hello')

system.time(
  result_part2 <- FullGibbs(N, X = X, w_start = result_part1[[N]]$w, theta_start = result_part1[[N]]$theta, Z_start = result_part1[[N]]$Z, alpha = alpha, beta=beta)
)
save(result_part2, file=paste0('gibbs_result_part2_',K,'_',N, '.RData'))
print('hello2')