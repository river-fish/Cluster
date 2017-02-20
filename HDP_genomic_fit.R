#Fitting DP mixture model for genomic data as in paper, default parameters are those applied wihin the paper

library(CoxHD) # library(devtools); install_github("mg14/CoxHD/CoxHD")
library(mg14) # library(devtools); install_github("mg14/mg14")
library(hdp) # devtools::install_github("nicolaroberts/hdp", build_vignettes = FALSE)

HDP_genomic_fit <- function(Imputed_genotypes,
                            shape=1,invscale=1, #Prior parameters for concentration parameters
                            burnin = 5000, #Burnin for markov chain 
                            postsamples = 10000, #Number of posterior samples
                            spacebw = 20, #space between posterior samples
                            cpsamples = 10,
                            seed=42){ #Number of concentration mixing step between posterior samples

  n <- ncol(Imputed_genotypes)
  hdp <- hdp_init(ppindex=0, #index of the parent DP for initial DP
                  cpindex=1, #index of alphaa and alphab for initial DP
                  hh=rep(1/n,n), #params for base distn (uniform Dirichlet)
                  alphaa=shape,
                  alphab=invscale)
  
  hdp <- hdp_adddp(hdp,
                   numdp=nrow(Imputed_genotypes), # one DP for every sample in that cancer type
                   pp=1, # parent DP for group i is the i-th+1 overall DP because of the grandparent at position 1
                   cp=1) # index of alphaa and alphab for each DP
  
  # Assign the data from each patient to a child DP
  hdp <- hdp_setdata(hdp = hdp, dpindex=1:nrow(Imputed_genotypes)+1, data=Imputed_genotypes)
  
  
  # Activate the DPs with specified number of classes (signatures)
  hdp <- dp_activate(hdp, 1:(nrow(Imputed_genotypes)+1), 5)
  #' DP parameters
  set.seed(seed)
  output <- hdp_posterior(hdp = hdp,burnin = burnin,n = postsamples,space = spacebw,cpiter = cpsamples)
  return(output)
}
dim(genotypesImputed)
IBP_genomic_fit <- function(Imputed_genotypes,
                            shape=1,invscale=1, #Prior parameters for concentration parameters
                            burnin = 5000, #Burnin for markov chain 
                            postsamples = 10000, #Number of posterior samples
                            spacebw = 20, #space between posterior samples
                            cpsamples = 10,
                            seed=42){ #Number of concentration mixing step between posterior samples
  n <- ncol(Imputed_genotypes)
  hdp <- hdp_init(ppindex=0, #index of the parent DP for initial DP
                  cpindex=1, #index of alphaa and alphab for initial DP
                  hh=rep(1/n,n), #params for base distn (uniform Dirichlet)
                  alphaa=rep(shape,2),
                  alphab=rep(invscale,2)
                  )
  hdp <- hdp_adddp(hdp,
                   numdp=1,
                   pp=1, # parent DP for group i is the i-th+1 overall DP because of the grandparent at position 1
                   cp=1) # index of alphaa and alphab for each DP
  
  hdp <- hdp_adddp(hdp,
                   numdp=nrow(Imputed_genotypes), # one DP for every sample in that cancer type
                   pp=2, # parent DP for group i is the i-th+1 overall DP because of the grandparent at position 1
                   cp=2) # index of alphaa and alphab for each DP
  
  
  
  # Assign the data from each patient to a child DP
  hdp <- hdp_setdata(hdp = hdp, dpindex=1:nrow(Imputed_genotypes)+2, data=Imputed_genotypes)
  
  
  # Activate the DPs with specified number of classes (signatures)
  hdp <- dp_activate(hdp, 1:(nrow(Imputed_genotypes)+2), 5)
  #' DP parameters
  set.seed(seed)
  output <- hdp_posterior(hdp = hdp,burnin = burnin,n = postsamples,space = spacebw,cpiter = cpsamples)
  return(output)
}

