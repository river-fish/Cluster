library(parallel)
source("HDP_genomic_fit.R")
load("genotypesImputed.Rdata")

start = Sys.time()

times = 4

set.seed(1)

X = mclapply(1:times,function(t){
  if (t==1) {
   
    list_a <- lapply(1:25, function(i){
      bootstrapped_genotypesImputed <- genotypesImputed[sample(1:nrow(genotypesImputed), nrow(genotypesImputed),  replace = TRUE),]
      a = HDP_genomic_fit(bootstrapped_genotypesImputed,
                          shape=1,invscale=1, #Prior parameters for concentration parameters
                          burnin = 5000, #Burnin for markov chain 
                          postsamples = 10000, #Number of posterior samples
                          spacebw = 20, #space between posterior samples
                          cpsamples = 10,
                          seed=1)
      
      save(a, file = paste0(i, "_a.RData"))
      
      return(a)
      
    })
    return(list_a)
    
  }
  
  if (t==2) {
   
    list_b <- lapply(1:25, function(i){
      bootstrapped_genotypesImputed <- genotypesImputed[sample(1:nrow(genotypesImputed), nrow(genotypesImputed),  replace = TRUE),]
      b = HDP_genomic_fit(bootstrapped_genotypesImputed,
                          shape=1,invscale=1, #Prior parameters for concentration parameters
                          burnin = 5000, #Burnin for markov chain 
                          postsamples = 10000, #Number of posterior samples
                          spacebw = 20, #space between posterior samples
                          cpsamples = 10,
                          seed=2)
      
      save(b, file = paste0(i, "_b.RData"))
     
      return(b)
    })
    return(list_b)
  }
  
  if (t==3) {
    list_c <- lapply(1:25, function(i){
      bootstrapped_genotypesImputed <- genotypesImputed[sample(1:nrow(genotypesImputed), nrow(genotypesImputed),  replace = TRUE),]
      c = HDP_genomic_fit(bootstrapped_genotypesImputed,
                          shape=1,invscale=1, #Prior parameters for concentration parameters
                          burnin = 5000, #Burnin for markov chain 
                          postsamples = 10000, #Number of posterior samples
                          spacebw = 20, #space between posterior samples
                          cpsamples = 10,
                          seed=3)
      
      save(c, file = paste0(i, "_c.RData"))
      return(c)
    })
    
    
    
    return(list_c)
  }
  if (t==4) {
    list_d <- lapply(1:25, function(i){
      bootstrapped_genotypesImputed <- genotypesImputed[sample(1:nrow(genotypesImputed), nrow(genotypesImputed),  replace = TRUE),]
      d = HDP_genomic_fit(bootstrapped_genotypesImputed,
                          shape=1,invscale=1, #Prior parameters for concentration parameters
                          burnin = 5000, #Burnin for markov chain 
                          postsamples = 10000, #Number of posterior samples
                          spacebw = 20, #space between posterior samples
                          cpsamples = 10,
                          seed=4)
      
      save(d, file = paste0(i, "_d.RData"))
      
      return(d)
    })
    
    return(list_d)
}}
  
,mc.cores = times)

now = Sys.time()
print(now - start)