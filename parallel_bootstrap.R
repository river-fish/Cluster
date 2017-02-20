library(parallel)
source("HDP_genomic_fit.R")
load("genotypesImputed.Rdata")

start = Sys.time()

times = 4

X = mclapply(1:times,function(t){
  if (t==1) {
    set.seed(1)
    bootstrapped_genotypesImputed <- genotypesImputed[sample(1:nrow(genotypesImputed), nrow(genotypesImputed),  replace = TRUE),]
    a = HDP_genomic_fit(bootstrapped_genotypesImputed,
                                shape=1,invscale=1, #Prior parameters for concentration parameters
                                burnin = 5000, #Burnin for markov chain 
                                postsamples = 10000, #Number of posterior samples
                                spacebw = 20, #space between posterior samples
                                cpsamples = 10,
                                seed=1)
    save(a, "a.Rdata")
    return(a)
  }
  
  if (t==2) {
    set.seed(2)
    bootstrapped_genotypesImputed <- genotypesImputed[sample(1:nrow(genotypesImputed), nrow(genotypesImputed),  replace = TRUE),]
    b = HDP_genomic_fit(bootstrapped_genotypesImputed,
                        shape=1,invscale=1, #Prior parameters for concentration parameters
                        burnin = 5000, #Burnin for markov chain 
                        postsamples = 10000, #Number of posterior samples
                        spacebw = 20, #space between posterior samples
                        cpsamples = 10,
                        seed=2)
    save(b, "b.Rdata")
    return(b)
  }
  
  if (t==3) {
    set.seed(3)
    bootstrapped_genotypesImputed <- genotypesImputed[sample(1:nrow(genotypesImputed), nrow(genotypesImputed),  replace = TRUE),]
    c = HDP_genomic_fit(bootstrapped_genotypesImputed,
                        shape=1,invscale=1, #Prior parameters for concentration parameters
                        burnin = 5000, #Burnin for markov chain 
                        postsamples = 10000, #Number of posterior samples
                        spacebw = 20, #space between posterior samples
                        cpsamples = 10,
                        seed=3)
    save(c, "c.Rdata")
    return(c)
  }
  
  if (t==4) {
    set.seed(4)
    bootstrapped_genotypesImputed <- genotypesImputed[sample(1:nrow(genotypesImputed), nrow(genotypesImputed),  replace = TRUE),]
    d = HDP_genomic_fit(bootstrapped_genotypesImputed,
                        shape=1,invscale=1, #Prior parameters for concentration parameters
                        burnin = 5000, #Burnin for markov chain 
                        postsamples = 10000, #Number of posterior samples
                        spacebw = 20, #space between posterior samples
                        cpsamples = 10,
                        seed=4)
    save(d, "d.Rdata")
    return(d)
  }
}
  
,mc.cores = times)

now = Sys.time()
print(now - start)