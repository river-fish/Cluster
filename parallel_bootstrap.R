library(parallel)
source("HDP_genomic_fit.R")
load("genotypesImputed.Rdata")
source("Output_to_Posterior.R")

start = Sys.time()

times = 15
N = 8

set.seed(1)

X = mclapply(1:times,function(t){
  
  for(j in 1:times){
    
    if(t==j){
      list_runs <- lapply(1:N, function(i){
        bootstrapped_genotypesImputed <- genotypesImputed[sample(1:nrow(genotypesImputed), nrow(genotypesImputed),  replace = TRUE),]
        a = HDP_genomic_fit(bootstrapped_genotypesImputed,
                            shape=1,invscale=1, #Prior parameters for concentration parameters
                            burnin = 5000, #Burnin for markov chain 
                            postsamples = 10000, #Number of posterior samples
                            spacebw = 20, #space between posterior samples
                            cpsamples = 10,
                            seed=1)
        
        a_merged = hdp_extract_components(a, cos.merge = 0.95)
        
        a_post = posterior_quantities(a_merged)
        
        vec <- apply(a_post$posteriorProbability,2, which.max)
        print(length(rownames(bootstrapped_genotypesImputed)))
        df_post_clust <- data.frame(cluster_id=vec, patient_id = rownames(bootstrapped_genotypesImputed))
        
        save(df_post_clust, file = paste0("bootstrapped_data/",j, "_" , i,"_df_post_clust.RData"))
      }
      )
      return("success")
    }
  }
}

,mc.cores = times)

now = Sys.time()
print(now - start)
