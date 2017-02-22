setwd("/homes/loach/Cluster")
wd = getwd()
print(wd)
library(parallel)
source("HDP_genomic_fit.R")
print("The file is running")
load("genotypesImputed.Rdata")
print("genotypes loaded")
source("Output_to_Posterior.R")

start = Sys.time()

times = 2
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

load("genotypesImputed.Rdata")
source("alternative_clusterings.R")

#load in the cluster assignments from the paper and create dataframe
load("Simply_cluster.Rdata")
df2 = data.frame(cluster_id=as.numeric(as.character(dpClass)), patient_id=rownames(genotypesImputed))

result_concordance_boot <- lapply(1:times, function(j){
  return(sapply(1:N, function(i){
    load(file=paste0("bootstrapped_data/",j, "_" , i,"_df_post_clust.RData")) #load the cluster assigments from a bootstrapped dataset
    ConcordanceFunction(df_post_clust, df2, same_patients=FALSE)
  }))
})

save(unlist(result_concordance_boot), "bootstrap_concordances.RData")

hist(unlist(result_concordance_boot))
