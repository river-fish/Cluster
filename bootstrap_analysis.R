load("genotypesImputed.Rdata")
source("alternative_clusterings.R")
times = 2 #nr of cores
N = 2 #nr of runs in lapply

#load in the cluster assignments from the paper and create dataframe
load("Simply_cluster.Rdata")
df2 = data.frame(cluster_id=as.numeric(as.character(dpClass)), patient_id=rownames(genotypesImputed))

result_concordance_boot <- lapply(1:times, function(j){
  return(sapply(1:N, function(i){
    load(file=paste0(j, "_" , i,"_df_post_clust.RData")) #load the cluster assigments from a bootstrapped dataset
    ConcordanceFunction(df_post_clust, df2, same_patients=FALSE)
  }))
})

# 
# 
# 
# concordance_vect = sapply(1:(ncol(df_post_clust)-1), function(i) {
#   df1 <- df_post_clust[i]
#   colnames(df1) <- c('cluster_id')
#   ConcordanceFunction(df1, df2)
# })
