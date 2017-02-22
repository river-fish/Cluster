# for (k in c(9,10,11)) {
#   h_df <- hierarchical(genotypesImputed, k=k)
#   path = paste0("clusterings/",k,"_h_df.RData")
#   save(h_df, file=path)
#   
#   km_df <- km(genotypesImputed, k=k)
#   path = paste0("clusterings/",k,"_km_df.RData")
#   save(km_df, file=path)
# }

load("genotypesImputed.Rdata")
load("Simply_cluster.Rdata")
source("alternative_clusterings.R")
paper_df <- data.frame(cluster_id=as.numeric(as.character(dpClass)), patient_id=rownames(genotypesImputed))
rownames(paper_df) <- rownames(genotypesImputed)

l_h <- list()
l_k <- list()

for (k in c(9,10,11)) {
  path <- paste0("clusterings/",k,"_h_df.RData")
  load(path)
  
  path <- paste0("clusterings/",k,"_km_df.RData")
  load(path)
  
  h_conc <- ConcordanceFunction(h_df, paper_df)
  k_conc <- ConcordanceFunction(km_df, paper_df)
  h_conc
  k_conc
  l_k[[length(l_k)+1]] <- k_conc
  print(length(l_k))
  l_h[[length(l_h)+1]] <- h_conc
  print(length(l_h))
}

random_df = paper_df ##with same size clusters

random_concs = lapply(1:10, function(k){
  random_df$cluster_id = sample(1:11,nrow(random_df), replace=TRUE)
  ConcordanceFunction(random_df, paper_df)
})



