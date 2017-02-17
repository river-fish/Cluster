# hierarchical clustering adjusted to our needs
hierarchical <- function(gene_patient_data){
  d = dist(gene_patient_data, method="binary")
  hc = data.frame(cutree(hclust(d, method="average"), k=11))
  colnames(hc) <- 'cluster_id'
  return(hc)
}

# kmeans function adjusted to our needs
km <- function(gene_patient_data, k){
  data.frame(cluster_id = kmeans(gene_patient_data, k)$cluster)
}


df <- hierarchical(genotypesImputed)

N <- 100
for (i in 1:2){
  genotypesImputed_boot <- genotypesImputed[sample(1:nrow(genotypesImputed), size = nrow(genotypesImputed), replace = TRUE),]
  res <- hierarchical(genotypesImputed_boot)
  res[,2] <- rownames(genotypesImputed_boot)
  colnames(res) = c(paste0('clusterization_', i), 'patient_id')
  df <- merge(df, res)
}



