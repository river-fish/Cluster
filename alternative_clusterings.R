binary_distance <- function(gene_patient_data){
  return(dist(gene_patient_data, method="binary"))
}

# hierachical clustering on binary data
hierarchical <- function(gene_patient_data, k=11){
  hc = data.frame(cutree(hclust(binary_distance(gene_patient_data), method="average"), k=k))
  colnames(hc) <- 'cluster_id'
  hc$patient_id <- rownames(gene_patient_data)
  return(hc)
}

# kmeans function adjusted to our needs
km <- function(gene_patient_data, k=11){
  data.frame(cluster_id = kmeans(gene_patient_data, k)$cluster, patient_id = rownames(gene_patient_data))
}


df1 <- hierarchical(genotypesImputed)
head(df1)
df2 <- km(genotypesImputed)

table(df1[,1], df2[,1])

k <- length(table(summary(df1$cluster_id)))
print(k)


ConcordanceFunction<- function(df1, df2){
 
  k=11
  concordances_in_clusters <- sapply(1:k, function(i){
    cluster_patient_ids <- df1[df1$cluster_id==i,]$patient_id
    df_temp <- df2[df2$patient_id %in% cluster_patient_ids,]
    people_together <- unname(summary(as.factor(df_temp$cluster_id)))
    sum(sapply(people_together, function(n) choose(n,2)))
  })
  return(sum(concordances_in_clusters))
}

ConcordanceFunction(df2,df1)

summary(as.factor(df1$cluster_id))

N <- 100
for (i in 1:2){
  genotypesImputed_boot <- genotypesImputed[sample(1:nrow(genotypesImputed), size = nrow(genotypesImputed), replace = TRUE),]
  res <- hierarchical(genotypesImputed_boot)
  res[,2] <- rownames(genotypesImputed_boot)
  colnames(res) = c(paste0('clusterization_', i), 'patient_id')
  df <- merge(df, res)
}



