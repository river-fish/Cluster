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


# ConcordanceFunction_ver2<- function(df1, df2){
#  
#   k=11
#   concordances_in_clusters <- sapply(1:k, function(i){
#     cluster_patient_ids <- df1[df1$cluster_id==i,]$patient_id
#     df_temp <- df2[df2$patient_id %in% cluster_patient_ids,]
#     people_together <- unname(summary(as.factor(df_temp$cluster_id)))
#     sum(sapply(people_together, function(n) choose(n,2)))
#   })
#   return(sum(concordances_in_clusters))
# }
# 
# ConcordanceFunction_ver2(df2,df1)
# 
# summary(as.factor(df1$cluster_id))

# N <- 100
# for (i in 1:2){
#   genotypesImputed_boot <- genotypesImputed[sample(1:nrow(genotypesImputed), size = nrow(genotypesImputed), replace = TRUE),]
#   res <- hierarchical(genotypesImputed_boot)
#   res[,2] <- rownames(genotypesImputed_boot)
#   colnames(res) = c(paste0('clusterization_', i), 'patient_id')
#   df <- merge(df, res)
# }


ConcordanceFunction <- function(df1, df2){
  #df1, df2 - data frames containing a column cluster_id
  if(nrow(df1)!=nrow(df2)){
    stop('The number of rows in the data frames are not equal')
  }
  df2_ordered <- df2[rownames(df1),]
  counts <- sapply(seq_len(nrow(df1)), function(k){
    count = 0
    for (i in seq_len(k-1)){
      if((df1$cluster_id[k] == df1$cluster[i]) & (df2_ordered$cluster_id[k] == df2_ordered$cluster[i])) count = count +1
      else if ((df1$cluster_id[k] != df1$cluster[i]) & (df2_ordered$cluster_id[k] != df2_ordered$cluster[i])) count = count +1
    }
    return(count)
  })
  return(sum(counts)/choose(nrow(df1),2))
}

ConcordanceFunction(df1[1:1000,], df1[1:1000,])



# ConcordanceFunctionF <- function(df1, df2){
#   if(nrow(df1)!=nrow(df2)){
#     stop('The number of rows in the data frames are not equal')
#   }
#   df2_ordered <- df2[rownames(df1),]
#   counts <- sapply(2:nrow(df1), function(k){
#     
#     single_count <- sapply(seq_len(k-1), function(i) {
#       if((df1$cluster_id[k] == df1$cluster[i]) & (df2_ordered$cluster_id[k] == df2_ordered$cluster[i])) return(1)
#       else if ((df1$cluster_id[k] != df1$cluster[i]) & (df2_ordered$cluster_id[k] != df2_ordered$cluster[i])) return(1)
#       else return(0)
#     })
#     return(sum(single_count))
#   })
#   return(sum(counts)/choose(nrow(df1),2))
# }
# 
# 
# system.time(ConcordanceFunctionF(df1[1:1000,], df1[1:1000,]))
system.time(ConcordanceFunction(df1[1:1000,], df1[1:1000,]))
