hc = hclust(dist(genotypesImputed, method = 'binary'), method = 'average')
plot(hc, 5)
plot(hc)
rect.hclust(hc, k = 3, border = "red")
result <- data.frame(cutree(hc, 11))
hist(result[,1])

set.seed(1)
genotypesImputed_boot <- genotypesImputed[sample(1:nrow(genotypesImputed), size = nrow(genotypesImputed), replace = TRUE),]
hc_boot <- hclust()
head(genotypesImputed_boot)


hierarchical <- function(gene_patient_data){
  d = dist(gene_patient_data, method="binary")
  hc = data.frame(cutree(hclust(d, method="average"), k=11))
  return(hc)
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

