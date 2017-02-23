
nr_clust <- 11
comparison_11 <- c()
for(j in 1:6){
  load(paste0('~/Project7/cl/comparison_final_gibbs',j,'_nr_clust_',nr_clust,'.RData'))
  comparison_11 <- c(comparison_11, concordances_gibbs)
}

hist(comparison_11)
mean(comparison_11)


nr_clust <- 10
comparison_10 <- c()
for(j in 1:6){
  load(paste0('~/Project7/cl/comparison_final_gibbs',j,'_nr_clust_',nr_clust,'.RData'))
  comparison_10 <- c(comparison_10, concordances_gibbs)
}
hist(comparison_10)
mean(comparison_10)


nr_clust <- 9
comparison_9 <- c()
for(j in 1:6){
  load(paste0('~/Project7/cl/comparison_final_gibbs',j,'_nr_clust_',nr_clust,'.RData'))
  comparison_9 <- c(comparison_9, concordances_gibbs)
}
hist(comparison_9)
mean(comparison_9)

#-----------------
load("~/Project7/Cluster/genotypesImputed.Rdata")
source('Clust_Image.R')
load("~/Project7/cl/gibbs_result_part2_10_50000.RData")
load("Simply_cluster.Rdata")
library(ggplot2)
library(RColorBrewer)
Cluster_Image(genotypesImputed, result_part2[[100]]$Z)

Cluster_Image(genotypesImputed, result_part2[[1000]]$Z)
Cluster_Image(genotypesImputed, result_part2[[10000]]$Z)
Cluster_Image(genotypesImputed, result_part2[[7]]$Z)
Cluster_Image(genotypesImputed, result_part2[[1300]]$Z)
Cluster_Image(genotypesImputed, result_part2[[666]]$Z)


Cluster_Image(genotypesImputed, as.numeric(as.character(dpClass)))
 
              Impose_NPM1_Clust = T,
              Impose_Fusion_clust=T,
              Impose_CEBPA_bi =T,
              Impose_IDH2_p172 =T,
              Impose_Chromatin_Splice = T,
              Impose_TP53_aneuploidy = TRUE)

load("~/Project7/cl/gibbs_result_part2_9_50000.RData")
clust_data <- data.frame(clust_9 = result_part2[[100]]$Z)

load("~/Project7/cl/gibbs_result_part2_10_50000.RData")
clust_data <- data.frame(clust_10 = result_part2[[100]]$Z)

load("~/Project7/cl/gibbs_result_part2_11_50000.RData")
clust_data <- data.frame(clust_11 = result_part2[[100]]$Z)

clust_data_bayesian <- clust_data
save(clust_data_bayesian, file = 'clust_data_bayesian.RData')
