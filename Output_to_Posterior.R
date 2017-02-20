#Takes complex out of Merged clusters and and provides a list of matrices and arrays
posterior_quantities <- function(Merged_clusters){
  
  #Mutation count per cluster
  posteriorSamples <- array(unlist(Merged_clusters@comp_categ_counts), dim=c(dim(Merged_clusters@comp_categ_counts[[1]]),
                                                                             length(Merged_clusters@comp_categ_counts)))
  posteriorSamples <- aperm(posteriorSamples,c(2,3,1))
  
  rownames(posteriorSamples) <- colnames(genotypesImputed)
  colnames(posteriorSamples) <- 1:ncol(posteriorSamples) -1
  
  posteriorMeans <- rowMeans(posteriorSamples, dim=2)
  dim(posteriorMeans)
  posteriorQuantiles <- apply(posteriorSamples, 1:2, quantile, c(0.025,.5,0.975))
  posteriorMode <- apply(posteriorSamples, 1:2, function(x) {t <- table(x); as.numeric(names(t)[which.max(t)])})
  
  
  #Posterior probabilities of patients being members of particlar clusters
  posteriorProbability <- apply(sapply(Merged_clusters@comp_dp_counts[-1], function(y) colMeans(as.matrix(y))),2,function(x) (x+.Machine$double.eps)/sum(x+.Machine$double.eps))
  
  return(list(posteriorSamples = posteriorSamples, #Array, Gene count x Clust_count x Markov_iterations
              posteriorMeans=posteriorMeans, #Matrix, Gene count x Clust_count
              posteriorQuantiles=posteriorQuantiles, #Matrix, Gene count x Clust_count
              posteriorMode=posteriorMode,#Matrix, Gene count x Clust_count
              posteriorProbability=posteriorProbability) #Matrix, Clust_count x obs_count
         ) 
}
posterior_quantities_IBP <- function(Merged_clusters){
  
  #Mutation count per cluster
  posteriorSamples <- array(unlist(Merged_clusters@comp_categ_counts), dim=c(dim(Merged_clusters@comp_categ_counts[[1]]),
                                                                             length(Merged_clusters@comp_categ_counts)))
  posteriorSamples <- aperm(posteriorSamples,c(2,3,1))
  
  rownames(posteriorSamples) <- colnames(genotypesImputed)
  colnames(posteriorSamples) <- 1:ncol(posteriorSamples) -1
  
  posteriorMeans <- rowMeans(posteriorSamples, dim=2)
  dim(posteriorMeans)
  posteriorQuantiles <- apply(posteriorSamples, 1:2, quantile, c(0.025,.5,0.975))
  posteriorMode <- apply(posteriorSamples, 1:2, function(x) {t <- table(x); as.numeric(names(t)[which.max(t)])})
  
  
  #Posterior probabilities of patients being members of particlar clusters
  posteriorProbability <- apply(sapply(Merged_clusters@comp_dp_counts[-(1:2)], function(y) colMeans(as.matrix(y))),2,function(x) (x+.Machine$double.eps)/sum(x+.Machine$double.eps))
  
  return(list(posteriorSamples = posteriorSamples, #Array, Gene count x Clust_count x Markov_iterations
              posteriorMeans=posteriorMeans, #Matrix, Gene count x Clust_count
              posteriorQuantiles=posteriorQuantiles, #Matrix, Gene count x Clust_count
              posteriorMode=posteriorMode,#Matrix, Gene count x Clust_count
              posteriorProbability=posteriorProbability) #Matrix, Clust_count x obs_count
  ) 
}