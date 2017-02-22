n_cores = 2
gibbs_output <- result

set.seed(1)

X = mclapply(1:n_cores,function(t){
  
  if (t==1) {
    j = t
    concordances_gibbs <- c()
    
    vec <- seq(from=j*1000, to = (j+1)*1000, by=500)
    for(nr_gibb in vec[-1]){
      
      cluster_id <- gibbs_output[[nr_gibb-1]]$Z
      df1 <- data.frame(cluster_id = cluster_id)
      rownames(df1) <- rownames(genotypesImputed)
     
      cluster_id <- gibbs_output[[nr_gibb]]$Z
      df2 <- data.frame(cluster_id = cluster_id)
      rownames(df2) <- rownames(genotypesImputed)
      
      concordances_gibbs <- c(concordances_gibbs, ConcordanceFunction(df1, df2))
    }
    
    save(concordances_gibbs, file=paste0('concordances_gibbs_',j,'.RData'))
    return(concordances_gibbs)
    
  }
  
  if (t==2) {
    j = t
    concordances_gibbs <- c()
    
    vec <- seq(from=j*1000, to = (j+1)*1000, by=500)
    for(nr_gibb in vec[-1]){
      
      cluster_id <- gibbs_output[[nr_gibb-1]]$Z
      df1 <- data.frame(cluster_id = cluster_id)
      rownames(df1) <- rownames(genotypesImputed)
      
      cluster_id <- gibbs_output[[nr_gibb]]$Z
      df2 <- data.frame(cluster_id = cluster_id)
      rownames(df2) <- rownames(genotypesImputed)
      
      concordances_gibbs <- c(concordances_gibbs, ConcordanceFunction(df1, df2))
    }
    
    save(concordances_gibbs, file=paste0('concordances_gibbs_',j,'.RData'))
    return(concordances_gibbs)

  }
}
  
  ,mc.cores = n_cores)