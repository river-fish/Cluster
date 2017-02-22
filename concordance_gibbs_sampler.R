<<<<<<< HEAD
n_cores = 6
library(parallel)

getwd()

set.seed(1)

for (nr_clust in c(9)){
  
  load(paste0('gibbs_result_part2_',nr_clust,'_50000.RData'))
  gibbs_output <- result_part1
  
  
  
  start = Sys.time()
  
  X = mclapply(1:n_cores,function(t){
    for(j in 1:n_cores){
      if (t==j) {
        
        concordances_gibbs <- c()
        
        vec <- 1+seq(from=(j-1)*1000, to = (j)*1000, by=100)
        for(nr_gibb in 2:length(vec)){
          
          cluster_id <- gibbs_output[[(vec[nr_gibb-1])]]$Z
          df1 <- data.frame(cluster_id = cluster_id)
          rownames(df1) <- rownames(genotypesImputed)
          
          cluster_id <- gibbs_output[[(vec[nr_gibb])]]$Z
          df2 <- data.frame(cluster_id = cluster_id)
          rownames(df2) <- rownames(genotypesImputed)
          
          concordances_gibbs <- c(concordances_gibbs, ConcordanceFunction(df1, df2))
        }
        
        save(concordances_gibbs, file=paste0('concordances_gibbs_part2',j,'_nr_clust_',nr_clust,'.RData'))
        return(concordances_gibbs)
        
      }
    }
    
  }
  
  ,mc.cores = n_cores)
}

now = Sys.time()
print(now - start)
=======
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
>>>>>>> ca2114da42dbf81cc49ba45c39dc51ec041c8afd
