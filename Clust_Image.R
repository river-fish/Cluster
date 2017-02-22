Cluster_Image <- function(genotypesImputed,class_1,
                          Impose_NPM1_Clust = FALSE,
                          Impose_Fusion_clust=FALSE,
                          Impose_CEBPA_bi =FALSE,
                          Impose_IDH2_p172 =FALSE,
                          Impose_Chromatin_Splice = FALSE,
                          Impose_TP53_aneuploidy = FALSE){
  
  Imposed_included <- Impose_NPM1_Clust | Impose_Fusion_clust | Impose_CEBPA_bi |Impose_IDH2_p172 | 
    Impose_Chromatin_Splice | Impose_TP53_aneuploidy
  
  par(bty="n", mgp = c(2,.5,0), mar=c(4,4,4,4)+.1, las=2, tcl=-.33)
  #current classification
  #class_1 <- factor(k_means2$cluster -1)
  #class_1 <- dpClass
  
  custom_class <- as.numeric(as.character(class_1))
  total_clust_count <- max(custom_class)
  gene_names <- colnames(genotypesImputed)
  
  
  ##Order known clusters by quantity
  #Clusters found within study
  Fusion_clusts <- c("inv16_t16_16","t_15_17","t_8_21","inv3_t3_3","t_6_9","t_v_11")
  Fusion_clusts <- 
    Fusion_clusts[order(
      as.numeric(
        colSums(genotypesImputed[,Fusion_clusts])
      )
      ,decreasing = FALSE)]
  
  #All of these to be in a path way, with NPM1 being final stage
  NPM1_clust <- "NPM1"
  NPM1_path <- c("DNMT3A","TET2","ASXL1","IDH1","IDH2_p140")
  #Stand alone as mutally exclusive 
  IDH2_p172_clust <- "IDH2_p172"
  CEBPA_biallele <- "CEBPA_bi"
  
  ##Groups to find within our clusters !
  Chromatin_Splice_clust <- c("SFRS2","SF3B1","U2AF1","ZRSR2","STAG2","BCOR",
                              "MLL", "EZH2", "PHF6","RUNX1")
  #CNA
  TP53_aneuploidy_clust <- c("TP53","minus5_5q","minus7","minus7q","abn7other","plus8_8q","minus9q",
                             "mono12_12p_abn12p","plus13","mono17_17p_abn17p","minus18_18q","minus20_20q","plus21",
                             "plus22","minusY","abn3q_other","plus11_11q","mono4_4q_abn4q","complex"
  )
  
  Chromatin_Splice_clust %in% gene_names
  TP53_aneuploidy_clust %in% gene_names
  
  if(Impose_Chromatin_Splice){
    custom_class[rowSums(genotypesImputed[,colnames(genotypesImputed) %in% Chromatin_Splice_clust]) >= 1] <- total_clust_count + 1
    total_clust_count <- Chroma_clust_lab <- total_clust_count + 1
    
  }
  if(Impose_TP53_aneuploidy){
    custom_class[rowSums(genotypesImputed[,colnames(genotypesImputed) %in% TP53_aneuploidy_clust]) >= 1] <- total_clust_count + 1
    total_clust_count <- TP53_cust_labs <-  total_clust_count + 1
  }
  
  if(Impose_IDH2_p172){
    custom_class[genotypesImputed[,colnames(genotypesImputed) == IDH2_p172_clust] == 1] <- total_clust_count + 1
    total_clust_count <- total_clust_count + 1  
  }
  if(Impose_CEBPA_bi){
    custom_class[genotypesImputed[,colnames(genotypesImputed) == CEBPA_biallele] == 1] <- total_clust_count + 1
    total_clust_count <- total_clust_count + 1
  }
  if(Impose_NPM1_Clust){
    custom_class[genotypesImputed[,colnames(genotypesImputed) == NPM1_clust] == 1] <- total_clust_count + 1
    total_clust_count <- total_clust_count + 1
  }
  if(Impose_Fusion_clust){
    for(k in 1:length(Fusion_clusts)){
      custom_class[genotypesImputed[,colnames(genotypesImputed) == Fusion_clusts[k]] == 1] <- total_clust_count + 1
      total_clust_count <- total_clust_count + 1
    }
  }
  
  
  
  
  
  print(paste("Total Cluster Count: ",length(unique(custom_class))))
  
  #Update cluster count incase some become empty
  o <- order(custom_class,decreasing = TRUE)
  
  custom_class <- factor(custom_class)
  total_clust_count <- length(unique(custom_class))
  levels(custom_class) <- 0:(total_clust_count -1)
  
  
  # 
  
  
  
  ##Order genes by their cluster assigments within paper
  
  gene_o <- c(which( (!colnames(genotypesImputed) %in% c(Fusion_clusts,NPM1_clust,NPM1_path,IDH2_p172_clust,CEBPA_biallele,
                                                         Chromatin_Splice_clust,TP53_aneuploidy_clust) )  &
                       as.numeric(colSums(genotypesImputed) > 20 )
  ),
  which(colnames(genotypesImputed) %in% Chromatin_Splice_clust)[
    order(colSums(genotypesImputed[,which(colnames(genotypesImputed) %in% Chromatin_Splice_clust)]))
    ],
  which(colnames(genotypesImputed) %in% TP53_aneuploidy_clust)[
    order(colSums(genotypesImputed[,which(colnames(genotypesImputed) %in% TP53_aneuploidy_clust)]))
    ],
  which(colnames(genotypesImputed) %in% CEBPA_biallele),
  which(colnames(genotypesImputed) %in% IDH2_p172_clust),
  which(colnames(genotypesImputed) %in% NPM1_clust),
  which(colnames(genotypesImputed) %in% NPM1_path[-1]),
  which(colnames(genotypesImputed) == "DNMT3A"), 
  which(colnames(genotypesImputed) %in% Fusion_clusts)[
    order(colSums(genotypesImputed[,which(colnames(genotypesImputed) %in% Fusion_clusts)]))
    ]
  )
  
  xlength <- dim(genotypesImputed)[1]
  ylength <- length(gene_o)
  
  print(length(o))
  
  brewer_cols <- c(brewer.pal(3,"Blues"))
  image(genotypesImputed[o,gene_o],
        col=brewer_cols[c(1,length(brewer_cols))],
        xaxt="n", yaxt="n", xlab="",ylab="",
        x=1:xlength, y=1:ylength
  )
  ylabs <- colnames(genotypesImputed)[gene_o]
  ylabs_cols <- rep(1,length(colnames(genotypesImputed)))
  names(ylabs_cols) <- colnames(genotypesImputed)
  ##Colour genes by their known properties i.e Fusion or not
  ylabs_cols[Fusion_clusts] <- "Brown"
  ylabs_cols[NPM1_clust] <- "Red"
  ylabs_cols[IDH2_p172_clust] <- "Blue"
  ylabs_cols[CEBPA_biallele] <- "Green"
  ylabs_cols[Chromatin_Splice_clust] <- "Purple"
  ylabs_cols[TP53_aneuploidy_clust] <- "Orange"
  ylabs_cols[NPM1_path] <- "DarkGreen"
  
  #Verticle labels
  mtext(side=2, at=1:ylength, ylabs, cex=.66, col=ylabs_cols[gene_o]) 
  #Vertical lines
  hori_change_points <- which(diff(as.numeric(as.character(custom_class[o]))) != 0)
  abline(v=hori_change_points, col="black", lwd=.5)
  
  #Add line to show cutoff from where imposed clusters stop
  if(Imposed_included){
    number_of_changepoints <- length(Fusion_clusts)*Impose_Fusion_clust + Impose_NPM1_Clust + Impose_CEBPA_bi +
      Impose_IDH2_p172 + Impose_TP53_aneuploidy + Impose_Chromatin_Splice
    abline(v=hori_change_points[number_of_changepoints], col="purple", lwd=2)
  }
  
  #Hori Lines - marking from when clusters must be inferred
  abline(h=ylength - length(Fusion_clusts) + 0.5)
  abline(h=ylength - length(Fusion_clusts) - length(NPM1_path) - length(NPM1_clust) + 0.5)
  
  #Axis labels
  labs <- c(rev(Fusion_clusts),NPM1_clust,CEBPA_biallele,IDH2_p172_clust,"Chromatin_Splice","TP53_aneuploidy")
  labs <- labs[c(rep(Impose_Fusion_clust,length(Fusion_clusts)),Impose_NPM1_Clust,Impose_CEBPA_bi,Impose_IDH2_p172,
                 Impose_Chromatin_Splice,Impose_TP53_aneuploidy)]
  
  labs <- c(labs,paste("Cluster",(1:(total_clust_count  - length(labs))  + length(labs))  ))
  mtext(side=1, at= c(0,hori_change_points),
        labs, cex=.66)
}
