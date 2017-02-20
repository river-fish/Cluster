#Clustering functions

#Creating an R file to focus on the clustering steps taken in paper name et al. (lol)
#Simply taken their code (DirichletClasses.R and data_clean.R) and tried to put into functions so not as verbose
#

#Packages
library(CoxHD) # library(devtools); install_github("mg14/CoxHD/CoxHD")
library(mg14) # library(devtools); install_github("mg14/mg14")
library(hdp) # devtools::install_github("nicolaroberts/hdp", build_vignettes = FALSE)

#Data - cleaned as in
#source("Our Analysis/data_clean.R")
load("AML-2017-02-16.RData") #load ready cleaned

#see  - "Our Analysis/Data_Description.R" for description of groups of variables
dataFrame #Data set 
groups #variables assigned to groups 
levels(groups)


#Impute genes missing genes - is this robust?
d <- as.matrix(dataFrame[groups %in% c("Genetics","Fusions","CNA")])
d[!d %in% c(0,1)] <- NA
#Number of genes don't take values within 0 or 1
sum(is.na(d))
genotypesImputed <- as.matrix(round(ImputeMissing(d)))


####################################################################
#Clustering methods  - 
####################################################################

##replicating what was done in paper name et al. (lol)
source("Our Analysis/HDP_genomic_fit.R")

##Fits Dirchlet process mixture model  as in paper, default options are their settings
#output <- HDP_genomic_fit(genotypesImputed) ---- Runs markov chain (approx 25 mins),

##Merge clusters that are close togther
#posteriorMerged <- hdp_extract_components(output,cos.merge=0.95) ---- (approx 30 mins)

##load data from previous code  
load("HDP_MCMC_1.Rdata")

###############################################################
#Investigating Clustering Methods------------------------
###############################################################

#First:
#Cleaning of HDP data

#output  & posteriorMerged are S3 objects 
#Two ways of viewing output:
#Mutation count per cluster: matrix(gene_count x clust_count )
#Observations distribution of being in each cluster: matrix(clust_count x patient_count)


#Function taking S3 object and making a list of nice things!
source("Our Analysis/Output_to_Posterior.R")
#Dimension of output (genecount x clust_count)
posteriorList <- posterior_quantities(posteriorMerged)




###############################################################
#Visualisation of output--------------------
#will likely need to differing comparision plots 
###############################################################


#Assign each observation to the cluster which gives it highest posterior probability 
dpClass <- factor(apply(posteriorList$posteriorProbability, 2, which.max)-1)
table(dpClass) #Number of observations in each Cluster

#Save these clusters 

save(posteriorList,dpClass,dataFrame,genotypesImputed,file="Our Analysis/Outputted_clusters/DP_clust_1.Rdata")
save(dpClass,file="Our Analysis/Outputted_clusters/Simply_cluster.Rdata")




#Load plot functions
source("Our Analysis/Descriptive_functions.R")
#Return those genes with highest values in each cluster
Genes_in_clusts(posteriorList$posteriorMeans)

Boxplot_proportions(posteriorList$posteriorProbability)


Prevelance_plot(genotypesImputed,dpClass)

Barplot_gene(genotypesImputed,posteriorList$posteriorQuantiles,dpClass)



#' ### Clinical associations
#+ DPclinical
boxplot(clinicalData$wbc ~ factor(dpClass), log='y', xlab="Class",ylab="wbc", col=col)
boxplot(clinicalData$BM_Blasts ~ factor(dpClass), xlab="Class",ylab="Blast %", col=col)
boxplot(clinicalData$AOD ~ factor(dpClass), xlab="Class",ylab="Age", col=col)
boxplot(rowSums(genotypesImputed) ~ factor(dpClass), xlab="Class",ylab="# Mutations", col=col)









































