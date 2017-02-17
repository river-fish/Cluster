load("genotypesImputed.Rdata")
dim(genotypesImputed)

heatmap(genotypesImputed, scale = "none")
km = kmeans(genotypesImputed,11)

