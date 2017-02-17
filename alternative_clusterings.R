load("genotypesImputed.Rdata")
dim(genotypesImputed)

heatmap(genotypesImputed, scale = "none")
km = kmeans(genotypesImputed,11)

d = dist(genotypesImputed, method="binary")
hc = data.frame(cutree(hclust(d, method="average"), k=11))
