#clean all objects
rm(list=ls())
setwd("~/bin/fingerprint_clustering")
#global variables
nb_metabolites=9
max_cluster=5
#max_cluster=nb_metabolites/1.5
  
library(gclus)
dist
#Pseudo-random settings: 
#set.seed(1)
#milisec * PID
set.seed(as.numeric(format(Sys.time(), "%OS2"))*100 * Sys.getpid())

################################################
#     Data test: Random distance matrix   
################################################
rand_distances=floor(runif(nb_metabolites*nb_metabolites, 1, 12)) #genration of number between 1 and 12
#conversion into matrix
data = matrix(rand_distances,nb_metabolites, nb_metabolites)
#conversion into triangular matrix
data[upper.tri(data)] = 0
#conversion into distance
distance_matrix=as.dist(data)

#library(vegan)
#distance_matrix=vegdist(data,"jaccard")
################################
#          Clustering
################################
classif=hclust(distance_matrix,method="ward.D2")
#dendro2=hclust(distance_matrix,method="complete")
#dendro3=hclust(distance_matrix,method="average")
#dendro4=hclust(distance_matrix,method="single")

#par(mfrow=c(2,2)); plot(dendro);plot(dendro2);plot(dendro3);plot(dendro4);par(mfrow=c(1,1))
#x11();

nb_clusters=5
#automaticly ordering by clusters
classif = reorder.hclust(classif, data)
#plot dendrogram
plot(classif, hang=-1, xlab=paste(nb_clusters," groups"), sub="",ylab="Height",main="Dendrogram of Ward",labels=cutree(classif, k=nb_clusters))
#projection of the clusters
rect.hclust(classif, k=nb_clusters, border=rainbow(nb_clusters))

clusters= function() {
  cutree(classif,nb_clusters)
}
#clusters<-as.factor(cutree(dendro3,nb_clusters))
table(clusters())

library(cluster)
si = silhouette(clusters(),distance_matrix)
plot(si)

################################
#          Fusion graph
################################

#Show differences between nodes levels
height_classif=as.matrix(classif$height)
height_diff=matrix(0, length(height_classif), 1)
for (i in 2:(length(height_classif))){
  height_diff[i,]=height_classif[i,]-height_classif[i-1,]
}
#matrix_height=cbind(height_classif,height_diff)
#colnames(matrix_height)=c("node height","difference")
#rownames(matrix_height)=c((length(height_classif)+1):2)
#matrix_height

rownames(height_diff)=c((length(height_classif)+1):2)
d=data.frame(height_diff)
height_diff[order(-d$data), , drop = FALSE]
d[order(-d), , drop = FALSE]

#(optimal_clusters=which.max(height_diff))
#max(height_diff)

#Plot fusion graph
plot(classif$height, nrow(data):2, type="S",main="Fusion levels - Ward",ylab="Number of clusters", xlab="Node height", col="grey")
text(classif$height, nrow(data):2, nrow(data):2, col="red", cex=0.8)
