setwd("~/bin/fingerprint_clustering")
#global variables
nb_metabolites=9


#Pseudo-random settings: 
#set.seed(1)
#milisec * PID
set.seed(as.numeric(format(Sys.time(), "%OS2"))*100 * Sys.getpid())

################################################
#     Data test: Random distance matrix   
################################################
rand_distances=floor(runif(nb_metabolites*nb_metabolites, 1, 12)) #genration of number between 1 and 12
#conversion into matrix
distance_matrix = matrix(rand_distances,nb_metabolites, nb_metabolites)
#conversion into triangular matrix
distance_matrix[upper.tri(distance_matrix)] = 0
#conversion into distance
distance_matrix=as.dist(distance_matrix)

#library(vegan)
#distance_matrix=vegdist(data,"jaccard")
################################
#          Clustering
################################
dendro1=hclust(distance_matrix,method="ward.D2")
dendro2=hclust(distance_matrix,method="complete")
dendro3=hclust(distance_matrix,method="average")
dendro4=hclust(distance_matrix,method="single")

par(mfrow=c(2,2)); plot(dendro1);plot(dendro2);plot(dendro3);plot(dendro4);par(mfrow=c(1,1))
x11();plot(dendro1)
rect.hclust(spe.chwo, k=k, border=rainbow(k))

clusters<-as.factor(cutree(dendro3,3))
table(clusters)

