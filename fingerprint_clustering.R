#clean all objects
rm(list=ls())
setwd("~/bin/fingerprint_clustering")
#global variables
font_size=2
nb_metabolites=9
max_cluster=5
#max_cluster=nb_metabolites/1.5
#choix du niveau de coupure

library(gclus)
#Pseudo-random settings: 
#set.seed(1)
#milisec * PID
set.seed(as.numeric(format(Sys.time(), "%OS2"))*100 * Sys.getpid())

################################################
#     Data test: Random distance matrix   
################################################
rand_distances=ceiling(runif(nb_metabolites*nb_metabolites, 0, 11)) #genration of number between 1 and 12
#conversion into matrix
data = matrix(rand_distances,nb_metabolites, nb_metabolites)
#conversion into triangular matrix
data[upper.tri(data)] = 0
#conversion into distance
distance_matrix=as.dist(data)
#plotcolors(dmat.color(distance_matrix))

#library(vegan)
#distance_matrix=vegdist(data,"jaccard")
################################
#          Clustering
################################
classif=hclust(distance_matrix,method="ward.D2")
#dendro2=hclust(distance_matrix,method="complete")
#dendro3=hclust(distance_matrix,method="median")

#par(mfrow=c(2,2)); plot(dendro);plot(dendro2);plot(dendro3);plot(dendro4);par(mfrow=c(1,1))
#x11();


#automaticly ordering by clusters
classif = reorder.hclust(classif, data)
#plot dendrogram
plot(classif, hang=-1, lwd=font_size,xlab="Metabolites", sub="",ylab="Distance before each fusion",main="Dendrogram of Ward",font.lab=font_size,axes=F)
axis(2, seq(0,ceiling(max(classif$height))),lwd=font_size,font.axis=font_size,cex.axis=0.8)
abline(h=c(classif$height), lty=3, col="grey")
#projection of the clusters
nb_clusters=5
rect.hclust(classif, k=nb_clusters, border=rainbow(nb_clusters))


clusters= function() {
  cutree(classif,nb_clusters)
}
#clusters<-as.factor(cutree(dendro3,nb_clusters))
table(clusters())

"
$merge # negative values: singleton fusion; positive values: cluster fusion
[,1] [,2] 
[1,]   -1   -3 # P1 = {1,3} #fusion of singleton 1 and 3
[2,]   -2    1 # P2 = {1,3,2} 
[3,]   -4    2 # P3 = {1,3,2,4} 
[4,]   -5    3 # P4 = {1,3,2,4,5} 
[5,]   -6   -7 # P5 = {6,7} 
[6,]   -8    5 # P6 = {6,7,8} 
[7,]   -9  -10 # P7 = {9,10} 
[8,]  -11    7 # P8 = {9,10,11} 
[9,]  -12    8 # P9 = {9,10,11,12} 
[10,]  -13    9 # P10 = {9,10,11,12,13} 
[11,]    6   10 # P11 = {6,7,8,9,10,11,12,13} #fusion of cluster at line 6 and the one at line 10
[12,]    4   11 # P12 = {1,3,2,4,5,6,7,8,9,10,11,12,13} "

library(cluster)
si = silhouette(clusters(),distance_matrix)
plot(si)

################################
#          Fusion graph
################################

#Show differences between nodes levels (distance between clusters)
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
d[order(-d), , drop = FALSE]

#(optimal_clusters=which.max(height_diff))
#max(height_diff)
max(classif$height)

label_text=round(height_diff,digits=2)

#Plot fusion graph
#plot_fusion = function() {
  plot(classif$height, nrow(data):2, type="S",main="Distance before each fusion",lwd=font_size,xlim=c(0,max(classif$height)+1),ylim=c(2,nrow(data)),font.lab=2,ylab="Number of clusters", xlab="Node height", col="grey", axes=F)
  axis(2, seq(2,nrow(data)),lwd=2,font.axis=font_size)
  axis(1, seq(0:max(classif$height)),lwd=2,font.axis=font_size)
  text(x=classif$height[-1], y=(nrow(data)-1):2, adj=c(-0.5,-0.5),labels=round(height_diff[-1],digits=2), col="red", cex=0.8)
#}

setGrpContent=function(x,y){
  return (paste("(",x,",",y,")",sep=""))
}
plot_fusion()

group_content=paste("{",paste(seq(asb(i):abs(j)),collapse=","),"}",sep="")

group_correspondance=c(1:(nrow(data)-1))
group_number=0
for (i in 1:(nrow(data)-1)){
  element1=classif$merge[i,1]
  element2=classif$merge[i,2]
  if( element1 < 0 && element2 < 0){
    group_number=group_number+1
    group_name=paste("G",group_number,sep="")
    ordered_singletons=sort(c(abs(element1),abs(element2)))
    group_content=setGrpContent(ordered_singletons[1],ordered_singletons[2])
    assign(group_name,group_content)
    group_correspondance[i]=group_name
  }else if( element1 < 0 && element2 > 0){
    #group_name=paste("G",element2,sep="")
    #group_content=setGrpContent(group_correspondance[i],abs(element1))
    "group_name=group_correspondance[element2]
    group_content=setGrpContent(get(group_name),abs(element1))
    assign(group_name,group_content)
    group_correspondance[i]=group_name"
  }else if( element1 > 0 && element2 > 0){
    "group_number=group_number+1
    ordered_singletons=sort(c(abs(element1),abs(element2)))
    group_correspondance[i]=setGrpContent(group_correspondance[ordered_singletons[1]],group_correspondance[ordered_singletons[2]])
  "}
}

group_correspondance
