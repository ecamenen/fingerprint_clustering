#clean all objects
rm(list=ls())
setwd("~/bin/fingerprint_clustering")

#global variables
font_size=2
nb_metabolites=9
max_cluster=8
margin=par(mar=c(5, 4, 4, 2) + 1.1)
interval=1
#max_cluster=nb_metabolites/1.5
#choix du niveau de coupure
def=par()

library(gclus)
library(cluster)

#Pseudo-random settings: 
#set.seed(1)
#milisec * PID
set.seed(as.numeric(format(Sys.time(), "%OS2"))*100 * Sys.getpid())

################################################
#     Data test: random distance matrix   
################################################
rand_distances=round(runif(nb_metabolites*nb_metabolites, 0, 1),2) #genration of number between 1 and 12
#conversion into matrix
data_test = matrix(rand_distances,nb_metabolites, nb_metabolites)
labels=paste("met",seq(1:nb_metabolites))
rownames(data_test)=labels
colnames(data_test)=labels
#conversion into triangular matrix
#data_test[upper.tri(data_test)] = 0


data=read.table("matrix.txt",header=F,sep="\t",dec=".",row.names=1)
colnames(data)=rownames(data)

#conversion into distance
distance_matrix=as.dist(data)
#plotcolors(dmat.color(distance_matrix))

#library(vegan)
#distance_matrix=vegdist(data,"jaccard")
################################
#          Clustering
################################
classif=hclust(distance_matrix,method="ward.D2")
#classif=hclust(distance_matrix,method="complete")
#classif=hclust(distance_matrix,method="median")

#par(mfrow=c(2,2)); plot(dendro);plot(dendro2);plot(dendro3);plot(dendro4);par(mfrow=c(1,1))

#automaticly ordering by clusters
classif = reorder.hclust(classif, data)

################################
#          Cophenetic
################################

coph_matrix=cophenetic(classif)
cor_coph=cor(distance_matrix,coph_matrix)
#print(paste("% of explained variance:",round(cor_coph^2,3)))
plotCohenetic=function(){
  plot(distance_matrix, coph_matrix, pch=19,col="red",cex=0.5,axes=F, cex.lab=1.2,cex.main=1.1,font.lab=font_size,xlim=c(0,max(distance_matrix)), ylim=c(0,max(coph_matrix)),xlab="Distance between metabolites",ylab="Cophenetic distance", asp=1, main=paste("Cophenetic correlation: ",round(cor_coph,3)))
  axis(2, seq(0.0,max(coph_matrix),interval),lwd=2,font.axis=font_size,cex.axis=0.8)
  axis(1, seq(0,max(distance_matrix),interval),lwd=2,font.axis=font_size,cex.axis=0.8)
  #lines(lowess(distance_matrix, coph_matrix), col="red",lwd=font_size)
  abline(0,1,col="grey",lwd=font_size,lty=2)
}

plotCohenetic()

################################
#          Height difference
################################

#Show differences between nodes levels (distance between clusters)
getHeightDifference=function(){
  height_classif=as.matrix(classif$height)
  height_diff=matrix(0, length(height_classif), 1)
  for (i in 2:(length(height_classif))){
    height_diff[i,]=height_classif[i,]-height_classif[i-1,]
  }
  rownames(height_diff)=c((length(height_classif)+1):2)
  return(height_diff)
}
height_diff=getHeightDifference()
#matrix_height=cbind(height_classif,height_diff)
#colnames(matrix_height)=c("node height","difference")

getRankedHeight=function(){
  ranked_height_diff=data.frame(height_diff)
  ranked_height_diff=ranked_height_diff[order(-ranked_height_diff), , drop = FALSE]
  rownames(ranked_height_diff)
  return(ranked_height_diff)
}
ranked_height_diff=getRankedHeight()

optimal_nb_clusters=rownames(ranked_height_diff)[1]
#optimal_nb_clusters=nrow(data)-(which.max(getHeightDifference()[-(nrow(data)-1)])-1)

getClusters= function(nb_clusters) {
  cutree(classif,nb_clusters)
}

clusters=getClusters(optimal_nb_clusters)
#clusters<-as.factor(cutree(dendro3,nb_clusters))

getClusterSizes=function(){
  cluster_sizes=table(clusters)
  cluster_sizes=data.frame(cluster_sizes)[,2]
  names(cluster_sizes)=paste("G",seq(1:optimal_nb_clusters),sep="")
  cluster_sizes=data.frame(cluster_sizes)
  names(cluster_sizes)="Effectif"
  return(cluster_sizes)
}
cluster_sizes=getClusterSizes()
print(cluster_size)

################################
#          Fusion levels
################################

getGroupContent=function(){
  
  result=matrix(0,nrow=(nrow(data)-1),ncol=2)
  
  setGrpContent=function(x,y){
    return (paste("(",x,",",y,")",sep=""))
  }
  
  getGrpNumber=function(x){
    strsplit(result[x,1],"G")[[1]][2] 
  }
  
  group_number=0
  for (i in 1:(nrow(data)-1)){
    element1=classif$merge[i,1]
    element2=classif$merge[i,2]
    if( element1 < 0 && element2 < 0){
      group_number=group_number+1
      group_name=paste("G",group_number,sep="") #ex. G1 for the first
      ordered_singletons=sort(c(abs(element1),abs(element2))) #if element1=-2 and element2=-1, reorder them
      group_content=setGrpContent(ordered_singletons[1],ordered_singletons[2]) #ex. (1,2)
    }else if( element1 < 0 && element2 > 0){
      group_name=result[element2,1]
      group_content=setGrpContent(get(group_name),abs(element1))
    }else if( element1 > 0 && element2 > 0){
      #the first group is printedfirst
      ordered_singletons=sort(c(getGrpNumber(element1),getGrpNumber(element2)))
      group_name=paste("G",ordered_singletons[1],sep="")
      group_content=setGrpContent(get(group_name),get(paste("G",ordered_singletons[2],sep="")))
    }
    assign(group_name,group_content) #create a variable named "group_name" with the content of "group_content"
    result[i,1]=group_name
    result[i,2]=group_content
    #print(group_content)
  }
  
  return (result)
}

#Plot fusion graph
plot_fusion_levels = function() {
  par(margin);#x11()
  subset_height=classif$height[(nrow(data)-max_cluster):(nrow(data)-1)]
  plot(classif$height, nrow(data):2, type="S",main="Distance before each fusion level",cex.main=2,cex.lab=1.5,lwd=font_size,xlim=c(min(subset_height),max(subset_height)),ylim=c(2,max_cluster),font.lab=2,ylab="Number of clusters", xlab="Node height", sub="(in red, difference in height with the previous fusion level)",col="grey", axes=F)
  axis(2, seq(2,max_cluster),lwd=2,font.axis=font_size)
  axis(1,seq(round(min(subset_height)),round(max(subset_height))),lwd=2,font.axis=font_size)
  result=matrix(0,nrow=(nrow(data)-1),ncol=2) #inialize an array for the result
  result=getGroupContent()
  text(x=classif$height[-1], y=(nrow(data)-1):2, labels=round(height_diff[-1],3), pos=3,col="red", cex=0.8)
  #catch_printing=identify(x=classif$height[-1], y=(nrow(data)-1):2,labels=paste(round(height_diff[-1],digits=2), result[-(nrow(data)-1),2], sep="\n"),col="red", cex=0.8,plot=T)
}

plot_fusion_levels()

################################
#          Dendrogram
################################

plotDendrogram=function(nb_clusters){
  par(margin);#x11()
  plot(classif, ylim=c(0,max(classif$height)),xlim=c(0,nrow(data)),hang=-1, cex.main=2, cex.lab=1.5,lwd=font_size,xlab="Metabolites", sub="",ylab="Distance before each fusion",main="Dendrogram",font.lab=font_size,axes=F)
  axis(2, seq(0,max(classif$height)),lwd=font_size,font.axis=font_size,cex.axis=0.8)
  #abline(h=seq(0.0,max(classif$height),0.1), lty=3, col="grey")
  #abline(h=c(classif$height), lty=3, col="grey")
  #projection of the clusters
  rect.hclust(classif, k=nb_clusters, border=rainbow(nb_clusters))
}

x11();plotDendrogram(as.numeric(optimal_nb_clusters))

################################
#          Silhouette
################################

plotAllSilhouette=function(max_cluster){
  asw <- numeric(max_cluster)
  for (k in 2:(max_cluster-1)) {
    sil <- silhouette(cutree(classif, k=k), distance_matrix)
    asw[k] <- summary(sil)$avg.width
  }
  
  k.best <- which.max(asw)
  # The plot is produced by function plot.silhouette {cluster}
  plot(1:nrow(data), asw, type="b", xlim=c(2,nrow(data)/3),ylim=c(0,max(asw)+0.1),col="grey",main="Silhouette-optimal number of clusters, Ward",xlab="Number of groups", ylab="Average silhouette width", axes=F)
  text(k.best,max(asw),paste("optimum",k.best,sep="\n \n"),col="red")
  axis(1, seq(2,nrow(data)/3),lwd=font_size,font.axis=font_size,cex.axis=0.8)
  axis(2, seq(0.0,(max(asw)+0.1),0.1),lwd=font_size,font.axis=font_size,cex.axis=0.8)
  points(k.best, max(asw), pch=21, col="red", cex=1)
  cat("","Silhouette-optimal number of clusters k =", k.best, "\n","with an average silhouette width of", round(max(asw),4), "\n")
  abline(v=k.best,lty=2,col="red")
}

plotSmallClusterSilhouette=function(max_cluster){
  par(mfrow=c(2,2))
  for (k in 3:max_cluster) {
    sil <- silhouette(cutree(classif, k=k), distance_matrix)
    silo <- sortSilhouette(sil)
    rownames(silo) <- row.names(data)[attr(silo,"iOrd")]
    plot(silo, main=paste("Silhouette plot for",k ,"groups"),cex.names=0.8, col=silo[,1]+1, nmax.lab=100)
  }
  par(mfrow=c(1,1))
}

plotAllSilhouette(max_cluster)
plotSmallClusterSilhouette(max_cluster)

################################
#          Inertie inter-
################################

getInertieInter=function(classif,k) {
  sum_inertia <- 0;
  element=length(classif$label)-1 ; imax <- k-1
  for (i in 1:imax) {
    sum_inertia=sum_inertia+classif$height[element]
    element=element-1
  }
  inertia=1000*sum_inertia/sum(classif$height) ; inertia=round(inertia)/10
  return(inertia)
}

plotInertiaInter=function(max_cluster){
  inertia=vector(mode="numeric",max_cluster-1)
  for (k in 2:max_cluster){
    inertia[k-1] = getInertieInter(classif,k)
  }
  k.best=which.max(inertia)
  plot(inertia,pch=19,type="b",ylim=c(0,(max(inertia)+5)),cex.lab=1.2,col="grey",axes=F,xlab="Nb. of cluster", ylab="Inertia inter-cluster")
  axis(1, seq(2,(nrow(data)/4)),lwd=font_size,font.axis=font_size,cex.axis=0.8)
  axis(2, seq(0,max(inertia)+5,10),lwd=font_size,font.axis=font_size,cex.axis=0.8)
  text(k.best,max(inertia),paste("optimum",k.best,sep="\n \n"),col="red")
  points(k.best, max(inertia), pch=21, col="red", cex=1)
  abline(v=k.best,lty=2,col="red")
  cat("","Silhouette-optimal number of clusters k =", k.best, "\n","with an inertia of", round(max(inertia),4), "\n")
}

plotInertiaInter(max_cluster)