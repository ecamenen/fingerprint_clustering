#clean all objects
rm(list=ls())
setwd("~/bin/fingerprint_clustering")

#global variables
nb_cluster=2
font_size=3
nb_metabolites=9
max_cluster=6
margin=par(mar=c(5, 4, 4, 2) + 1.1)
interval=1
typeClassif=3
#max_cluster=nb_metabolites/1.5
#choix du niveau de coupure
default_graph_par=par()

library(cluster)
library(gclus)
library(ade4)

#Pseudo-random settings: 
#set.seed(1)
#milisec * PID
set.seed(as.numeric(format(Sys.time(), "%OS2"))*100 * Sys.getpid())

################################################
#     Data test: random distance matrix   
################################################
rand_distances=ceiling(runif(nb_metabolites*nb_metabolites, 0, 11)) #generation of number between 1 and 11
#conversion into matrix
data_test = matrix(rand_distances,nb_metabolites, nb_metabolites)
labels=paste("met",seq(1:nb_metabolites))
rownames(data_test)=labels
colnames(data_test)=labels
data_test[cbind(1:nrow(data_test),1:nrow(data_test))] = 0
#conversion into symmetric matrix
data_test[lower.tri(data_test)] = t(data_test)[lower.tri(data_test)]
#write(data_test,"data_test.tsv",ncolumns=nb_metabolites,sep="\t")

#data=read.table("data_test.tsv",header=F,sep="\t",dec=".")
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

if (typeClassif==1){
  classif=pam(data,nb_cluster,diss=F,stand=F)
}else if (typeClassif==2){
  classif=kmeans(data,centers=nb_cluster,nstart=100)
}else if (typeClassif==3){
  classif=hclust(distance_matrix,method="ward.D2")
}else if (typeClassif==4){
  classif=hclust(distance_matrix,method="complete")
}else if (typeClassif==5){
  classif=hclust(distance_matrix,method="median")
}

#par(mfrow=c(2,2)); plot(dendro);plot(dendro2);plot(dendro3);plot(dendro4);par(mfrow=c(1,1))

#automaticly ordering by clusters
if(typeClassif>2) classif = reorder.hclust(classif, data)

colPers = colorRampPalette(c(rgb(0.6,0.1,0.5,1), rgb(1,0,0,1), rgb(0.9,0.6,0,1), rgb(0.1,0.6,0.3,1), rgb(0.1,0.6,0.5,1), rgb(0,0,1,1)), alpha = TRUE)

################################
#          Cophenetic
################################

plotCohenetic=function(){
  coph_matrix=cophenetic(classif)
  cor_coph=cor(distance_matrix,coph_matrix)
  #print(paste("% of explained variance:",round(cor_coph^2,3)))
  #x11()
  pdf("shepard_graph.pdf")
  plot(distance_matrix, coph_matrix, pch=19,col="red",cex=1,axes=F, cex.lab=1.5,cex.main=2,font.lab=font_size,xlim=c(0,max(distance_matrix)), ylim=c(0,max(coph_matrix)),xlab="Distance between metabolites",ylab="Cophenetic distance", asp=1, main=paste("Cophenetic correlation: ",round(cor_coph,3)))
  axis(2, seq(0.0,max(coph_matrix),interval),lwd=2,font.axis=font_size,cex.axis=0.8)
  axis(1, seq(0,max(distance_matrix),interval),lwd=2,font.axis=font_size,cex.axis=0.8)
  #lines(lowess(distance_matrix, coph_matrix), col="red",lwd=font_size)
  abline(0,1,col="grey",lwd=font_size,lty=2)
  dev.off()
}

if(typeClassif>2) plotCohenetic()

################################
#          Inertie inter-
################################

getCumulatedInertiaInter = function(classif, nb_cluster) {
  if (typeClassif==2) {
    classif=kmeans(data,centers=nb_cluster,nstart=100)
    return (round(classif$betweenss/classif$totss,3)*100)
  }else{
    sum_inertia = 0
    element = length(classif$label) - 1
    imax = nb_cluster - 1
    for (i in 1:imax) {
      sum_inertia = sum_inertia + classif$height[element]
      element = element-1
    }
    inertia = 100 * sum_inertia / sum(classif$height)
    inertia = round(inertia,2)
    return(inertia)
  }
}

getInertiaInter = function(classif, max_cluster) {
  inertia=vector(mode="numeric",max_cluster)
  if (typeClassif==2) {
    for (k in 2:(max_cluster+1)){
      inertia[k-1] = (kmeans(data,centers=k,nstart=100))$betweenss
    }
    return (inertia)
  }else{
    max_lenght=length(classif$height)
    return (classif$height[(max_lenght-max_cluster+1):max_lenght])
  }
}

getCumulatedInertiaInterPerCluster=function(max_cluster){
  inertia=vector(mode="numeric",max_cluster)
  for (k in 2:(max_cluster+1)){
    inertia[k-1] = getCumulatedInertiaInter(classif,k)
  }
  return (inertia)
}

plotCumulatedInertiaInter=function(max_cluster){
  #x11()
  inertia=getCumulatedInertiaInterPerCluster(max_cluster)
  k.best=which.max(inertia)
  pdf("cumulated_between.pdf")
  plot(inertia,type="b",ylim=c(0,(max(inertia)+5)),cex.lab=1.2,col="grey",axes=F,xlab="Nb. of cluster", ylab="Cumulated inertia inter-cluster")
  axis(1, seq(2,max_cluster),lwd=font_size,font.axis=font_size,cex.axis=0.8)
  axis(2, seq(0,max(inertia)+5,10),lwd=font_size,font.axis=font_size,cex.axis=0.8)
  text(k.best,max(inertia),paste("",round(max(inertia),4),sep="\n \n"),col="red",pos=2)
  points(k.best, max(inertia), pch=21, col="red", cex=1)
  abline(v=k.best,lty=2,col="red")
  cat("","Between-inertia optimal number of clusters k =", k.best, "\n","with a cumulated inter-inertia of", round(max(inertia),4), "\n")
  dev.off()
}

if(typeClassif>1) plotCumulatedInertiaInter(max_cluster)

################################
#          Height difference
################################

#Show distance differences between nodes levels (distance between clusters)
getHeightDifference=function(){
  height_classif=getInertiaInter(classif,max_cluster)
  height_diff=matrix(0, length(height_classif), 1)
  for (i in 2:(length(height_classif))){
    height_diff[i,]=height_classif[i]-height_classif[i-1]
  }
  rownames(height_diff)=c((length(height_classif)+1):2)
  return(height_diff[-1])
}

height_diff=getHeightDifference()

#matrix_height=cbind(height_classif,height_diff)
#colnames(matrix_height)=c("node height","difference")

getRankedHeight=function(){
  ranked_height_diff=data.frame(getHeightDifference())
  ranked_height_diff=ranked_height_diff[order(-ranked_height_diff), , drop = FALSE]
  rownames(ranked_height_diff)=max_cluster-as.numeric(rownames(ranked_height_diff))+1
  return(ranked_height_diff)
}

printTableInertia=function(){
  ranked_height_diff=getRankedHeight()
  matrix_output1=cbind(getInertiaInter(classif,max_cluster)[-1],getHeightDifference())
  for(i in 1:ncol(matrix_output1)) {matrix_output1[,i] = rev(matrix_output1[,i])}
  matrix_output1=cbind(matrix_output1,getCumulatedInertiaInterPerCluster(max_cluster-1))
  rownames(matrix_output1)=seq(2,max_cluster)
  colnames(matrix_output1)=c("Branch height", "Differences","Cumulated inertia")
  matrix_output1=round(matrix_output1,2)
  return (matrix_output1)
}

writeTsv=function(x,f){
  print(x)
  x[is.na(x)] <-""
  output=rbind(c("",colnames(x)), cbind(rownames(x),x))
  write(t(output),file=f,ncolumns=ncol(output),sep="\t")
}

summary_between=printTableInertia()
writeTsv(summary_between,"summary_between.tsv")


optimal_nb_clusters=as.numeric(rownames(getRankedHeight())[1])
#optimal_nb_clusters=nrow(data)-(which.max(getHeightDifference()[-(nrow(data)-1)])-1)

getClusters= function(nb_clusters) {
  if(typeClassif> 2) cutree(classif,nb_clusters)
  else if (typeClassif ==2) (kmeans(data,centers=nb_clusters,nstart=100))$cluster
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
writeTsv(cluster_sizes,"cluster_sizes.tsv")

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
  subset_height=getInertiaInter(classif,max_cluster)
  pdf("fusion_levels.pdf")
  plot(2:max_cluster, rev(subset_height[-1]), font.lab=3,type="b",ylab="Distance inter-group",cex.main=2,cex.lab=1.5,lwd=font_size,ylim=c(min(subset_height),max(subset_height)),xlim=c(2,max_cluster),xlab="Number of groups", main="Fusion levels plot",col="grey", axes=F)
  legend("top",legend="(in red, distance difference with the previous fusion level)",bty="n")
  axis(1, seq(2,max_cluster),lwd=2)
  if(typeClassif==2){ 
    interval=100
  }else{ 
    interval=1
    }
  axis(2,seq(round(min(subset_height)),round(max(subset_height)),by=interval),lwd=2)
  text(y=rev(subset_height[-1]), x=2:max_cluster,labels=rev(round(height_diff,2)), cex=1.2,pos=4,col="red")
  points(optimal_nb_clusters, subset_height[max_cluster+2-optimal_nb_clusters], pch=19, col="red", cex=1.8)
  abline(v=optimal_nb_clusters,col="red",lty=2,lwd=2)
  #catch_printing=identify(x=classif$height[-1], y=(nrow(data)-1):2,labels=paste(round(height_diff[-1],digits=2), result[-(nrow(data)-1),2], sep="\n"),col="red", cex=0.8,plot=T)
  dev.off()
  }

if(typeClassif>1) plot_fusion_levels()

################################
#            PDIS
################################
centreduire <- function(T) { 
  N <- nrow(T) ; T1 <- scale(T); 
  return(T1*sqrt(N/(N-1))) 
} 

#Partitionner la classification
#Sortie: partitionnement contenant k clusters
cutTree=function(classif,k){
  if(typeClassif>2) cutree(classif, k=k)
  else if (typeClassif==1) pam(data,k,diss=F,stand=F)
  else if (typeClassif==2) (kmeans(data,centers=k,nstart=100))$cluster
}

# Pouvoir discriminant des variables
# Parametres :	table des donnees,	classification hierarchique,
#		nombre de classes
# Sortie : les carres des distances (PDIS)
pdis = function(T1,H,k) {
  T = centreduire(T1)
  N = nrow(T) ; M = ncol(T)
  C = cutTree(H,k)
  ctr = matrix(0,nrow=k,ncol=M)
  for (i in 1:N) {
    cli = C[i]
    for (j in 1:M) ctr[cli,j] = ctr[cli,j] + T[i,j]
  }
  r = vector(mode="numeric",M)
  for (i in 1:k)
    for (j in 1:M) ctr[i,j] = ctr[i,j]^2/(N*length(C[C==i]))
    for (i in 1:M) r[i] = sum(ctr[,i])
    return(round(1000*r)/10)
}

pdis_classif=matrix(0,max_cluster-1,ncol(data))
colnames(pdis_classif)=colnames(data)
rownames(pdis_classif)=seq(2,max_cluster)

for (k in 2:max_cluster){
  res = pdis(data,classif,k)
  for(i in 1:length(res)){
    pdis_classif[k-1,i]=round(res[i],2)
  }
}

writeTsv(pdis_classif,"discriminant_power.tsv")

################################
#            RHO2
################################

# Distance**2 des classes au centre du nuage
# Parametres :	table des donnees,
#		classement hierarchique,
#		nombre de classes
# Sortie : les carres des distances (souvent notes RHO2)
rho2 <- function(T1,H,k) {
  T <- centreduire(T1);
  N <- nrow(T) ; M <- ncol(T);
  C <- cutTree(H,k);
  cdg <- matrix(data=0,nrow=k,ncol=M);
  for (i in 1:N) {
    cli <- C[i];
    for (j in 1:M) cdg[cli,j] <- cdg[cli,j] + T[i,j];
  };
  for (i in 1:k)
    for (j in 1:M) cdg[i,j] <- cdg[i,j]/length(C[C==i]);
    r <- vector(mode="numeric",k);
    for (i in 1:k) r[i] <- sum(cdg[i,]^2);
    return(r)
}

excentricity=matrix(0,max_cluster-1,max_cluster)
rownames(excentricity)=seq(2,max_cluster)
colnames(excentricity)=paste("G",seq(1,max_cluster),sep="")
for (k in 2:max_cluster){
  res=rho2(data,classif,k);print(res)
  for(i in 1:length(res)){
    excentricity[k-1,i]=round(res[i],2)
  }
}
excentricity[excentricity==0] <-NA

writeTsv(excentricity,"excentricity.tsv")
################################
#            CTR
################################

# Contribution relative des classes a inertie du nuage
# Parametres :	table des donn?es,
#		classfication hierarchique,
#		nombre de classes
# Sortie : les contributions (souvent notees CTR)
ctrng <- function(T1,H,k) {
  T <- centreduire(T1)
  N <- nrow(T) ; M <- ncol(T)
  C <- cutTree(H,k)
  ctr <- matrix(0,nrow=k,ncol=M)
  for (i in 1:N) {
    cli <- C[i]
    for (j in 1:M) ctr[cli,j] <- ctr[cli,j] + T[i,j]
  }
  for (i in 1:k)
    for (j in 1:M) ctr[i,j] <- ctr[i,j]^2/(N*length(C[C==i]))
    ctrframe <- as.data.frame(ctr)
    colnames(ctrframe) <- colnames(T1)
    return(round(1000*ctrframe)/10)
}

relative_ctr = ctrng(data,classif,3)
writeTsv(relative_ctr,"relative_ctr.tsv")

################################
#          Silhouette
################################

plotAllSilhouette=function(max_cluster){
  
  asw <- numeric(max_cluster-1)
  for (k in 2:(max_cluster-1)) {
    si = silhouette(cutTree(classif, k=k), distance_matrix)
    asw[k] <- summary(si)$avg.width
  }
  
  x11()
  k.best <- which.max(asw)
  # The plot is produced by function plot.silhouette {cluster}
  plot(1:(max_cluster-1), asw, type="b", lwd=2,cex=1.2,font.lab=3,xlim=c(2,(max_cluster-1)),cex.main=2, cex.lab=1.5,ylim=c(0,max(asw)+0.1),col="grey",main="Silhouette plot for k groups",xlab="Number of groups", ylab="Average silhouette width", axes=F)
  text(k.best,max(asw),round(max(asw),3),col="red", pos=4,cex=1.2)
  axis(1, seq(2,(max_cluster-1)),lwd=font_size,font.axis=font_size)
  axis(2, seq(0.0,(max(asw)+0.1),0.1),lwd=font_size,font.axis=font_size)
  points(k.best, max(asw), pch=19, col="red", cex=1.5)
  cat("","Silhouette-optimal number of clusters k =", k.best, "\n","with an average silhouette width of", round(max(asw),4), "\n")
  abline(v=k.best,lty=2,col="red",lwd=2)
  return (k.best)
}

plotSmallClusterSilhouette=function(max_cluster){
  x11()
  par(mfrow=c(2,2))
  for (k in 2:4) {
    sil <- silhouette(cutTree(classif, k=k), distance_matrix)
    silo <- sortSilhouette(sil)
    rownames(silo) <- row.names(data)[attr(silo,"iOrd")]
    plot(silo, main=paste("Silhouette plot for",k ,"groups"),cex.names=0.8, col=silo[,1]+1, nmax.lab=100)
  }
  par(default_graph_par)
}

optimal_nb_clusters=plotAllSilhouette(max_cluster+1)
plotSmallClusterSilhouette(max_cluster)

################################
#          Dendrogram
################################

plotDendrogram=function(nb_clusters){
  par(margin);#x11()
  pdf("dendrogram.pdf")
  plot(classif, ylim=c(0,max(classif$height)),xlim=c(0,nrow(data)),hang=-1, cex.main=2, cex.lab=1.5,lwd=font_size,xlab="Metabolites", sub="",ylab="Distance inter-group",main="Dendrogram",font.lab=font_size,axes=F)
  axis(2, seq(0,max(classif$height)),lwd=font_size,font.axis=font_size,cex.axis=0.8)
  #abline(h=seq(0.0,max(classif$height),0.1), lty=3, col="grey")
  #abline(h=c(classif$height), lty=3, lwd=2, col="grey")
  #projection of the clusters
  rect.hclust(classif, k=nb_clusters, border=colPers(nb_clusters))
  dev.off()
}

if(typeClassif>2) plotDendrogram(as.numeric(optimal_nb_clusters))

################################
#            PCA
################################

acp=dudi.pca(data,scannf=F)
#x11()

pdf("pca.pdf")
par(mar=c(0,0,0,0))
clusters=getClusters(optimal_nb_clusters)
title=paste("Cumulated inertia: ",round((acp$eig[1]+acp$eig[2])/sum(acp$eig),4)*100,"%",sep="")
#s.class(addaxes=F,acp$li,sub=title,possub="topright",csub=1.5,as.factor(getClusters(optimal_nb_clusters)),grid=F,col=colPers(optimal_nb_clusters))
s.class(addaxes=F,acp$li,ylim=c(min(acp$li[,2])-3,max(acp$li[,2])+3),xlim=c(min(acp$li[,1])-3,max(acp$li[,1])+3),sub=title,possub="topright",csub=1.5,as.factor(clusters),grid=F,col=colPers(optimal_nb_clusters))
abline(h=0,v=0,lty=2,lwd=2,col="grey")
for (i in 1:optimal_nb_clusters){
  clusters[clusters==i] = colPers(optimal_nb_clusters)[i]
}
text(x=acp$li[,1],y=acp$li[,2],labels=rownames(acp$li),col=clusters,cex=0.5)
#s.label(acp$li,add.plot = T,boxes=F,clabel=0.8,col="red")
dev.off()


#heatmap.2(as.matrix(distance_matrix),distfun = function(x) dist(x,method = 'euclidean'),hclustfun = function(x) hclust(x,method = 'complete'),key=FALSE, density.info="none",lhei = c(0.05,0.95) ,cexRow=2,cexCol=2,margins=c(20,26),trace="none",srtCol=45,dendrogram="row")
#plot(data[,1],data[,2],type="p",pch="+",xlab="Axis 1", ylab="Axis 2",col=rainbow(optimal_nb_clusters)[clusters])
#points(classif$medoids[,1],classif$medoids[,2], cex=1.5,pch=16,col=c("red","blue","green")[1:3])
#s.label(data,boxes=F,clabel=0.8,neig = F,add.plot = T)
#points(data[,1],data[,2],col=c("lightcoral","skyblue","greenyellow")[classif$clustering])

#classif$cluster[order(-classif$cluster,decreasing=T)]