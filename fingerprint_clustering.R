#clean all objects
rm(list=ls())
setwd("~/bin/fingerprint_clustering")

#global variables
font_size=2
nb_metabolites=9
max_cluster=6
margin=par(mar=c(5, 4, 4, 2) + 1.1)
interval=1
typeClassif=3
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

if (typeClassif==1){
  classif=pam(data,3,diss=F,stand=F)
}else if (typeClassif==2){
  classif=kmeans(data,centers=4,nstart=100)
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

if(typeClassif>2) plotCohenetic()

################################
#          Inertie inter-
################################

getCumulatedInertieInter = function(classif, nb_cluster) {
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

getCumulatedInertiaInterPerCluster=function(max_cluster){
  inertia=vector(mode="numeric",max_cluster)
  for (k in 2:(max_cluster+1)){
    inertia[k-1] = getCumulatedInertieInter(classif,k)
  }
  return (inertia)
}

plotCumulatedInertiaInter=function(max_cluster){
  inertia=getCumulatedInertiaInterPerCluster(max_cluster)
  k.best=which.max(inertia)
  plot(inertia,type="b",ylim=c(0,(max(inertia)+5)),cex.lab=1.2,col="grey",axes=F,xlab="Nb. of cluster", ylab="Cumulated inertia inter-cluster")
  axis(1, seq(2,max_cluster),lwd=font_size,font.axis=font_size,cex.axis=0.8)
  axis(2, seq(0,max(inertia)+5,10),lwd=font_size,font.axis=font_size,cex.axis=0.8)
  text(k.best,max(inertia),paste("",round(max(inertia),4),sep="\n \n"),col="red",pos=2)
  points(k.best, max(inertia), pch=21, col="red", cex=1)
  abline(v=k.best,lty=2,col="red")
  cat("","Silhouette-optimal number of clusters k =", k.best, "\n","with a cumulated inter-inertia of", round(max(inertia),4), "\n")
}

plotCumulatedInertiaInter(max_cluster)

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
  ranked_height_diff=data.frame(getHeightDifference())
  ranked_height_diff=ranked_height_diff[order(-ranked_height_diff), , drop = FALSE]
  rownames(ranked_height_diff)
  return(ranked_height_diff)
}

printTableInertia=function(){
  ranked_height_diff=getRankedHeight()
  matrix_output1=cbind(classif$height,getHeightDifference())
  for(i in 1:ncol(matrix_output1)) {matrix_output1[,i] = rev(matrix_output1[,i])}
  matrix_output1=cbind(matrix_output1,getCumulatedInertiaInterPerCluster(nrow(data)-1))
  rownames(matrix_output1)=seq(2,(nrow(data)))
  colnames(matrix_output1)=c("Branch height", "Differences","Cumulated inertia")
  matrix_output1=round(matrix_output1,2)
  print (matrix_output1)
}

printTableInertia()

optimal_nb_clusters=as.numeric(rownames(getRankedHeight())[1])
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
print(cluster_sizes)

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
  plot(classif$height, nrow(data):2, type="b",main="Distance between each clustering event (=fusion)",cex.main=2,cex.lab=1.5,lwd=font_size,xlim=c(min(subset_height),max(subset_height)),ylim=c(2,max_cluster),font.lab=2,ylab="Number of clusters", xlab="Node height", sub="(in red, difference in height with the previous fusion level)",col="grey", axes=F)
  axis(2, seq(2,max_cluster),lwd=2,font.axis=font_size)
  axis(1,seq(round(min(subset_height)),round(max(subset_height))),lwd=2,font.axis=font_size)
  result=matrix(0,nrow=(nrow(data)-1),ncol=2) #inialize an array for the result
  result=getGroupContent()
  text(x=classif$height[-1], y=(nrow(data)-1):2, labels=round(height_diff[-1],3), pos=3,col="red", cex=0.8)
  #catch_printing=identify(x=classif$height[-1], y=(nrow(data)-1):2,labels=paste(round(height_diff[-1],digits=2), result[-(nrow(data)-1),2], sep="\n"),col="red", cex=0.8,plot=T)
}

plot_fusion_levels()

################################
#          COR
################################
centreduire <- function(T) {
  N <- nrow(T) ; T1 <- scale(T);
  return(T1*sqrt(N/(N-1)))
  }

# Calcule les cdg  centres reduits des classes.
# Paramatres :	table des donnees,
#		hierarchie,
#		nombre de classes
# Sortie : les coordonnees des centres de gravite des classes
cdgcl <- function(T1,H,k) {
  T <- centreduire(T1);
  N <- nrow(T) ; M <- ncol(T);
  C <- cutree(H,k);
  cdg <- matrix(data=0,nrow=k,ncol=M);
  for (i in 1:N) {
    cli <- C[i];
    for (j in 1:M) cdg[cli,j] <- cdg[cli,j] + T[i,j];
  };
  for (i in 1:k)
    for (j in 1:M) cdg[i,j] <- cdg[i,j]/length(C[C==i]);
    cdgframe <- as.data.frame(cdg);
    names(cdgframe) <- names(T1);
    return(cdgframe)}

# Contribution relative des variables a l'eloignement des classes
# Parametres :	table des donnees,
#		classement hierarchique,
#		nombre de classes
# Sortie : les contributions relatives (isouvent notees COR)

ctrcl <- function(T,H,k) {
  N <- nrow(T) ; M <- ncol(T);
  cdg <- cdgcl(T,H,k);
  ctr <- cdg;
  for (i in 1:k) {
    s2 <- sum(cdg[i,]^2);
    for (j in 1:M) ctr[i,j] <- cdg[i,j]*abs(cdg[i,j])/s2;
  };
  ctr <- round(1000*ctr);
  return(ctr/10)}

correlation_var=ctrcl(data,classif,optimal_nb_clusters)
correlation_var

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
  C <- cutree(H,k);
  cdg <- matrix(data=0,nrow=k,ncol=M);
  for (i in 1:N) {
    cli <- C[i];
    for (j in 1:M) cdg[cli,j] <- cdg[cli,j] + T[i,j];
  };
  for (i in 1:k)
    for (j in 1:M) cdg[i,j] <- cdg[i,j]/length(C[C==i]);
    r <- vector(mode="numeric",k);
    for (i in 1:k) r[i] <- sum(cdg[i,]^2);
    return(r)}

excentricity=matrix(0,max_cluster-1,max_cluster-1)
rownames(excentricity)=seq(2,max_cluster)
colnames(excentricity)=paste("G",seq(1,max_cluster-1),sep="")
for (k in 2:max_cluster){
  res=rho2(data,classif,k);print(res)
  for(i in 1:length(res)){
    excentricity[k-1,i]=round(res[i],2)
  }
}
#excentricity[excentricity==0] <-" "
excentricity[excentricity==0] <-NA
excentricity

################################
#            CTR
################################

# Contribution relative des classes a inertie du nuage
# Parametres :	table des donn?es,
#		classfication hierarchique,
#		nombre de classes
# Sortie : les contributions (souvent notees CTR)
ctrng <- function(T1,H,k) {
  T <- centreduire(T1);
  N <- nrow(T) ; M <- ncol(T);
  C <- cutree(H,k);
  ctr <- matrix(0,nrow=k,ncol=M);
  for (i in 1:N) {
    cli <- C[i];
    for (j in 1:M) ctr[cli,j] <- ctr[cli,j] + T[i,j];
  };
  for (i in 1:k)
    for (j in 1:M) ctr[i,j] <- ctr[i,j]^2/(N*length(C[C==i]));
    ctrframe <- as.data.frame(ctr);
    names(ctrframe) <- names(T1);
    return(round(1000*ctrframe)/10)}

relative_ctr = ctrng(data,classif,optimal_nb_clusters)
relative_ctr

################################
#            CTR
################################

# Pouvoir discriminant des variables
# Parametres :	table des donnees,
#		tableau hierarchique,
#		nombre de classes
# Sortie : les carres des distances 

pdis <- function(T1,H,k) {
  T <- centreduire(T1);
  N <- nrow(T) ; M <- ncol(T);
  C <- cutree(H,k);
  ctr <- matrix(data=0,nrow=k,ncol=M);
  for (i in 1:N) {
    cli <- C[i];
    for (j in 1:M) ctr[cli,j] <- ctr[cli,j] + T[i,j];
  };
  r <- vector(mode="numeric",M);
  for (i in 1:k)
    for (j in 1:M) ctr[i,j] <- ctr[i,j]^2/(N*length(C[C==i]));
    for (i in 1:M) r[i] <- sum(ctr[,i]) ;
    return(round(1000*r)/10)}

pvdiscriminant=matrix(0,max_cluster-1,nrow(data))
colnames(pvdiscriminant)=colnames(data)
rownames(pvdiscriminant)=seq(1,max_cluster-1)

for (k in 2:max_cluster){
  res=pdis(data,classif,k)
  for(i in 1:length(res)){
    pvdiscriminant[k-1,i]=round(res[i],2)
  }
}


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
  asw <- numeric(max_cluster-1)
  for (k in 2:(max_cluster-1)) {
    sil <- silhouette(cutree(classif, k=k), distance_matrix)
    asw[k] <- summary(sil)$avg.width
  }
  
  k.best <- which.max(asw)
  # The plot is produced by function plot.silhouette {cluster}
  plot(1:(max_cluster-1), asw, type="b", xlim=c(2,(max_cluster-1)),ylim=c(0,max(asw)+0.1),col="grey",main="Silhouette-optimal number of clusters, Ward",xlab="Number of groups", ylab="Average silhouette width", axes=F)
  text(k.best,max(asw),paste("optimum",k.best,sep=":"),col="red", pos=2)
  axis(1, seq(2,(max_cluster-1)),lwd=font_size,font.axis=font_size,cex.axis=0.8)
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

plotAllSilhouette(40)
plotSmallClusterSilhouette(max_cluster)

plot(data[,1],data[,2],type="p",pch="+",xlab="X1", ylab="X2",col=c("lightcoral","skyblue","greenyellow")[classif$clustering])
points(classif$medoids[,1],classif$medoids[,2], cex=1.5,pch=16,col=c("red","blue","green")[1:3])
s.label(data,boxes=F,clabel=0.8,neig = F,add.plot = T)
points(data[,1],data[,2],col=c("lightcoral","skyblue","greenyellow")[classif$clustering])

classif$cluster[order(-classif$cluster,decreasing=T)]