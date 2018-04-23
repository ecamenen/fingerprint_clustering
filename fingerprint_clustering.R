#clean all objects
rm(list=ls())
setwd("~/bin/fingerprint_clustering")

#global variables
#choix du niveau de coupure
nb_clusters=2
font_size=3
nb_metabolites=9
max_cluster=6
#margin=par(mar=c(5, 4, 4, 2) + 1.1)
Betweenval=1
typeClassif=4


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

#Output: a random distance matrix (symetric, with a diagonal of 0)
setRandomDataSet = function(){
  #generation of number between 1 and 11
  rand_distances = ceiling(runif(nb_metabolites * nb_metabolites, 0, 11)) 
  #conversion into matrix
  data_test = matrix(rand_distances, nb_metabolites, nb_metabolites)
  #label met1, met2,...
  labels=paste("met", seq(1:nb_metabolites))
  rownames(data_test) = labels
  colnames(data_test) = labels
  #conversion of diagonal into 0
  data_test[cbind(1:nrow(data_test), 1:nrow(data_test))] = 0
  #conversion into symmetric matrix
  data_test[lower.tri(data_test)] = t(data_test)[lower.tri(data_test)]
  return (data_test)
}

#data=setRandomDataSet()
data = read.table("matrix.txt", header=F, sep="\t", dec=".", row.names=1)
colnames(data) = rownames(data)

#conversion into distance
#distance_matrix=as.dist(data)
getDistance = function(x) dist(x, method = "euclidian")
#library(vegan)
#distance_matrix=vegdist(data,"jaccard")
################################
#          Clustering
################################

#Inputs:
# t: number of type of classification
# d: data (or distance matrix for hierarchic)
#Ouput: Hierarchical classification
getCAH = function(d, t){
  if(t>2){
  #dis: distance matrix
  dis = dist(d, method = "euclidian")
    if (t==3) meth="ward.D2"
    else if (t==4) meth="complete"
    else if (t==5) meth="average"
    else if (t==6) meth="mcquitty"
    #cah: classification hierarchic ascending
    cah = hclust(dis, method=meth)
  #automaticly ordering by clusters
  return (reorder.hclust(cah, d))
  }
  #TODO: exit if 0 < t < 6
}

classif = getCAH(data, typeClassif)

#Inputs: 
# t: number of type of classification
# d: data (or distance matrix for hierarchic)
# k: number of clusterting
#Ouput: Non-hierarchical classification
getCNH = function(t, d, k){
  if (t==1) return (pam(d, k, diss=F, stand=F))
  else if (t==2) return (kmeans(d, centers=k, nstart=100))
}

#Usage: colPers(x), x a number of colours in output
#Gradient of color
colPers = colorRampPalette(c(rgb(0.6,0.1,0.5,1), rgb(1,0,0,1), rgb(0.9,0.6,0,1), rgb(0.1,0.6,0.3,1), rgb(0.1,0.6,0.5,1), rgb(0,0,1,1)), alpha = TRUE)

#Inputs: x : a matrix
#filename of the saved file
#Prints the matrix, save the matrix
writeTsv = function(x,f){
  print(x)
  output=as.matrix(rbind(c("", colnames(x)), cbind(rownames(x),x)))
  output[is.na(output)] = ""
  write(t(output), file=f, ncolumns=ncol(output), sep="\t")
  #write.table(x, f, na = "",col.names = colnames(x),row.names = rownames(x),append = F,sep = "\t")
}
############################################################
#          Cophenetic (dendrogram distance matrix)
############################################################

#Inputs:
# d : data
# cah : hierarchical classification
plotCohenetic=function(d,cah){
  library(scales)
  dis = getDistance(d)
  coph_matrix = cophenetic(cah)
  cor_coph = cor(dis, coph_matrix)
  print(paste("% of explained variance by Cophenetic:", round(cor_coph^2,3)))
  x11()
  par(mar=c(5.1,5.1,4.1,2.1))
  #pdf("shepard_graph.pdf")
  plot(dis, coph_matrix, pch=19,col=alpha("red",0.2), cex=1, axes=F, cex.lab=1.5, cex.main=2, font.lab=font_size, xlim=c(0,max(dis)), ylim=c(0,max(coph_matrix)), xlab="Distance between metabolites",ylab="Cophenetic distance", asp=1, main=paste("Cophenetic correlation: ",round(cor_coph,3)))
  axis(2, seq(0.0,max(coph_matrix),1), lwd=font_size, font.axis=3, cex.axis=0.8)
  axis(1, seq(0,max(dis),1), lwd=font_size, font.axis=3, cex.axis=0.8)
  abline(0,1, col="grey", lwd=font_size, lty=2)
  #dev.off()
}

if(typeClassif>2) plotCohenetic(data, classif)

##############################################
#          Between-group inertia
##############################################

# Inputs: 
# t: number of type of classification
# k: number of clusters
# c: hierarchical classification
# d: dataframe
#Ouput: cumulated between-group inertia of a classification
getCumulatedBetweenInertia = function(t, k, c=NULL, d=NULL) {
  if (t==2) {
    c = getKmeans(d, k)
    return (round(c$betweenss/c$totss,3) * 100)
  }else{
    sum_inertia = 0
    #begin with the last element
    element = length(c$label) - 1
    imax = k - 1
    for (i in 1:imax) {
      #sum until k elements (nb of partitionning)
      sum_inertia = sum_inertia + c$height[element]
      element = element-1
      #decremential loop
    }
    return(round(100 * sum_inertia / sum(c$height),2))
  }
}

# Inputs:
# t: number of type of classification
# n: maximum number of clusters
# c: hierarchical classification
# d: dataframe
# Output: between-group inertia for all clusters
getBetweenInertia = function(t, n, c=NULL, d=NULL) {
  inertia = vector(mode="numeric", n)
  if (t==2) {
    for (k in 2:(n+1)){
      inertia[k-1] = getKmeans(d, k)$betweenss
    }
    return (inertia)
  }else{
    max_lenght = length(c$height)
    #get all the elements until the max of partitionned fixed
    return (c$height[(max_lenght-n+1):max_lenght])
  }
}

# Inputs:
# t: number of type of classification
# n: maximum number of clusters
# c:  hierarchical classification
# d: data
getCumulatedBetweenInertiaPerCluster = function(t, n, c=NULL, d=NULL){
  inertia = vector(mode="numeric", n)
  for (k in 2:(n+1)){
    inertia[k-1] = getCumulatedBetweenInertia(t, k, c, d)
  }
  return (inertia)
}

# Inputs:
# t: number of type of classification
# n: maximum number of clusters
# k: number of clusters
# c: hierarchical classification
# d: data
plotCumulatedBetweenInertia = function(t, n, c=NULL, d=NULL){
  x11()
  par(mar=c(5.1,5.1,4.1,2.1))
  if(t==2) inertia = getCumulatedBetweenInertiaPerCluster(t, n, d=d)
  if(t>2) inertia = getCumulatedBetweenInertiaPerCluster(t, n, c)
  k.best = which.max(inertia)
  #pdf("cumulated_between.pdf")
  plot(inertia, type="b", lwd=font_size, font.lab=3, cex.lab=font_size/2, cex.main=font_size/1.5, col="grey", xlim=c(2, n), ylim=c(0,(max(inertia)+5)), axes=F, xlab="Nb. of clusters", ylab="Cumulated between-group inertia")
  axis(1, seq(2,n), lwd=font_size, font.axis=3, cex.axis=0.8)
  axis(2, seq(0,max(inertia)+5,10),lwd=font_size, font.axis=3, cex.axis=0.8)
  text(k.best, max(inertia), paste("",round(max(inertia),4),sep="\n \n"), col="red", pos=2, cex=font_size/2.5)
  points(k.best, max(inertia), pch=20, col="red", cex=font_size)
  abline(v=k.best, lty=2 ,col="red", lwd=font_size/1.5)
  cat("","between-inertia optimal number of clusters k =", k.best, "\n","with a cumulated between-inertia of", round(max(inertia),4), "\n")
  #dev.off()
}

plotCumulatedBetweenInertia(typeClassif, max_cluster, classif, data)

############################################################
#          Between-cluster differences
############################################################

# Inputs: 
# t: number of type of classification
# n: maximum number of clusters
# c: hierarchical classification
# d: data
getBetweenDifference = function(t, n, c=NULL, d=NULL){
  if(t==2) inertia = getBetweenInertia(t, n, d=d)
  if(t>2) inertia = getBetweenInertia(t, n, c)
  inertia_diff = matrix(0, length(inertia), 1)
  for (i in 2:(length(inertia))){
    inertia_diff[i,] = inertia[i] - inertia[i-1]
  }
  rownames(inertia_diff) = c((length(inertia) + 1):2)
  return(inertia_diff[-1])
}

getRankedInertia = function(t, n, c=NULL, d=NULL){
  ranked_inertia_diff = data.frame(getBetweenDifference(t, n, c, d))
  ranked_inertia_diff = ranked_inertia_diff[order(-ranked_inertia_diff), , drop = FALSE]
  rownames(ranked_inertia_diff) = n - as.numeric(rownames(ranked_inertia_diff)) + 1
  return(ranked_inertia_diff)
}

printTableInertia = function(t, n, c=NULL, d=NULL){
  table_inertia = cbind(getBetweenInertia(t, n, c, d)[-1], getBetweenDifference(t, n, c, d))
  #outputs are reversed comparatively to CumulatedBetween outputs
  for(i in 1:ncol(table_inertia)) {table_inertia[,i] = rev(table_inertia[,i])}
  table_inertia = cbind(table_inertia, getCumulatedBetweenInertiaPerCluster(t, n - 1, c, d))
  rownames(table_inertia) = seq(2, n)
  colnames(table_inertia) = c("Branch height", "Differences","Cumulated inertia")
  table_inertia = round(table_inertia, 2)
  return (table_inertia)
}
summary_between = printTableInertia(typeClassif, max_cluster, classif, data)
writeTsv(summary_between,"summary_between.tsv")

optimal_nb_clusters = as.numeric(rownames(getRankedInertia(typeClassif, max_cluster, classif, data))[1])

# Inputs: 
# t: number of type of classification
# k: number of clusters
# c: hierarchical classification
# d: data
# Output: partitionning contening k clusters
getClusters = function(t, k, c=NULL, d=NULL) {
  if (t > 2) cutree(c, k)
  else { 
    cnh = getCNH(t, d, k)
    if (t == 1) cnh$clustering
    else cnh$cluster
  }
}
clusters = getClusters(typeClassif, optimal_nb_clusters, classif, data)

#Input:
# cl: clusters
getClusterSizes = function(cl){
  cluster_sizes = table(cl)
  cluster_sizes = data.frame(cluster_sizes)[,2]
  cluster_sizes = data.frame(cluster_sizes)
  names(cluster_sizes) = "Effectif"
  return(cluster_sizes)
}
cluster_sizes=getClusterSizes(clusters)
writeTsv(cluster_sizes,"cluster_sizes.tsv")

################################
#          Fusion levels
################################

#Plot fusion graph
plot_fusion_levels = function(t, n, c=NULL, d=NULL) {
  x11()
  subset_height = getBetweenInertia(t, n, c, d)
  height_diff = getBetweenDifference(t, n, c, d)
  #pdf("fusion_levels.pdf")
  par(mar=c(5.1,5.1,5.1,2.1))
  plot(2:n, rev(subset_height[-1]), type="b", cex.lab=font_size/2, lwd=font_size, font.lab=3, ylim=c(min(subset_height),max(subset_height)), xlim=c(2,n), xlab="Nb. of clusters", ylab="Between-group inertia", col="grey", axes=F)
  title(main="Fusion levels", line=2,cex.main=font_size/1.5)
  mtext("(in red, inter-group differences with the previous clustering)", side=3, line=1)
  axis(1, seq(2,n), lwd=font_size, font.axis=3, cex.axis=0.8)
  if (typeClassif == 2) interval = 100
  else interval = 1
  axis(2, seq(round(min(subset_height)),round(max(subset_height)), by=interval), lwd=font_size, font.axis=3, cex.axis=0.8)
  text(y=rev(subset_height[-1]), x=2:max_cluster, labels=rev(round(height_diff,2)), cex=1.2, pos=4, col="red")
  points(optimal_nb_clusters, subset_height[max_cluster+2-optimal_nb_clusters], pch=19, col="red", cex=font_size/1.5)
  abline(v=optimal_nb_clusters, col="red", lty=2, lwd=font_size/1.5)
  #catch_printing=identify(x=classif$height[-1], y=(nrow(data)-1):2,labels=paste(round(height_diff[-1],digits=2), result[-(nrow(data)-1),2], sep="\n"),col="red", cex=0.8,plot=T)
  #dev.off()
  }

if(typeClassif > 1) plot_fusion_levels(typeClassif, max_cluster, classif, data)

################################
#            PDIS
################################
centreduire <- function(T) { 
  N <- nrow(T) ; T1 <- scale(T); 
  return(T1*sqrt(N/(N-1))) 
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

relative_ctr = ctrng(data,classif,optimal_nb_clusters)
writeTsv(relative_ctr,"relative_ctr.tsv")

################################
#          Silhouette
################################

plotAllSilhouette = function(t, n, c=NULL, d=NULL){
  
  mean_silhouette = numeric(n - 1)
  for (k in 2:(n - 1)) {
    if(t > 1){
      si = silhouette(getClusters(t, k, c, d), getDistance(d))
      mean_silhouette[k] = summary(si)$avg.width
    }else{
      mean_silhouette[k] = (pam(data,k,diss=F,stand=F))$silinfo$avg.width
    }
    print(mean_silhouette[k])
  }
  
  x11()
  par(mar=c(5.1,5.1,5.1,2.1))
  k.best = which.max(mean_silhouette)
  plot(1:(n-1), mean_silhouette, type="b", lwd=2, cex=1.2, font.lab=3, xlim=c(2,(n - 1)), cex.main=2, cex.lab=1.5, ylim=c(0,max(mean_silhouette)+0.1), col="grey", main="Silhouette plot for k groups", xlab="Nb. of clusters", ylab="Average silhouette width", axes=F)
  text(k.best, max(mean_silhouette), round(max(mean_silhouette),3), col="red", pos=4, cex=1.2)
  axis(1, seq(2,(max_cluster)), lwd=font_size, font.axis=font_size)
  axis(2, seq(0.0,(max(mean_silhouette)+0.1),0.1), lwd=font_size, font.axis=font_size)
  points(k.best, max(mean_silhouette), pch=19, col="red", cex=1.5)
  cat("","Silhouette-optimal number of clusters k =", k.best, "\n","with an average silhouette width of", round(max(mean_silhouette),4), "\n")
  abline(v=k.best, lty=2, col="red", lwd=2)
  return (k.best)
}
optimal_nb_clusters = plotAllSilhouette(typeClassif, max_cluster + 1, classif, data)

plotSmallClusterSilhouette = function(t, n, c=NULL, d=NULL){
  x11()
  par(mfrow=c(2,2),mar=c(5, 10, 4, 2))
  par(mar=c(5, 10, 2, 2))
  for (k in 2:4) {
    clusters = getClusters(t, k , c, d)
    if(t > 1){
      sil = silhouette(clusters, getDistance(d))
    }else{
      classif =  pam(data,k,diss=F,stand=F)
      sil = silhouette(clusters, classif$diss)
    }
    silo = sortSilhouette(sil)
    rownames(silo) = row.names(d)[attr(silo,"iOrd")]
    clusters=silo[,1]
    for (i in 1:4){
      clusters[clusters==i] = colPers(k)[i]
    }
    plot(silo, max.strlen=20, main=" ", cex.names=0.8, col=clusters, nmax.lab=100)
  }
  par(mfrow=c(1,1), mar=c(5, 4, 4, 2) + 0.1 )
}
plotSmallClusterSilhouette(typeClassif, max_cluster, classif, data)

#Plot a heatMap
#Input: s: an ordonned silhouette object
heatMap = function(s){
  matrix=as.matrix(distance_matrix)
  matrix=matrix[attr(s,"iOrd"),attr(s,"iOrd")]
  rownames(matrix) = rownames(data)[attr(s,"iOrd")]
  plotcolors(dmat.color(as.dist(matrix),colors=heat.colors(1000)),na.color="red",rlabels=rownames(data)[attr(silo,"iOrd")],clabels=rownames(data)[attr(silo,"iOrd")],border=0)
}

heatMap(silo)
################################
#          Dendrogram
################################

#cah
#Input: n, nb of clusters
plotDendrogram=function(n){
  par(margin);x11()
  #pdf("dendrogram.pdf")
  plot(classif, ylim=c(0,max(classif$height)),xlim=c(0,nrow(data)),hang=-1, cex.main=2, cex.lab=1.5,lwd=font_size,xlab="Metabolites", sub="",ylab="Distance Between-group",main="Dendrogram",font.lab=font_size,axes=F)
  axis(2, seq(0,max(classif$height)),lwd=font_size,font.axis=font_size,cex.axis=0.8)
  #abline(h=seq(0.0,max(classif$height),0.1), lty=3, col="grey")
  #abline(h=c(classif$height), lty=3, lwd=2, col="grey")
  #projection of the clusters
  rect.hclust(classif, k=n, border=colPers(n))
  #dev.off()
}

if(typeClassif>2) plotDendrogram(as.numeric(optimal_nb_clusters))

################################
#            PCA
################################

plotAcp=function(optimal_nb_clusters){
  acp=dudi.pca(data,scannf=F)
  x11()
  #pdf("pca.pdf")
  #par(mar=c(0,0,0,0))
  clusters=getClusters(optimal_nb_clusters)
  title=paste("Cumulated inertia: ",round((acp$eig[1]+acp$eig[2])/sum(acp$eig),4)*100,"%",sep="")
  #s.class(addaxes=F,acp$li,sub=title,possub="topright",csub=1.5,as.factor(getClusters(optimal_nb_clusters)),grid=F,col=colPers(optimal_nb_clusters))
  s.class(addaxes=F,acp$li,ylim=c(min(acp$li[,2])-3,max(acp$li[,2])+3),xlim=c(min(acp$li[,1])-3,max(acp$li[,1])+3),sub=title,possub="topright",csub=1.5,as.factor(clusters),grid=F,col=colPers(optimal_nb_clusters))
  abline(h=0,v=0,lty=2,lwd=2,col="grey")
  for (i in 1:optimal_nb_clusters){
    clusters[clusters==i] = colPers(optimal_nb_clusters)[i]
  }
  text(x=acp$li[,1],y=acp$li[,2],labels=rownames(acp$li),col=clusters,cex=1)
  #s.label(acp$li,add.plot = T,boxes=F,clabel=0.8,col="red")
  #dev.off()
}

plotAcp(2)


#heatmap.2(as.matrix(distance_matrix),distfun = function(x) dist(x,method = 'euclidean'),hclustfun = function(x) hclust(x,method = 'complete'),key=FALSE, density.info="none",lhei = c(0.05,0.95) ,cexRow=2,cexCol=2,margins=c(20,26),trace="none",srtCol=45,dendrogram="row")
#plot(data[,1],data[,2],type="p",pch="+",xlab="Axis 1", ylab="Axis 2",col=rainbow(optimal_nb_clusters)[clusters])
#points(classif$medoids[,1],classif$medoids[,2], cex=1.5,pch=16,col=c("red","blue","green")[1:3])
#s.label(data,boxes=F,clabel=0.8,neig = F,add.plot = T)
#points(data[,1],data[,2],col=c("lightcoral","skyblue","greenyellow")[classif$clustering])

#classif$cluster[order(-classif$cluster,decreasing=T)]