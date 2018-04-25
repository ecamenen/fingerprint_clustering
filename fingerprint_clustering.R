#!/usr/bin/env Rscript

#clean all objects
rm(list=ls())
setwd("~/bin/fingerprint_clustering")

library(cluster)
library(gclus)
library(ade4)

suppressPackageStartupMessages(require(optparse))

option_list = list(
  make_option(c("-i", "--infile"), type="character", default="matrix.txt", 
              help="Fingerprint file name [default: %default]"),
  make_option(c("-m", "--maxCluster"), type="integer", default=6, 
              help="Maximum number of clusters [default: Complete link"),
  make_option(c("-t", "--typeClassif"), type="integer", default=4, 
              help="Type of classifation [default: %default] (1: K-menoids; 2: K-means; 3: Ward; 4: Complete link; 5: UPGMA; 6: WPGMA"),
  make_option(c("-adv", "--advanced"), type="logical", action="store_true", 
              help="Activate advanced mode (print more outputs)"),
  make_option(c("-n", "--nbCluster"), type="integer", help="Fix the number of clusters")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

#Pseudo-random settings: 
#set.seed(1)
#milisec * PID
set.seed(as.numeric(format(Sys.time(), "%OS2"))*100 * Sys.getpid())

#global variables
#choix du niveau de coupure
nb_clusters=opt$nbCluster
font_size=3
max_cluster=opt$maxCluster
typeClassif=opt$typeClassif
advanced = "advanced" %in% names(opt)


################################################
#     Data test: random distance matrix   
################################################

#data=setRandomDataSet()
data = read.table(opt$infile, header=F, sep="\t", dec=".", row.names=1)
colnames(data) = rownames(data)

#conversion into distance
#distance_matrix=as.dist(data)
getDistance = function(d, t, k=NULL){
  if (t > 1) dist(d, method = "euclidian")
  else getCNH(t,d,k)$diss
}
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


#Inputs: x : a matrix
#filename of the saved file
#Prints the matrix, save the matrix
writeTsv = function(x,f, h=TRUE){
  options(warn = -1)
  if(h==TRUE) output=as.matrix(rbind(c("", colnames(x)), cbind(rownames(x),x)))
  else output = x
  output[is.na(output)] = ""
  colnames(output)=rep("", ncol(output)); rownames(output)=rep("", nrow(output))
  print(output,row.names=FALSE, col.names=FALSE, quote=F)
  write(t(output), file=f, ncolumns=ncol(output), sep="\t")
  #write.table(x, f, na = "",col.names = colnames(x),row.names = rownames(x),append = F,sep = "\t")
  options(warn = 0)
}
############################################################
#          Cophenetic (dendrogram distance matrix)
############################################################

#Inputs:
# d : data
# cah : hierarchical classification
plotCohenetic=function(t, d,cah){
  library(scales)
  dis = getDistance(d, t)
  coph_matrix = cophenetic(cah)
  cor_coph = cor(dis, coph_matrix)
  cat(paste("\nCOPHENETIC:\n% of explained variance:", round(cor_coph^2,3)))
  #x11()
  pdf("shepard_graph.pdf")
  par(mar=c(5.1,5.1,4.1,2.1))
  plot(dis, coph_matrix, pch=19,col=alpha("red",0.2), cex=1, axes=F, cex.lab=1.5, cex.main=2, font.lab=font_size, xlim=c(0,max(dis)), ylim=c(0,max(coph_matrix)), xlab="Distance between metabolites",ylab="Cophenetic distance", asp=1, main=paste("Cophenetic correlation: ",round(cor_coph,3)))
  axis(2, seq(0.0,max(coph_matrix),1), lwd=font_size, font.axis=3, cex.axis=0.8)
  axis(1, seq(0,max(dis),1), lwd=font_size, font.axis=3, cex.axis=0.8)
  abline(0,1, col="grey", lwd=font_size, lty=2)
  suprLog = dev.off()
}

if(typeClassif>2) plotCohenetic(typeClassif, data, classif)

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
  if(t==2) inertia = getCumulatedBetweenInertiaPerCluster(t, n, d=d)
  if(t>2) inertia = getCumulatedBetweenInertiaPerCluster(t, n, c)
  k.best = which.max(inertia)
  pdf("cumulated_between.pdf")
  #x11()
  par(mar=c(5.1,5.1,4.1,2.1))
  plot(inertia, type="b", lwd=font_size, font.lab=3, cex.lab=font_size/2, cex.main=font_size/1.5, col="grey", xlim=c(2, n), ylim=c(0,(max(inertia)+5)), axes=F, xlab="Nb. of clusters", ylab="Cumulated between-group inertia")
  axis(1, seq(2,n), lwd=font_size, font.axis=3, cex.axis=0.8)
  axis(2, seq(0,max(inertia)+5,10),lwd=font_size, font.axis=3, cex.axis=0.8)
  text(k.best, max(inertia), paste("",round(max(inertia),4),sep="\n \n"), col="red", pos=2, cex=font_size/2.5)
  points(k.best, max(inertia), pch=20, col="red", cex=font_size)
  abline(v=k.best, lty=2 ,col="red", lwd=font_size/1.5)
  cat("\n\nBETWEEN-INERTIA:\nOptimal number of clusters k = ", k.best, "\n","With a value of ", round(max(inertia),4), "\n", sep="")
  suprLog = dev.off()
}
#plotCumulatedBetweenInertia(typeClassif, max_cluster, classif, data)

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

################################
#          Fusion levels
################################

#Plot fusion graph
plot_fusion_levels = function(t, n, c=NULL, d=NULL) {
  subset_height = getBetweenInertia(t, n, c, d)
  height_diff = getBetweenDifference(t, n, c, d)
  #x11()
  pdf("fusion_levels.pdf")
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
  cat("\n\nINERTIA BREAKS:\nOptimal number of clusters k = ", optimal_nb_clusters, "\n","With a value of ", max(rev(round(height_diff,2))), "\n", sep="")
  #catch_printing=identify(x=classif$height[-1], y=(nrow(data)-1):2,labels=paste(round(height_diff[-1],digits=2), result[-(nrow(data)-1),2], sep="\n"),col="red", cex=0.8,plot=T)
  suprLog = dev.off()
}

if(typeClassif > 1) plot_fusion_levels(typeClassif, max_cluster, classif, data)

################################
#          Silhouette
################################

#Ouput: an organised silhouette object
getSilhouette = function(t, k , c, d){
  clusters = getClusters(t, k , c, d)
  diss = getDistance(d,t,k)
  sil = sortSilhouette(silhouette(clusters, diss))
  rownames(sil) = row.names(d)[attr(sil,"iOrd")]
  return (sil)
}

cat("\nSILHOUETTE:\n")
# Plot the best average silhouette width for all clustering possible
plotAverageSilhouette = function(t, n, c=NULL, d=NULL){
  
  mean_silhouette = numeric(n - 1)
  for (k in 2:(n - 1)) {
    si = getSilhouette(t, k , c, d)
    mean_silhouette[k] = summary(si)$avg.width
    cat(paste("G",k, ": ", round(mean_silhouette[k],3), "\n",sep=""))
  }
  
  #x11()
  pdf("average_silhouettes.pdf")
  par(mar=c(5.1,5.1,5.1,2.1))
  k.best = which.max(mean_silhouette)
  plot(1:(n-1), mean_silhouette, type="b", lwd=2, cex=1.2, font.lab=3, xlim=c(2,(n - 1)), cex.main=2, cex.lab=1.5, ylim=c(0,max(mean_silhouette)+0.1), col="grey", main="Silhouette plot for k groups", xlab="Nb. of clusters", ylab="Average silhouette width", axes=F)
  text(k.best, max(mean_silhouette), round(max(mean_silhouette),3), col="red", pos=4, cex=1.2)
  axis(1, seq(2,(max_cluster)), lwd=font_size, font.axis=font_size)
  axis(2, seq(0.0,(max(mean_silhouette)+0.1),0.1), lwd=font_size, font.axis=font_size)
  points(k.best, max(mean_silhouette), pch=19, col="red", cex=1.5)
  cat("Optimal number of clusters k = ", k.best, "\n","With an average silhouette width of ", round(max(mean_silhouette),4), "\n", sep="")
  abline(v=k.best, lty=2, col="red", lwd=2)
  suprLog = dev.off()
  return (k.best)
}
optimal_nb_clusters = plotAverageSilhouette(typeClassif, max_cluster + 1, classif, data)
sil = getSilhouette(typeClassif, optimal_nb_clusters, classif, data)

colorClusters = function(cl){
  nb_clusters = length(levels(as.factor(cl)))
  for (i in 1:nb_clusters){
    cl[cl==i] = colPers(nb_clusters)[i]
  }
  return (cl)
}

plotSilhouette = function(s){
  #x11()
  pdf("silhouette.pdf")
  par(mar=c(4, 8, 3, 2))
  plot(s, max.strlen=20, main=" ", sub= "", do.clus.stat=FALSE, cex.lab=font_size/2, font.lab=3, xlab="Silhouette width", cex.names=0.8, col=colorClusters(s[,1]), nmax.lab=100, do.n.k = FALSE, axes=F)
  mtext(paste("Average silhouette width:", round(summary(s)$avg.width,3)), font=2, cex=font_size/2, line=1)
  axis(1, seq(0,1,by=0.2), lwd=font_size, font.axis=3, cex.axis=0.8)
  suprLog = dev.off()
}
plotSilhouette(sil)

dis = getDistance(data, typeClassif, optimal_nb_clusters)

clusters = getClusters(typeClassif, optimal_nb_clusters, classif, data)

#Plot a heatMap
#Input:
# d: a distance object
# s: an organised silhouette object
heatMap = function(d, s){
  #x11()
  pdf("heat_map.pdf")
  par(mar=c(1, 8, 8, 1))
  matrix=as.matrix(d)
  matrix=matrix[attr(s,"iOrd"),attr(s,"iOrd")]
  rownames(matrix) = rownames(data)[attr(s,"iOrd")]
  labels = attr(d, "Labels")[attr(s,"iOrd")]
  plotcolors(dmat.color(as.dist(matrix), colors=heat.colors(1000)), na.color="red", rlabels=labels, clabels=labels, border=0)
  suprLog = dev.off()
}
heatMap(dis, sil)

################################
#          Dendrogram
################################

# Inputs:
# k: number of clusters
# c: hierarchical classification
plotDendrogram = function(c, k){
  #x11()
  pdf("dendrogram.pdf")
  par(mar=c(2,5,5,1))
  plot(c, ylim=c(0,max(c$height)), xlim=c(0,length(c$labels)), hang=-1, cex.main=2, cex.lab=1.5, lwd=font_size, sub="", ylab="Distance Between-group", main="Dendrogram", font.lab=font_size, axes=F)
  axis(2, seq(0,max(c$height)), lwd=font_size, font.axis=font_size, cex.axis=0.8)
  #projection of the clusters
  rect.hclust(c, k=as.numeric(k), border=colPers(k))
  suprLog = dev.off()
}
if(typeClassif > 2) plotDendrogram(classif, optimal_nb_clusters)

################################
#            PCA
################################

plotPca = function(t, k, c, d){
  pca = dudi.pca(d, scannf=F)
  #x11()
  pdf("pca.pdf")
  #par(mar=c(0,0,0,0))
  clusters = getClusters(t, k, c, d)
  title = paste("Cumulated inertia:", round((pca$eig[1]+pca$eig[2])/sum(pca$eig),4)*100, "%")
  s.class(addaxes=F, pca$li ,ylim=c(min(pca$li[,2])-3, max(pca$li[,2])+3), xlim=c(min(pca$li[,1])-3, max(pca$li[,1])+3), csub=1.5, as.factor(clusters), grid=F, col=colPers(optimal_nb_clusters))
  mtext(title, font=2, cex=font_size/2, line=1)
  abline(h=0, v=0, lty=2, lwd=2, col="grey")
  text(x=pca$li[,1], y=pca$li[,2], labels=rownames(pca$li), col=colorClusters(clusters), cex=1)
  suprLog = dev.off()
}

plotPca(typeClassif, optimal_nb_clusters, classif, data)


################################
#            CTR
################################

scalecenter <- function(d) {
  N = nrow(d) ; d = scale(d);
  return(d * sqrt(N/(N-1)))
}

# Inputs:
# d: data
# cl: clusters object
getBetweenPerVariable = function(d, cl){
  #get percent values in output
  d = scalecenter(d)
  nb_cl = length(levels(as.factor(cl)))
  nb_met = length(cl)
  ctr = matrix(0, nrow=nb_cl, ncol=nb_met)
  for (i in 1:nb_met) {
    #get the group number for each row
    cli = cl[i]
    #in the dataset, for a metabolite row, loop an each metadabolite column
    #values are affected the corresponding cluster row and metabolite column in ctr
    for (j in 1:nb_met) ctr[cli,j] = ctr[cli,j] + d[i,j]
  }
  return (ctr)
}

# Relative contributions of the metabolites to inertia of each clusters (CTR)
# Inputs:
# t: number of type of classification
# k: number of clusters
# c: hierarchical classification
# d: data
getCtrVar = function(t, k, c, d) {
  cl = getClusters(t, k, c, d)
  nb_cl = length(levels(as.factor(cl)))
  nb_met = length(cl)
  
  ctr = getBetweenPerVariable(d, cl)
  rownames(ctr) = paste("G", seq(1, k), sep=""); colnames(ctr) = colnames(d)
  for (i in 1:nb_cl)
    for (j in 1:nb_met) ctr[i,j] = ctr[i,j]^2 / (nb_met * length(cl[cl==i]))
  
  return(ctr)
}

if (advanced == TRUE){
  ctrVar = round(1000 * getCtrVar(typeClassif, optimal_nb_clusters, classif, data) / 10)
  cat("\nCONTRIBUTION:")
  writeTsv(ctrVar,"relative_ctr.tsv")
}

################################
#            PDIS
################################

# Discriminant power (PDIS)
# Relative contributions of the metabolites to inertia of a partitionning
# Inputs: 
# t: number of type of classification
# n: number max of clusters
# c: hierarchical classification
# d: data
getPdis = function(t, k, c, d) {
  
  nb_met = nrow(d)
  ctrVar = getCtrVar(t, k, c, d)
  
  #for each metabolite contribution (in column), sum the k clusters values
  pdis = vector(mode="numeric", nb_met)
  for (i in 1:nb_met) pdis[i] = sum(ctrVar[,i])
  return(pdis)
}

# Inputs: 
# t: number of type of classification
# n: number max of clusters
# c: hierarchical classification
# d: data
# index: pdis or rho2 calculation
getIndexPerPartition = function(t, n, c, d, index){
  
  if (index == "pdis") nb_col = ncol(d)
  else nb_col = n
  index_per_partition = matrix(NA, n-1, nb_col)
  rownames(index_per_partition) = seq(2, n)
  
  for (k in 2:n){
    if (index == "pdis"){
      colnames(index_per_partition) = colnames(d)
      res = getPdis(t, k, c, d)
    }else{
      colnames(index_per_partition) = paste("G", seq(1, n), sep="")
      res = getRho2(t, k, c, d)
    }
    for(i in 1:length(res)){
      index_per_partition[k-1, i] = round(res[i], 2)
    }
  }
  return (index_per_partition)
}

if (advanced == TRUE){
  pdis_per_partition = getIndexPerPartition(typeClassif, max_cluster, classif, data, "pdis")
  cat("\nDISCRIMINANT POWER:")
  writeTsv(pdis_per_partition, "discriminant_power.tsv")
}

#########################################
#            Excentricity (RHO2)
#########################################

# Distance**2 of each cluster from the data center (RHO2)
# Inputs: 
# t: number of type of classification
# n: number max of clusters
# c: hierarchical classification
# d: data
getRho2 = function(t, k, c, d) {
  
  cl = getClusters(t, k, c, d)
  nb_cl = length(levels(as.factor(cl)))
  nb_met = length(cl)
  ctr = getBetweenPerVariable(d, cl)
  
  for (i in 1:nb_cl) {
    cli = cl[i]
    for (j in 1:nb_met) ctr[cli,j] = ctr[cli,j] + d[i,j]
  }
  
  for (i in 1:k)
    for (j in 1:nb_met) ctr[i,j] = ctr[i,j]/length(cl[cl==i])
  
  rho2 = vector(mode="numeric", k)
  for (i in 1:k) rho2[i] = sum(ctr[i,]^2)
  
  return(rho2)
}

if (advanced == TRUE){
  excentricity = getIndexPerPartition(typeClassif, max_cluster, classif, data, "rho2")
  cat("\nEXCENTRICIY:")
  writeTsv(excentricity,"excentricity.tsv")
}

################################
#            CLUSTERS
################################

writeClusters = function(c, f){
  nb_cl = length(levels(as.factor(c)))
  output = matrix(NA, length(c), nb_cl)
  for (i in 1:nb_cl ){
    output[c(1:length(c[c==i])),i] <- names(c[c==i])
  }
  writeTsv(output, f, h=FALSE)
}
cat("\nCLUSTERS:")
writeClusters(clusters, "clusters.tsv")

#Input:
# cl: clusters
getClusterSizes = function(cl){
  cluster_sizes = table(cl)
  cluster_sizes = data.frame(cluster_sizes)[,2]
  cluster_sizes = data.frame(cluster_sizes)
  names(cluster_sizes) = "Effectif"
  return(cluster_sizes)
}
#cluster_sizes=getClusterSizes(clusters)
#writeTsv(cluster_sizes,"cluster_sizes.tsv")