getArgs = function(){
  option_list = list(
    make_option(c("-w", "--workdir"), type="character", metavar="character",
                help="Working directory path [default: the folder where the script is launched]"),
    make_option(c("-i", "--infile"), type="character", default="matrix.txt", 
                metavar="character",
                help="Fingerprint file name [default: %default]"),
    make_option(c("-m", "--maxCluster"), type="integer", default=6, metavar="integer",
                help="Maximum number of clusters [default: %default]"),
    make_option(c("-t", "--typeClassif"), type="integer", default=2, metavar="integer",
                help="Type of classifation [default: %default] (1: K-menoids; 2: K-means; 3: Ward; 4: Complete link; 5: UPGMA; 6: WPGMA"),
    make_option(c("-adv", "--advanced"), type="logical", action="store_true", 
                help="Activate advanced mode (print more outputs)"),
    make_option(c("-q", "--quiet"), type="logical", action="store_true",
                help="Activate quiet mode"), 
    make_option(c("-n", "--nbCluster"), type="integer", metavar="integer",
                help="Fix the number of clusters"),
    make_option(c("-r", "--ranked"), type="logical", action="store_true", 
                help="Rank the metabolites in clusters by silhouette scores instead of alphabetically")
    )
  
  return (OptionParser(option_list=option_list))
}

#Check the arguments validity
#Inputs:
# a: arguments (optionParser object)
checkArg = function(a){
  opt = parse_args(a)
  # o: one argument from the list of arguments
  # def: defaul message
  
  checkMinCluster = function (o, def="")
  if (opt[[o]] < 2){
    print_help(a)
    stop(paste("--",o ," must be upper or equal to 2",def,".\n",sep=""), call.=FALSE)
  }
  checkMinCluster("maxCluster"," [by default: 6]")
  if(!is.null(opt$nbCluster)) checkMinCluster("nbCluster")
  
  if ((opt$typeClassif < 1) || (opt$typeClassif > 6)){
    print_help(a)
    stop("--typeClassif must be comprise between 1 and 6 [by default: 2].\n", call.=FALSE)
  }
  
  checkFile = function (o){
    if(!file.exists(opt[[o]])){
      print_help(a)
      stop(paste("--", o, " name does not exist\n", sep=""), call.=FALSE)
    }
  }
  if(!is.null(opt$workdir)) checkFile("workdir")
  if(!is.null(opt$infile)) checkFile("infile")
}

#Checking clusters args after data loading
#Inputs:
# a: arguments (optionParser object)
# d: data
# o: one argument from the list of arguments
# def: defaul message
postChecking = function (a, d){
  
  opt = parse_args(a)
  
  checkMaxCluster = function (o, def="")
    if (opt[[o]] > nrow(d)){
      print_help(a)
      stop(paste("--", o," must be lower or equal to the fingerprint",def,".\n",sep=""), call.=FALSE)
    }
  
  checkMaxCluster("maxCluster"," [by default: 6]")
  if(!is.null(opt$nbCluster)) checkMaxCluster("nbCluster")
}

#Usage: colPers(x), x a number of colours in output
#Gradient of color
colPers = colorRampPalette(c(rgb(0.6,0.1,0.5,1), rgb(1,0,0,1), rgb(0.9,0.6,0,1), rgb(0.1,0.6,0.3,1), rgb(0.1,0.6,0.5,1), rgb(0,0,1,1)), alpha = TRUE)

scalecenter = function(d) {
  N = nrow(d) ; d = scale(d);
  return(d * sqrt(N/(N-1)))
}

getDistance = function(d, t, k=NULL){
  if (t > 1) dist(d, method = "euclidian")
  else getCNH(t,d,k)$diss
}

#Inputs: x : a matrix
#filename of the saved file
#Prints the matrix, save the matrix
writeTsv = function(x, h=TRUE){
  #print on stdout
  if (v==T) cat(paste("\n", gsub("_", " ", toupper(x)), ":", sep=""))
  #disabling warning
  options(warn = -1)
  #get variable
  tab = get(x)
  if(h==TRUE) output=as.matrix(rbind(c("", colnames(tab)), cbind(rownames(tab),tab)))
  else output = tab
  #discard empty rows
  output = output[rowSums(is.na(output)) != ncol(output),]
  #TODOD:
  #output = output[,colSums(is.na(output)) != nrow(output)]
  output[is.na(output)] = ""
  colnames(output)=rep("", ncol(output)); rownames(output)=rep("", nrow(output))
  if (v==T)  print(output, row.names=FALSE, col.names=FALSE, quote=F)
  write(t(output), paste(x,".tsv",sep=""), ncolumns=ncol(output), sep="\t")
  #write.table(x, f, na = "",col.names = colnames(x),row.names = rownames(x),append = F,sep = "\t")
  options(warn = 0)
}

################################
#          Graphic
################################


printAxis = function (side, min, max){
  axis(side, seq(min,max), lwd=3)
}


#f: filename
savePdf = function (f){
  pdf(f)
  setGraphic()
}

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

#Inputs: 
# t: number of type of classification
# d: data (or distance matrix for hierarchic)
# k: number of clusterting
#Ouput: Non-hierarchical classification
getCNH = function(t, d, k){
  if (t==1) return (pam(d, k, diss=F, stand=F))
  else if (t==2) return (kmeans(d, centers=k, nstart=100))
}

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

#Input:
# cl: clusters
colorClusters = function(cl){
  nb_clusters = length(levels(as.factor(cl)))
  for (i in 1:nb_clusters){
    cl[cl==i] = colPers(nb_clusters)[i]
  }
  return (cl)
}

#Inputs:
# cl: clusters
# f : filename
# r: ordered alphabetically
writeClusters = function(cl, r=FALSE){
  nb_cl = length(levels(as.factor(cl)))
  clusters = matrix(NA, length(cl), nb_cl)
  for (i in 1:nb_cl ){
    if (r == FALSE){
      #get metabolites from clusters and put into a column of the output matrix
      # from the begining of the column to the end of the vector of metabolites names
      clusters[c(1:length(cl[cl==i])),i] = names(cl[cl==i])
    }else if (r == TRUE){
      #ordering alphabetically
      clusters[c(1:length(cl[cl==i])),i] = sort(names(cl[cl==i]))
    }
    #ordering by clusters size
    length_cl = colSums(!is.na(clusters))
    for (i in 2:nb_cl) {
      #inversion if a column as more metabolites than the previous
      if (length_cl[i] > length_cl[i-1]){
        temp = clusters[,i-1]
        clusters[,i-1] = clusters[,i]
        clusters[,i] = temp
      }
    }
  }
  #dirty way to force saving a local variable
  # (because writeTsv use only global variables)
  assign("clusters", clusters,.GlobalEnv)
  writeTsv("clusters", h=FALSE)
}

############################################################
#          Cophenetic (dendrogram distance matrix)
############################################################

#Inputs:
# d : data
# cah : hierarchical classification
plotCohenetic=function(t, d,cah){
  dis = getDistance(d, t)
  coph_matrix = cophenetic(cah)
  cor_coph = cor(dis, coph_matrix)
  if (v==T) cat(paste("\nCOPHENETIC:\nExplained variance (%):", round(cor_coph^2,3), "\nCorrelation:",round(cor_coph,3),"\n"))
  #x11()
  savePdf("shepard_graph.pdf")
  plot(dis, coph_matrix, pch=19, col=alpha("red",0.2), axes=F, xlim=c(0,max(dis)), ylim=c(0,max(coph_matrix)), xlab="Distance between metabolites",ylab="Cophenetic distance", asp=1, main=paste("Cophenetic correlation: ",round(cor_coph,3)))
  printAxis(2, 0, max(coph_matrix))
  printAxis(1, 0, max(dis))
  abline(0, 1, col="grey", lwd=3, lty=2)
  suprLog = dev.off()
}

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
    c = getCNH(t, d, k)
    return (round(c$betweenss/c$totss,5) * 100)
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
    return(round(100 * sum_inertia / sum(c$height),3))
  }
}

getTotInertia = function(t, c, d) {
    if(t > 2 ) sum(getBetweenInertia(t, nrow(d), c, d))
    else if (t == 2) getCNH(typeClassif, data, 2)$totss
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
      inertia[k-1] = getCNH(t, d, k)$betweenss
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
  if (t > 2) return (inertia)
  else if (t==2) return (inertia[-1])
}

# Inputs: 
# t: number of type of classification
# n: maximum number of clusters
# c: hierarchical classification
# d: data
getBetweenDifference = function(t, n, c=NULL, d=NULL){
  if(t==2) inertia = getBetweenInertia(t, n+1, d=d)[-1]
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
  if (t==2) rownames(ranked_inertia_diff) = as.numeric(rownames(ranked_inertia_diff)) + 1
  else if (t >2) rownames(ranked_inertia_diff) = n - as.numeric(rownames(ranked_inertia_diff)) + 1
  return(ranked_inertia_diff)
}

printTableInertia = function(t, n, c=NULL, d=NULL){
  options(warn = -1)
  table_inertia = cbind(getBetweenInertia(t, n, c, d)[-1], getBetweenDifference(t, n, c, d)) / getTotInertia(t, c, d)
  #outputs are reversed comparatively to CumulatedBetween outputs
  if (t > 2) for(i in 1:ncol(table_inertia)) table_inertia[,i] = rev(table_inertia[,i])
  table_inertia = cbind(table_inertia*100, getCumulatedBetweenInertiaPerCluster(t, n, c, d))
  rownames(table_inertia) = seq(2, n)
  colnames(table_inertia) = c("Between-inertia (%)", "Differences (%)","Cumulated inertia (%)")
  table_inertia = round(table_inertia, 3)
  options(warn = 0)
  return (table_inertia)
}

# Between inertia differences between a partionning and the previous
plot_fusion_levels = function(t, n, c=NULL, d=NULL) {
  subset_height = rev((getBetweenInertia(t, n, c, d)[-1] / getTotInertia(t, c, d)) *100)
  height_diff = (getBetweenDifference(t, n, c, d) / getTotInertia(t, c, d))*100
  if (t==2) height_diff = rev(height_diff)
  #x11()
  savePdf("fusion_levels.pdf")
  plot(2:n, subset_height, type="b", ylim=c(round(min(subset_height))-1,round(max(subset_height))+1), xlim=c(2,n), xlab="Nb. of clusters", ylab="Between-group inertia (%)", col="grey", axes=F)
  title(main="Fusion levels", line=2, cex.main=3/1.5)
  mtext("(in red, between-group differences with the previous clustering)", side=3, line=1)
  printBestClustering(subset_height, 
                      "difference with the next partitionning",
                      rev(round(height_diff,3)))
  #catch_printing=identify(x=classif$height[-1], y=(nrow(data)-1):2,labels=paste(round(height_diff[-1],digits=2), result[-(nrow(data)-1),2], sep="\n"),col="red", cex=0.8,plot=T)
  suprLog = dev.off()
}


setGraphic = function(){
  par(font.axis=3, cex.axis=0.8, cex.lab=3/2, cex.main=1.5, cex=1, font.lab=3, lwd=3, mar=c(5.1,5.1,5.1,2.1))
}

printBestClustering = function(pointValue, valueType, textValue){
  printAxis(1, 2, max_cluster)
  printAxis(2, round(min(pointValue)), round(max(pointValue)))
  abline(v=optimal_nb_clusters, col="red", lty=2, lwd=3/1.5)
  points(optimal_nb_clusters, pointValue[1], pch=19, col="red", cex=3/1.5)
  text(y=pointValue, x=2:max_cluster, labels=textValue, cex=1.2, pos=4, col="red")
  if (v==T) cat("Optimal number of clusters k = ", optimal_nb_clusters, "\n","With a ", valueType, " of ", max(textValue), "\n", sep="")
}

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

# Plot the best average silhouette width for all clustering possible
plotAverageSilhouette = function(t, n, c=NULL, d=NULL){
  if (v==T) cat("\nSILHOUETTE:\n")
  mean_silhouette = numeric(n - 1)
  for (k in 2:(n - 1)) {
    si = getSilhouette(t, k , c, d)
    mean_silhouette[k] = summary(si)$avg.width
    if (v==T) cat(paste("", k, ": ", round(mean_silhouette[k],3), "\n",sep=""))
  }
  
  #x11()
  pdf("average_silhouettes.pdf")
  par(mar=c(5.1,5.1,5.1,2.1))
  k.best = which.max(mean_silhouette)
  plot(1:(n-1), mean_silhouette, type="b", lwd=2, cex=1.2, font.lab=3, xlim=c(2,(n - 1)), cex.main=2, cex.lab=1.5, ylim=c(0,max(mean_silhouette)+0.1), col="grey", main="Silhouette plot for k groups", xlab="Nb. of clusters", ylab="Average silhouette width", axes=F)
  text(k.best, max(mean_silhouette), round(max(mean_silhouette),3), col="red", pos=4, cex=1.2)
  axis(1, seq(2,(max_cluster)), lwd=3, font.axis=3)
  axis(2, seq(0.0,(max(mean_silhouette)+0.1),0.1), lwd=3, font.axis=3)
  points(k.best, max(mean_silhouette), pch=19, col="red", cex=1.5)
  if (v==T) cat("Optimal number of clusters k = ", k.best, "\n","With an average silhouette width of ", round(max(mean_silhouette),3), "\n", sep="")
  abline(v=k.best, lty=2, col="red", lwd=2)
  suprLog = dev.off()
  return (k.best)
}

plotSilhouette = function(s){
  #x11()
  pdf("silhouette.pdf")
  par(mar=c(4, 8, 3, 2))
  plot(s, max.strlen=20, main=" ", sub= "", do.clus.stat=FALSE, cex.lab=3/2, font.lab=3, xlab="Silhouette width", cex.names=0.8, col=colorClusters(s[,1]), nmax.lab=100, do.n.k = FALSE, axes=F)
  mtext(paste("Average silhouette width:", round(summary(s)$avg.width,3)), font=2, cex=3/2, line=1)
  axis(1, seq(0,1,by=0.2), lwd=3, font.axis=3, cex.axis=0.8)
  suprLog = dev.off()
}

################################
#          HEATMAP
################################
#Inputs:
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
  plot(c, ylim=c(0,max(c$height)), xlim=c(0,length(c$labels)), hang=-1, cex.main=2, cex.lab=1.5, lwd=3, sub="", ylab="Distance Between-group", main="Dendrogram", font.lab=3, axes=F)
  axis(2, seq(0,max(c$height)), lwd=3, font.axis=3, cex.axis=0.8)
  #projection of the clusters
  rect.hclust(c, k=as.numeric(k), border=colPers(k))
  suprLog = dev.off()
}

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
  mtext(title, font=2, cex=3/2, line=1)
  abline(h=0, v=0, lty=2, lwd=2, col="grey")
  text(x=pca$li[,1], y=pca$li[,2], labels=rownames(pca$li), col=colorClusters(clusters), cex=1)
  suprLog = dev.off()
}

################################
#            CTR
################################

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

################################
#            MAIN
################################

#Pseudo-random settings: 
#milisec * PID
set.seed(as.numeric(format(Sys.time(), "%OS2"))*100 * Sys.getpid())

#Loading librairies
librairies = c("cluster", "optparse", "gclus", "ade4", "scales")
for (l in librairies){
  if (! (l %in% installed.packages()[,"Package"])) install.packages(l, repos = "http://cran.us.r-project.org")
  library(l, character.only = TRUE)
}

#Global variables
args = getArgs()
checkArg(args)
opt = parse_args(args)
nb_clusters = opt$nbCluster
max_cluster = opt$maxCluster
typeClassif = opt$typeClassif
advanced = "advanced" %in% names(opt)
v = !("quiet" %in% names(opt))
ranked = !("ranked" %in% names(opt))
if (!is.null(opt$workdir)) setwd(opt$workdir)
#setwd("~/bin/fingerprint_clustering/")

#Loading data
data = read.table(opt$infile, header=F, sep="\t", dec=".", row.names=1)
colnames(data) <- substr(rownames(data), 1, 35) -> rownames(data)
postChecking(args, data)

#Perform classification
classif = getCAH(data, typeClassif)
if(typeClassif>2) plotCohenetic(typeClassif, data, classif)

#Between inertia analysis
between_inertia = printTableInertia(typeClassif, max_cluster, classif, data)
writeTsv("between_inertia")
optimal_nb_clusters = as.numeric(rownames(getRankedInertia(typeClassif, max_cluster, classif, data))[1])
if(typeClassif > 1) plot_fusion_levels(typeClassif, max_cluster, classif, data)

#Silhouette analysis
optimal_nb_clusters = plotAverageSilhouette(typeClassif, max_cluster + 1, classif, data)
if(!is.null(nb_clusters)) optimal_nb_clusters = nb_clusters
sil = getSilhouette(typeClassif, optimal_nb_clusters, classif, data)
plotSilhouette(sil)

#Global variables
dis = getDistance(data, typeClassif, optimal_nb_clusters)
clusters = getClusters(typeClassif, optimal_nb_clusters, classif, data)

#Plots
heatMap(dis, sil)
if(typeClassif > 2) plotDendrogram(classif, optimal_nb_clusters)
plotPca(typeClassif, optimal_nb_clusters, classif, data)

#Advanced indexes
if (advanced == TRUE){
  
  contribution = round(1000 * getCtrVar(typeClassif, optimal_nb_clusters, classif, data) / 10)
  discriminant_power = getIndexPerPartition(typeClassif, max_cluster, classif, data, "pdis")
  excentricity = getIndexPerPartition(typeClassif, max_cluster, classif, data, "rho2")
  
  for (i in c("contribution", "discriminant_power", "excentricity"))
    writeTsv(i)
}

#Final outputs
writeClusters(clusters, ranked)
if (v != T) cat(paste("Clustering done.\nOptimal number of clusters choosen:", optimal_nb_clusters,"\n"))
