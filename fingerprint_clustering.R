getArgs = function(){
  option_list = list(
    make_option(c("-w", "--workdir"), type="character", metavar="character",
                help="Working directory path [default: the folder where the script is launched]"),
    make_option(c("-i", "--infile"), type="character", default="matrix.txt", 
                metavar="character",
                help="Fingerprint file name [default: %default]"),
    make_option(c("-m", "--maxCluster"), type="integer", default=6, metavar="integer",
                help="Maximum number of clusters [default: %default]"),
    make_option(c("-t", "--typeClassif"), type="integer", default=4, metavar="integer",
                help="Type of classifation [default: Complete links] (1: K-menoids; 2: K-means; 3: Ward; 4: Complete links; 5: UPGMA; 6: WPGMA"),
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

#Get the normalized distance between each points and the center
#Outputs:
# for each column, the mean=0 and the variance is the same
scalecenter = function(d) {
  #output scale function: for each column, mean=0, sd=1
  return(scale(d) * sqrt(nrow(d)/(nrow(d)-1)))
  # without multiplying by this constante, for advanced outputs, total (max_cluster=nrow(data)) will be different from 1
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

setGraphic = function(){
  setGraphicBasic()
  par(mar=c(5.1,5.1,5.1,2.1))
}

setGraphicBasic = function(){
  par(cex.lab=1.5, font.lab=3, font.axis=3, cex.axis=0.8, cex.main=2, cex=1, lwd=3)
}

printAxis = function (side, min, max, interval = 1){
  axis(side, seq(min,max, interval), lwd=3)
}

printBestClustering = function(sub_title, values, values_type, optimal_nb_clusters, interval = 1){
  printAxis(1, 2, max_cluster)
  if (interval >= 1){ axisSeq=round(values)
  }else{ axisSeq = c(0, max(values) +0.1)}
  printAxis(2, min(axisSeq), max(axisSeq), interval)
  title(main="Optimal number of clusters", line=1, cex.main=2)
  mtext(text=sub_title, font=3, cex=1.2, line = -1)
  abline(v=optimal_nb_clusters, col="red", lty=2, lwd=2)
  points(optimal_nb_clusters, max(values), pch=19, col="red", cex=2)
  text(y=values, x=2:max_cluster, labels=round(values,2), cex=1.2, pos=4, col="red")
  if (v==T) cat("Optimal number of clusters k = ", optimal_nb_clusters, "\n","With a", values_type, " of ", round(max(values),2), "\n", sep="")
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
    #cah: classification hierarchic ascending
    cah = hclust(dis, method=getTypeClassif(t))
  #automaticly ordering by clusters
  return (reorder.hclust(cah, d))
  }
  #TODO: exit if 0 < t < 6
}

#Inputs:
# t: number of type of classification
getTypeClassif = function(t){
  if (t==3) "ward.D2"
  else if (t==4) "complete"
  else if (t==5) "average"
  else if (t==6) "mcquitty"
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

# Distance matrix between each leaf of the dendogramm
#Inputs:
# d : data
# cah : hierarchical classification
plotCohenetic=function(t, d,cah){
  dis = getDistance(d, t)
  coph_matrix = cophenetic(cah)
  cor_coph = cor(dis, coph_matrix)
  if (v==T) cat(paste("\nCOPHENETIC:\nExplained variance (%):", round(cor_coph^2,3), "\nCorrelation with the data:",round(cor_coph,3),"\n"))
  #x11()
  savePdf("shepard_graph.pdf")
  plot(dis, coph_matrix, pch=19, col=alpha("red",0.2), axes=F, xlim=c(0,max(dis)), ylim=c(0,max(coph_matrix)), xlab="Distance between metabolites",ylab="Cophenetic distance", asp=1, main=paste("Cophenetic correlation: ",round(cor_coph,3)))
  printAxis(2, 0, max(coph_matrix))
  printAxis(1, 0, max(dis))
  abline(0, 1, col="grey", lwd=3, lty=2)
  suprLog = dev.off()
}

##############################################
#          Inertia
##############################################
# Raw inter-group inertia from a CAH
# Inputs:
# k: number of clusters
# c: hierarchical classification
getBetween = function(t, k, c=NULL, d=NULL) {
  if (t == 2) return (getCNH(t, d, k)$betweenss)
  # get a vector of heights for the k first groups
  # sum of squares of the vector (height is a distanc, an inertia is dist^2)
  # begin to 2 clusters for CAH
  if (t > 2 ) sum(c$height[length(c$height):(length(c$height)-k+2)]^2)
}

# Same tot no matter k, so k =2 for CNH
getTotInertia = function(t, c, d) {
  #(height is a distanc, an inertia is dist^2)
  if(t > 2 ) sum(c$height^2)
  else if (t == 2) getCNH(t, d, 2)$totss
}

# Raw inter-group inertia for each partitionning
# Inputs:
# t: number of type of classification
# n: maximum number of clusters
# c: hierarchical classification
# d: dataframe
# Output: between-group inertia for all clusters
getBetweenPerPart = function(t, n, c=NULL, d=NULL) {
  #TODO: warnings if t< 2  and n = nrow(d)
  #TODO: warnings if t> 2 and no classif
  #warnings if t <= 2 and no data
  between = vector(mode="numeric", n-1)
  for (k in 2:n)
    between[k-1] = getBetween(t, k, c, d)
  return (between)
}

getRelativeBetweenPerPart = function(t, n, c=NULL, d=NULL){
  100*((getBetweenPerPart(t, n, c, d) / getTotInertia(t, c, d)))
}

getBetweenDifferences = function(t, n, c=NULL, d=NULL){
  between = round(getRelativeBetweenPerPart(t, n, c, d),2)
  diff = rep(0,length(between))
  #The difference between a uniq cluster and 2 is by default, the first inertia value
  diff[1] = between[1]
  for (i in 2:length(between)) diff[i] = between[i] - between[i-1]
  return(diff)
}

# Between inertia differences between a partionning and the previous
plotFusionLevels = function(t, n, c=NULL, d=NULL) {
  if (v==T) cat("\nBETWEEN DIFFERENCES:\n")
  between_diff = getBetweenDifferences(t, n, c, d)
  
  #x11()
  optimal_nb_clusters = which.max(between_diff)+1
  savePdf("fusion_levels.pdf")
  plot(2:n, between_diff, type="b", ylim=c(round(min(between_diff))-1,round(max(between_diff))+1), xlim=c(2,n+1), xlab="Nb. of clusters", ylab="Between-cluster differences (%)", col="grey", axes=F)
  printBestClustering("Inertia gap method", between_diff, " difference with the previous partitionning (%)", optimal_nb_clusters)
  suprLog = dev.off()
}

################################
#          Silhouette
################################

#Ouput: an ordered silhouette object
getSilhouette = function(t, k , c, d){
  clusters = getClusters(t, k , c, d)
  diss = getDistance(d,t,k)
  sil = sortSilhouette(silhouette(clusters, diss))
  rownames(sil) = row.names(d)[attr(sil,"iOrd")]
  return (sil)
}

getSilhouettePerPart =function(t, n, c=NULL, d=NULL){
  mean_silhouette = numeric(n - 1)
  for (k in 2:(n - 1)) {
    si = getSilhouette(t, k , c, d)
    mean_silhouette[k] = summary(si)$avg.width
  }
  return(mean_silhouette[-1])
}

# Plot the best average silhouette width for all clustering possible
plotSilhouettePerPart = function(t, n, c=NULL, d=NULL){
  if (v==T) cat("\nSILHOUETTE:\n")
  mean_silhouette = getSilhouettePerPart(t, n, c, d)
  
  #x11()
  savePdf("average_silhouettes.pdf")
  optimal_nb_clusters = which.max(mean_silhouette)+1
  plot(2:(n-1), mean_silhouette, type="b", xlim=c(2,n), ylim=c(0,max(mean_silhouette)+0.1), col="grey", xlab="Nb. of clusters", ylab="Average silhouette width", axes=F)
  printBestClustering("Silhouette method", mean_silhouette,"n average silhouette width", optimal_nb_clusters, 0.1)
  suprLog = dev.off()
  return (optimal_nb_clusters)
}

#TODO: here: setParam
plotSilhouette = function(s){
  #x11()
  pdf("silhouette.pdf")
  setGraphicBasic()
  par(mar=c(4, 12, 3, 2))
  plot(s, max.strlen=25, main=" ", sub= "", do.clus.stat=TRUE, xlab="Silhouette width", cex.names=0.8, col=colorClusters(s[,1]), nmax.lab=100, do.n.k = FALSE, axes=F)
  mtext(paste("Average silhouette width:", round(summary(s)$avg.width,3)), font=2, cex=1.5, line=1)
  printAxis(1, 0, 1, 0.2)
  suprLog = dev.off()
}

printSummary = function(t, n, c=NULL, d=NULL){ 
  #TODO: no n = nrow(data)
  between = getRelativeBetweenPerPart(t, n, c, d)
  summary = cbind(between, getBetweenDifferences(t, n, c, d), 100-between, getSilhouettePerPart(t,n+1,c,d))
  rownames(summary) = seq(2, n) 
  colnames(summary) = c("Between-inertia (%)", "Between-differences (%)", "Within-inertia (%)", "Silhouette index") 
  return (summary)
}

################################
#          HEATMAP
################################

#Inputs:
# s: an organised silhouette object
printRect = function (s){
  # size of each clusters
  cl_sizes = summary(s)$clus.sizes
  tempSize = 0
  #vector of color: one by cluster
  colors = colPers(length(cl_sizes))
  for (i in 1:length(cl_sizes)){
    #y begin at the top, so sum(cl_sizes) must be substracted to y coord.
    #rect(xleft, ybottom, xright, ytop)
    # +0.5 because x, y coord are shifted to 0.5 comparativly to plotcolors functions
    rect(tempSize + 0.5, sum(cl_sizes) -tempSize -cl_sizes[i] +0.5, cl_sizes[i] +tempSize +0.5, sum(cl_sizes) -tempSize +0.5, border = colors[i], lwd=3)
    #memorize the size of the cluster (for a bottom-right shift)
    tempSize = cl_sizes[i]
  }
}

#Inputs:
# d: a distance object
# s: an organised silhouette object
heatMap = function(d, s, text=FALSE){
  #x11()
  pdf("heat_map.pdf")
  options(warn = -1)
  matrix=as.matrix(d)
  matrix=matrix[attr(s,"iOrd"),attr(s,"iOrd")]
  rownames(matrix) = rownames(data)[attr(s,"iOrd")]
  labels = attr(d, "Labels")[attr(s,"iOrd")]
  #if(tri == TRUE) matrix[!lower.tri(matrix)] = NA
  #image(1:ncol(matrix), 1:ncol(matrix), t(matrix), axes=F, xlab="", ylab="")

  par(mar=c(1, 8, 8, 1))
  par(fig=c(0,0.9,0,1), new=TRUE)
  plotcolors(dmat.color(matrix, colors=heat.colors(1000),byrank = FALSE), ptype="image", na.color="red", rlabels=FALSE, clabels=FALSE, border=0)
  
  mtext('Distance matrix ordered by silhouette\'s scores', 3, line=6, font=4, cex=1.5)
  
  text(-0.5, 0:(ncol(matrix)-1)+1, rev(labels), xpd=NA, adj=1, cex=0.7)
  text(0.5:(ncol(matrix)-0.5), ncol(matrix)+1, substr(labels, 0, 20), xpd=NA, cex=0.7, srt=65, pos=4)
  printRect(s)
  
  if (text==TRUE)   text(expand.grid(1:ncol(matrix), ncol(matrix):1), sprintf("%d", matrix), cex=0.4)
  #axis(2, 1:ncol(matrix), labels, cex.axis = 0.5, las=1, tck=0, lwd=-1, font.axis=3)
  #axis(3, 1:ncol(matrix), labels, cex.axis = 0.5, las=1, tck=0, lwd=-1, font.axis=3)
  
  par(mar=c(5, 0, 4, 0) + 0.1)
  par(fig=c(0.85,1,0.3,0.8),new=TRUE)
  legend_image = as.raster(matrix(heat.colors(1000), ncol=1))
  plot(c(0,1),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = '')
  rasterImage(legend_image, 0.4, 0, 0.5, 1)
  mtext('   Distance', 3, line=0.5, cex=0.85, font=2)
  text(x=0.5, y = seq(0,1,l=3), labels = round(seq(max(matrix),0,l=3)),cex=0.7,pos=4)
  options(warn = 0)
  suprLog = dev.off()
}

################################
#          Dendrogram
################################

# Inputs:
# k: number of clusters
plotDendrogram = function(t, k, c, d){
  #x11()
  pdf("dendrogram.pdf")
  setGraphicBasic()
  par(mar=c(2,5,5,1))
  #having relative inertia instead of raw cophenetic distance
  #rev() beacause in cah function, the vector is inversed
  #c$height = rev(getBetweenDifferences(t, nrow(d), c, d))
  plot(c, ylim=c(0,max(c$height)), xlim=c(0,length(c$labels)), hang=-1, sub="", cex=0.8, font=3, ylab="Between-cluster differences (%)", main="Dendrogram", axes=F)
  #text(-0.5, 0:(ncol(matrix)-1)+1, rev(labels), xpd=NA, adj=1, cex=0.7)
  printAxis(2, 0, max(c$height))
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
  par(mar=c(0,0,4.1,0))
  clusters = getClusters(t, k, c, d)
  title = paste("Cumulated inertia:", round((pca$eig[1]+pca$eig[2])/sum(pca$eig),4)*100, "%")
  s.class(addaxes=F, pca$li ,ylim=c(min(pca$li[,2])-3, max(pca$li[,2])+3), xlim=c(min(pca$li[,1])-3, max(pca$li[,1])+3), csub=1.5, as.factor(clusters), grid=F, col=colPers(optimal_nb_clusters))
  mtext(title, font=2, cex=1.5, line=1)
  abline(h=0, v=0, lty=2, lwd=2, col="grey")
  text(x=pca$li[,1], y=pca$li[,2], labels=rownames(pca$li), col=colorClusters(clusters), cex=0.6)
  suprLog = dev.off()
}

################################
#            CTR
################################

# For a given partition, contribution of each metabolites to each clusters
# Inputs:
# d: data
# cl: clusters object
getBetweenPerVariable = function(d, cl){
  #Get the normalized distance between each points and the center
  # normalized to get percent values in output
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

# For a given partition, relative contributions of each metabolites to inertia of each clusters (CTR)
# The total of the clusters for each column corresponds to PDIS
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
# Relative contributions of the metabolites to inertia of a partitionning (in %)
# Inputs: 
# t: number of type of classification
# k: number of clusters
# c: hierarchical classification
# d: data
getPdis = function(t, k, c, d) {
  
  #for each metabolite contribution (in column), sum the k clusters values
  return(apply(getCtrVar(t, k, c, d), 2, sum))
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
# Rho=0, a cluster close to the center
# Rho= infinity, a cluster far from the center
# Inputs: 
# t: number of type of classification
# k: number of clusters
# c: hierarchical classification
# d: data
getRho2 = function(t, k, c, d) {
  
  cl = getClusters(t, k, c, d)
  nb_cl = length(levels(as.factor(cl)))
  nb_met = length(cl)
  ctr = getBetweenPerVariable(d, cl)
  
  for (i in 1:k)
    #each metabolites is weighted by the total number of metabolites in each clusters
    for (j in 1:nb_met) ctr[i,j] = ctr[i,j]/length(cl[cl==i])
  
  #sum of row ^2 (a row =  a cluster)
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
librairies = c("cluster", "optparse", "gclus", "ade4", "scales", "gplots")
for (l in librairies){
  if (! (l %in% installed.packages()[,"Package"])) install.packages(l, repos = "http://cran.us.r-project.org")
  library(l, character.only = TRUE)
}

#Get arguments
args = getArgs()
checkArg(args)
opt = parse_args(args)

#Global variables settings
nb_clusters = opt$nbCluster
max_cluster = opt$maxCluster
typeClassif = opt$typeClassif
advanced = "advanced" %in% names(opt)
v = !("quiet" %in% names(opt))
ranked = !("ranked" %in% names(opt))
if (!is.null(opt$workdir)) setwd(opt$workdir)
#to work under Rstudio
#setwd("~/bin/fingerprint_clustering/")

#Loading data
data = read.table(opt$infile, header=F, sep="\t", dec=".", row.names=1)
colnames(data) <- substr(rownames(data), 1, 25) -> rownames(data)
postChecking(args, data)

#Perform classification
classif = getCAH(data, typeClassif)
if(typeClassif>2) plotCohenetic(typeClassif, data, classif)
if(typeClassif>1) plotFusionLevels(typeClassif, max_cluster, classif, data)

#Silhouette analysis
optimal_nb_clusters = plotSilhouettePerPart(typeClassif, max_cluster + 1, classif, data)
if(!is.null(nb_clusters)) optimal_nb_clusters = nb_clusters
sil = getSilhouette(typeClassif, optimal_nb_clusters, classif, data)
plotSilhouette(sil)
if(typeClassif > 1){ 
  summary = round(printSummary(typeClassif, max_cluster, classif, data),2)
  writeTsv("summary")
}

#Global variables settings
dis = getDistance(data, typeClassif, optimal_nb_clusters)
clusters = getClusters(typeClassif, optimal_nb_clusters, classif, data)

#Plots
if(typeClassif > 2) plotDendrogram(typeClassif, optimal_nb_clusters, classif, data)
plotPca(typeClassif, optimal_nb_clusters, classif, data)

#Advanced indexes
if (advanced == TRUE){
  
  contribution = round(1000 * getCtrVar(typeClassif, optimal_nb_clusters, classif, data) / 10)
  discriminant_power = getIndexPerPartition(typeClassif, max_cluster, classif, data, "pdis")
  excentricity = getIndexPerPartition(typeClassif, max_cluster, classif, data, "rho2")
  
  for (i in c("contribution", "discriminant_power", "excentricity"))
    writeTsv(i)
  

}
if(typeClassif <= 2 | advanced == TRUE){
  heatMap(dis, sil, TRUE)
}else{
  pdf("heat_map.pdf")
  heatmap.2(as.matrix(data),
            distfun = function(x) dist(x,method = 'euclidean'),
            hclustfun = function(x) hclust(x,method = getTypeClassif(typeClassif)),
            #density.info="none",
            trace="none",
            #lhei = c(0.05,0.95),
            cexRow=0.6,
            cexCol=0.6,
            key=TRUE,
            key.xlab="Distance",
            key.ylab="",
            #keysize=1,
            key.title="",
            margins=c(6,8),
            revC=TRUE,
            srtCol=45,
            symm=TRUE,
            dendrogram="row",
            #col=heat.colors(1000)
            )
  title(main="Distance matrix\n ordered by CAH's distances", cex.main=1.7, line=-1)
  suprLog = dev.off()
}

#Final outputs
writeClusters(clusters, ranked)
if (v != T) cat(paste("Clustering done.\nOptimal number of clusters choosen:", optimal_nb_clusters,"\n"))