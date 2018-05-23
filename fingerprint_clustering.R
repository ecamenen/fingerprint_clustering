getArgs = function(){
  option_list = list(
    make_option(c("-w", "--workdir"), type="character", metavar="character",
                help="Working directory path [default: the folder where the script is launched]"),
    make_option(c("-i", "--infile"), type="character", default="matrix.txt", 
                metavar="character",
                help="Fingerprint file name [default: %default]"),
    make_option(c("-m", "--maxCluster"), type="integer", default=6, metavar="integer",
                help="Maximum number of clusters [default: %default]"),
    make_option(c("-t", "--classif_type"), type="integer", default=4, metavar="integer",
                help="Type of classifation [default: automatic selection of best CAH] (1: K-menoids; 2: K-means; 3: Ward; 4: Complete links; 5: Single links; 6: UPGMA; 7: WPGMA; 8: WPGMC; 9: UPGMC)"),
    make_option(c("-adv", "--advanced"), type="logical", action="store_true", 
                help="Activate advanced mode (print more outputs)"),
    make_option(c("-q", "--quiet"), type="logical", action="store_true",
                help="Activate quiet mode"),
    make_option(c("-T", "--text"), type="logical", action="store_true",
                help="DO NOT print values on graph"), 
    make_option(c("-n", "--nbCluster"), type="integer", metavar="integer",
                help="Fix the number of clusters"),
    make_option(c("-r", "--ranked"), type="logical", action="store_true", 
                help="Rank the metabolites in clusters by silhouette scores instead of alphabetically"),
    make_option(c("-b", "--bootstrap"), type="integer", default=500, metavar="integer",
                help="Number of bootstraps for Gap statistic (advanced mode)")
    
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
  if(opt$bootstrap < 100 || opt$bootstrap > 1000){
    print_help(a)
    stop("--bootstrap comprise between 100 and 1000", call.=FALSE)
  }
  
  checkMinCluster = function (o, def="")
  if (opt[[o]] < 2){
    print_help(a)
    stop(paste("--",o ," must be upper or equal to 2",def,".\n",sep=""), call.=FALSE)
  }
  checkMinCluster("maxCluster"," [by default: 6]")
  if(!is.null(opt$nbCluster)) checkMinCluster("nbCluster")
  
  if ((opt$classif_type < 1) || (opt$classif_type > 9)){
    print_help(a)
    stop("--classif_type must be comprise between 1 and 6 [by default: 2].\n", call.=FALSE)
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
  # ponderation for sampling index (var use n-1)
  # without this constante, for advanced outputs, total (max_cluster=nrow(data)) will be different from 1
}

getDistance = function(d, t, k=NULL){
  if (t > 1) dist(d, method = "euclidian")
  else getCNH(t,d,k)$diss
}

#Inputs: x : a matrix
#filename of the saved file
#Prints the matrix, save the matrix
writeTsv = function(x, cl=TRUE){
  #print on stdout
  if (isTRUE(verbose)) cat(paste("\n", gsub("_", " ", toupper(x)), ":\n", sep=""))
  #disabling warning
  options(warn = -1)
  #get variable
  tab = get(x)
  if(isTRUE(cl)) output=as.matrix(rbind(c("", colnames(tab)), cbind(rownames(tab),tab)))
  else output = tab
  #discard empty rows
  output = output[rowSums(is.na(output)) != ncol(output),]
  #TODOD:
  #output = output[,colSums(is.na(output)) != nrow(output)]
  output[is.na(output)] = ""
  colnames(output)=rep("", ncol(output)); rownames(output)=rep("", nrow(output))
  if (isTRUE(verbose)){
    if (isTRUE(cl)){
    printed = round(apply(output[-1,-1],2,as.numeric),2)
    rownames(printed) = rownames(tab)
    colnames(printed) = colnames(tab)
    }else{
      printed = output
    }
    print(printed,  quote=F)
  }
  write(t(output), paste(x,".tsv",sep=""), ncolumns=ncol(output), sep="\t")
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

plotAxis = function (side, min, max, interval = 1){
  axis(side, seq(min,max, interval), lwd=3)
}

plotBestClustering = function(sub_title, values, values_type, optimal_nb_clusters, interval = 1, min_x=2, best=NULL){
  plotAxis(1, 2, max_cluster)
  if (interval >= 1) axisSeq=round(values)
  else axisSeq = c(0, max(values) +0.1)
  #case of plotting gap statistics
  if (min_x < 2) best_y=values[optimal_nb_clusters]
  else best_y = max(values)
  #for non-elbow plots
  if (is.null(best)) best = round(max(values),2)
  plotAxis(2, min(axisSeq), max(axisSeq), interval)
  title(main="Optimal number of clusters", line=2, cex.main=2)
  mtext(text=sub_title, font=3, cex=1.2, line = 0.5)
  abline(v=optimal_nb_clusters, col="red", lty=2, lwd=2)
  points(optimal_nb_clusters, best_y, pch=19, col="red", cex=2)
  if(isTRUE(text)) text(y=values, x=min_x:max_cluster, labels=round(values,2), cex=1.2, pos=4, col="red")
  if (isTRUE(verbose)) cat("Optimal number of clusters k = ", optimal_nb_clusters, "\n","With a", values_type, " of ", best, "\n", sep="")
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
    dis = getDistance(d, t)
    #cah: classification hierarchic ascending
    cah = hclust(dis, method=getClassifType(t))
  #automaticly ordering by clusters
  return (reorder.hclust(cah, d))
  }
}

# Selects best algo based on cophenetic calculation
# d: data (or distance matrix for hierarchic)
selectBestCAH = function (d, v=F){
  dis = getDistance(d, 2)
  temp = 0
  for (i in 3:9){
    cah = getCAH(data, i)
    res = cor(dis, cophenetic(cah))
    if (isTRUE(v)) cat(paste(getClassifType(i), ":",round(res,3), "\n"))
    if (res > temp){ 
      temp = res
      t = i
    }
  }
  cat(paste("Selected:", getClassifType(t),"\n"))
  return (t)
}

#Inputs:
# t: number of type of classification
getClassifType = function(t){
  if (t==3) "ward.D2"
  else if (t==4) "complete"
  else if (t==5) "single"
  else if (t==6) "average"
  else if (t==7) "mcquitty"
  else if (t==8) "median"
  else if (t==9) "centroid"
}

#Agglomerative coefficient ()
getCoefAggl = function(c)
  coef.hclust(c)

#Inputs: 
# t: number of type of classification
# d: data (or distance for pam)
# k: number of clusterting
#Ouput: Non-hierarchical classification
getCNH = function(t, d, k){
  if (t==1) return (pam(d, k))
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
  writeTsv("clusters", F)
}

############################################################
#          Cophenetic (dendrogram distance matrix)
############################################################

# Distance matrix between each leaf of the dendogramm
#Inputs:
# d : data
# cah : hierarchical classification
plotCohenetic=function(t, d, cah){
  dis = getDistance(d, t)
  coph_matrix = cophenetic(cah)
  cor_coph = cor(dis, coph_matrix)
  if (isTRUE(verbose)) cat(paste("\nCOPHENETIC:\nExplained variance (%):", round(cor_coph^2,3), "\nCorrelation with the data:",round(cor_coph,3),"\n"))

  savePdf("shepard_graph.pdf")
  plot(dis, coph_matrix, pch=19, col=alpha("red",0.2), axes=F, xlim=c(0,max(dis)), ylim=c(0,max(coph_matrix)), xlab="Distance between metabolites",ylab="Cophenetic distance", asp=1, main=paste("Cophenetic correlation: ",round(cor_coph,3)))
  plotAxis(2, 0, max(coph_matrix))
  plotAxis(1, 0, max(dis))
  abline(0, 1, col="grey", lwd=3, lty=2)
  suprLog = dev.off()
}

##############################################
#          Inertia
##############################################

# Relative inter-group inertia for each partitionning
# Inputs:
# t: number of type of classification
# n: maximum number of clusters
# c: hierarchical classification
# d: dataframe
getRelativeBetweenPerPart = function(t, n, c = NULL, d = NULL){
  d=as.matrix(d)
  between = rep(0, n-1)
  # total sum of square
  TSS = sum(scale(d, scale = FALSE)^2)
  for (i in 2:n) {
    cl = as.factor(getClusters(t, i, c, d))
    # tapply(data[,i], Cla, mean) :
    # centroids of each clusters for a column i
    # sapply(1:ncol(data), function(i) tapply(data[,i], Cla, mean)) :
    # centroids of each clusters for each column
    # apply(d, 2, mean) : centroids for each column
    # as.vector(table(cl) : size of each clusters
    # t : vector rotation for arithmetic with other row or column vectors
    between[i-1] = sum(t((t(sapply(1:ncol(d), function(i) tapply(d[,i], cl, mean)))-
                            apply(d, 2, mean))^2) * as.vector(table(cl)))/TSS
  }
  return (100*between)
}

# Difference between each case of a vector
getBetweenDifferences = function(t, n, c=NULL, d=NULL){
  between = getRelativeBetweenPerPart(t, n, c, d)
  # apply produce a list, unlist convert in vector
  diff = unlist(sapply(1:n, function(i) between[i]-between[i-1]))
  return (as.vector(cbind(between[1], t(diff[-(n-1)]))))
  #-n-1 to remove the last NA value (pairwise comparison)
  #between[1] to get the difference with 1 cluster
}

# Between inertia differences between a partionning and the previous
plotFusionLevels = function(t, n, c=NULL, d=NULL) {
  if (isTRUE(verbose)) cat("\nBETWEEN DIFFERENCES:\n")
  between_diff = getBetweenDifferences(t, n, c, d)
  #between_diff = 100 - getRelativeBetweenPerPart(t, n, c, d)
  #relative.loss = intra[2:(nrow(data)-1)]/intra[1:(nrow(data) - 2)]
  
  optimal_nb_clusters = which.max(between_diff)+1
  savePdf("between_differences.pdf")
  plot(2:n, between_diff, type="b", ylim=c(round(min(between_diff))-1,round(max(between_diff))+1), xlim=c(2,n+1), xlab="Nb. of clusters", ylab="Between-cluster variation (%)", col="grey", axes=F)
  plotBestClustering("Fusion level method", between_diff, " variation with the previous partitionning (%)", optimal_nb_clusters)
  suprLog = dev.off()
}

plotElbow = function(t, n, c=NULL, d=NULL) {
  if (isTRUE(verbose)) cat("\nELBOW:\n")
  within = c(100, 100 - getRelativeBetweenPerPart(t, n, c, d))
  ratio = within[1:(n-1)] / within[2:n]
  best = round(min(ratio),2); optimal_nb_clusters = which.min(ratio)
  savePdf("elbow.pdf")
  plot(1:n, within, type="b", ylim=c(-1,101), xlim=c(1,n+1), xlab="Nb. of clusters", ylab="Relative within inertia (%)", col="grey", axes=F)
  plotBestClustering("Elbow method", within, " Wk/Wk+1 ratio ", optimal_nb_clusters, 5, 1, best)
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
  if (isTRUE(verbose)) cat("\nSILHOUETTE:\n")
  mean_silhouette = getSilhouettePerPart(t, n, c, d)
  
  savePdf("average_silhouettes.pdf")
  optimal_nb_clusters = which.max(mean_silhouette)+1
  plot(2:(n-1), mean_silhouette, type="b", xlim=c(2,n), ylim=c(0,max(mean_silhouette)+0.1), col="grey", xlab="Nb. of clusters", ylab="Average silhouette width", axes=F)
  plotBestClustering("Silhouette method", mean_silhouette,"n average width", optimal_nb_clusters, 0.1)
  suprLog = dev.off()
  return (optimal_nb_clusters)
}

#TODO: here: setParam
plotSilhouette = function(s){
  pdf("silhouette.pdf")
  setGraphicBasic()
  par(mar=c(4, 12, 3, 2))
  plot(s, max.strlen=25, main=" ", sub= "", do.clus.stat=TRUE, xlab="Silhouette width", cex.names=0.8, col=colorClusters(s[,1]), nmax.lab=100, do.n.k = FALSE, axes=F)
  mtext(paste("Average silhouette width:", round(summary(s)$avg.width,3)), font=2, cex=1.5, line=1)
  plotAxis(1, 0, 1, 0.2)
  suprLog = dev.off()
}

###################################
#          GAP STATISTICS
###################################

#B: nb of bootstrap
getGapPerPart = function(d, n, B=500){
  #FUN mus have only two args in this order and return a list with an object cluster
  if(classif_type>2) gapFun = function(x,k) list(cluster = getClusters(classif_type, k, getCAH(x,classif_type)))
  else gapFun = function(x,k) getCNH(classif_type,x,k)
  clusGap(d,FUN=gapFun,K.max=n, verbose=F, B=B)
}

#g: gap object
getGapBest = function (g, M="Tibs2001SEmax"){
  with(g, maxSE(Tab[,"gap"], Tab[,"SE.sim"], method=M))
}

# Plot the gap statistics width for all clustering possible
plotGapPerPart = function(t, n, d, B){
  if (isTRUE(verbose)) cat("\nGAP STATISTICS:\n")
  if(t==2 & (max_cluster > 10 | nrow(d)>=100)) plural=c("few ", "s")
  else plural=c("","")
  if(t==2 | nrow(d)>=100 ) cat(paste("It could take a ",plural[1], "minute",plural[2],"...\n",sep=""))
  gap = getGapPerPart(d, n, B)
  savePdf("gap_statistics.pdf")
  optimal_nb_clusters = getGapBest(gap)
  plot(gap, arrowArgs = list(col="gray", length=1/15, lwd=2, angle=90, code=3), type="b", xlim=c(1,n+1), ylim=c(0,max(gap$Tab[,"gap"])+0.1), col="grey", xlab="Nb. of clusters", ylab=expression(Gap[k]), main="",axes=F)
  plotBestClustering("Gap statistics method", gap$Tab[,"gap"]," gap value", optimal_nb_clusters, 0.1, 1)
  #cat(paste("With a corrected index, optimal number of clusters k =",getGapBest(gap,"firstSEmax"), "\n"))
  suprLog = dev.off()
  return (gap)
}

#Plot the gap between the two function: within and random within average
plotGapPerPart2 = function(g, n){
  savePdf("gap_statistics2.pdf")
  min_y=round(min(g$Tab[,c(1,2)]),1)
  max_y=round(max(g$Tab[,c(1,2)]),1)
  plot(0,0, xlim=c(1,n), ylim=c(min_y-0.1,max_y+0.1),type="n", xlab="Nb. of clusters", ylab="log(within-inertia)", axes=F)
  title(main="Optimal number of clusters", line=2, cex.main=2)
  mtext(text="Gap statistics method", font=3, cex=1.2, line = 0.5)
  optimal_nb_clusters = getGapBest(g)
  abline(v=optimal_nb_clusters, col="gray", lty=2, lwd=2)
  lines(seq(1:n),g$Tab[,1],type="b", col="red")
  lines(seq(1:n),g$Tab[,2],type="b", col="blue")
  plotAxis(1,1,n)
  plotAxis(2,min_y, max_y,0.1)
  legend("topright",c("log(W)", "E.log(W)"), col=c("red","blue"), lty=1, box.lwd=-1, bg = "white")
  suprLog = dev.off()
}

printSummary = function(t, n, adv=F, c=NULL, d=NULL){ 
  #TODO: no n = nrow(data)
  between = getRelativeBetweenPerPart(t, n, c, d)
  summary = cbind(between, getBetweenDifferences(t, n, c, d), 100- getRelativeBetweenPerPart(t,n,classif, data), getSilhouettePerPart(t,n+1,c,d))
  names = c("Between-inertia (%)", "Between-differences (%)", "Within-inertia (%)", "Silhouette index") 
  if(isTRUE(adv)) {
    summary = cbind(summary, gap$Tab[,"gap"][-1], gap$Tab[,"SE.sim"][-1])
    names=c(names, "Gap statistics", "Gap SE")
  }
  rownames(summary) = seq(2, n) 
  colnames(summary) = names  
  return (summary)
}

################################
#          HEATMAP
################################

#Inputs:
# cl_size: vector of size for each clusters
plotRect = function (cl_sizes, colors){
  # size of each clusters
  temp_size = 0
  for (i in 1:length(cl_sizes)){
    #y begin at the top, so sum(cl_sizes) must be substracted to y coord.
    #rect(xleft, ybottom, xright, ytop)
    # +0.5 because x, y coord are shifted to 0.5 comparativly to plotcolors functions
    rect(temp_size + 0.5, sum(cl_sizes) -temp_size -cl_sizes[i] +0.5, cl_sizes[i] +temp_size +0.5, sum(cl_sizes) -temp_size +0.5, border = colors[i], lwd=3)
    #memorize the size of the cluster (for a bottom-right shift)
    temp_size = temp_size + cl_sizes[i]
  }
}

#Outputs:
# lenght of clusters ordered by the clusters order
getOrderedClusterSize = function(cl){
  nb_cl =  length(levels(as.factor(cl))) 
  size_cl = rep(0, nb_cl)
  temp_cl = rep(0, length(cl))
  j = 0
  
  for (i in 1:length(cl)) {
    if (!cl[i] %in% temp_cl) j = j+1
    size_cl[j] = size_cl[j] + 1
    temp_cl[i] = cl[i]
  }
  return (size_cl)
}

#Inputs:
# d: a distance object
# s: an organised silhouette object
# c: CAH
# c: clusters from CAH
heatMap = function(d, s=NULL, c=NULL, cl=NULL, text=FALSE){
  
  if(!is.null(s)){
    order = attr(s,"iOrd")
    cl_sizes = summary(s)$clus.size
    title = "silhouette\'s scores"
    colors = colPers(length(levels(as.factor(sil[,1]))))
  }else{
    order = c$order
    cl_sizes = getOrderedClusterSize(cl[order])
    title="dendrogram"
    colors = orderColors(c, cl)
  }

  matrix=as.matrix(d)
  matrix=matrix[order, order]
  rownames(matrix) <- rownames(d)[order] -> labels
  #if(tri == TRUE) matrix[!lower.tri(matrix)] = NA
  #image(1:ncol(matrix), 1:ncol(matrix), t(matrix), axes=F, xlab="", ylab="")

  options(warn = -1)
  pdf("heat_map.pdf")
  
  par(fig=c(0,0.9,0,1), new=TRUE)
  par(mar=c(1, 8, 8, 1))
  plotcolors(dmat.color(matrix, colors=heat.colors(1000),byrank = FALSE), ptype="image", na.color="red", rlabels=FALSE, clabels=FALSE, border=0)
  mtext(paste('Distance matrix ordered by', title), 3, line=6, font=4, cex=1.5)
  text(-0.5, 0:(ncol(matrix)-1)+1, rev(labels), xpd=NA, adj=1, cex=0.7)
  text(0.5:(ncol(matrix)-0.5), ncol(matrix)+1, substr(labels, 0, 20), xpd=NA, cex=0.7, srt=65, pos=4)
  plotRect(cl_sizes, colors)
  if (isTRUE(text)) text(expand.grid(1:ncol(matrix), ncol(matrix):1), sprintf("%d", matrix), cex=0.4)

  par(fig=c(0.85,1,0.3,0.8),new=TRUE)
  par(mar=c(5, 0, 4, 0) + 0.1)
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
plotDendrogram = function(t, k, c, d, adv=FALSE){

  pdf("dendrogram.pdf")
  setGraphicBasic()
  par(mar=c(2,5,5,1))
  plot(c, hang=-1, ylim=c(0,max(c$height)), xlim=c(0,length(c$labels)), sub="", cex=0.8, font=3, ylab="Cophenetic distance", main="Dendrogram", axes=F)
  plotAxis(2, 0, max(c$height))
  #projection of the clusters
  rect.hclust(c, k=as.numeric(k), border=orderColors(c, clusters))
  suprLog = dev.off()
}

# Get colors ordered for dendrogram
orderColors = function(c, cl){
  col_in = colorClusters(cl)[c$order]
  j = 1
  col_ordered = rep(NA, length(table(clusters)))
  col_ordered[1] = col_in[1]
  for (i in 2:length(col_in)){
    if (col_in[i] != col_in[i-1]){
      j = j + 1
      col_ordered[j] = col_in[i]
    }
  }
  #vector of color: one by cluster
  return (col_ordered)
}

################################
#            PCA
################################

plotPca = function(t, k, c, d){
  pca = dudi.pca(d, scannf=F)
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

#########################################
#            Variables contribution
#########################################

#For a given partition (cl) and each variables (dataset columns)
#pondered distance between the centroid of each clusters and the global centroid of the cloud 
# Inputs:
# d: data
# cl: clusters object
getDistPerVariable = function(d, cl){
  #Distance between the centroid of each variables 
  #ponderation by the sd of the variable (=total inertia per var)
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
  
  ctr = getDistPerVariable(d, cl)
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
getPdisPerPartition = function(t, n, c, d){
  
  pdis_per_partition = matrix(NA, n-1, ncol(d))
  rownames(pdis_per_partition) = seq(2, n)
  
  for (k in 2:n){
    colnames(pdis_per_partition) = colnames(d)
    res = getPdis(t, k, c, d)
    for(i in 1:length(res)){
      pdis_per_partition[k-1, i] = res[i]
    }
  }
  return (pdis_per_partition)
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
  if (! (l %in% installed.packages()[,"Package"])) install.packages(l, repos = "http://cran.us.r-project.org", quiet = T)
  library(l, character.only = TRUE)
}

#Get arguments
args = getArgs()
checkArg(args)
opt = parse_args(args)

#Global variables settings
nb_clusters = opt$nbCluster
max_cluster = opt$maxCluster
classif_type = opt$classif_type
bootstrap = opt$bootstrap
advanced = "advanced" %in% names(opt)
verbose= !("quiet" %in% names(opt))
ranked = !("ranked" %in% names(opt))
text = !("text" %in% names(opt))
if (!is.null(opt$workdir)) setwd(opt$workdir)

#Loading data
data = read.table(opt$infile, header=F, sep="\t", dec=".", row.names=1)
colnames(data) <- substr(rownames(data), 1, 25) -> rownames(data)
postChecking(args, data)

#Perform classification
classif = getCAH(data, classif_type)
if(classif_type>2){
  plotCohenetic(classif_type, data, classif)
  if(isTRUE(advanced)) cat(paste("\nAGGLOMERATIVE COEFFICIENT: ", round(getCoefAggl(classif),3), "\n", sep=""))
}

plotFusionLevels(classif_type, max_cluster, classif, data)

#Silhouette analysis
optimal_nb_clusters = plotSilhouettePerPart(classif_type, max_cluster + 1, classif, data)
if(!is.null(nb_clusters)) optimal_nb_clusters = nb_clusters
sil = getSilhouette(classif_type, optimal_nb_clusters, classif, data)
plotSilhouette(sil)

#Global variables settings
dis = getDistance(data, classif_type, optimal_nb_clusters)
#cl_temp, because writeTsv(clusters) recreate an object
clusters <- getClusters(classif_type, optimal_nb_clusters, classif, data) -> cl_temp

#Advanced indexes
if (isTRUE(advanced)){
  
  elbow_k = plotElbow(classif_type, max_cluster, classif, data)
  gap = plotGapPerPart(classif_type, max_cluster, data, bootstrap)
  plotGapPerPart2(gap, max_cluster)
  
  contribution = 100 * getCtrVar(classif_type, optimal_nb_clusters, classif, data)
  discriminant_power = 100 * getPdisPerPartition(classif_type, max_cluster, classif, data)
  
  verbose = F
  for (i in c("contribution", "discriminant_power"))
    writeTsv(i)
  verbose = !("quiet" %in% names(opt))
}

#Plots
if(classif_type > 2) plotDendrogram(classif_type, optimal_nb_clusters, classif, data, advanced)
plotPca(classif_type, optimal_nb_clusters, classif, data)
if(classif_type <= 2 || isTRUE(advanced)){
  heatMap(data, sil, text=T)
}else{
  heatMap(data, c=classif, cl=clusters, text=T)
}

#Final outputs
summary = printSummary(classif_type, max_cluster, advanced, classif, data)
writeTsv("summary")
writeClusters(clusters, ranked)
if (!isTRUE(verbose)) cat(paste("Optimal number of clusters:", optimal_nb_clusters,"\n"))

#errors
if (optimal_nb_clusters==max_cluster) message("\n[WARNING] The optimal number of clusters equals the maximum number of clusters. \nNo cluster structure has been found.")
if(min(table(cl_temp))==1) message("\n[WARNING] A cluster with an only singleton biased the silhouette score.")
