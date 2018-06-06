getArgs = function(){
  option_list = list(
    make_option(c("-w", "--workdir"), type="character", metavar="character",
                help="Working directory path [default: the folder where the script is launched]"),
    make_option(c("-i", "--infile"), type="character", default="data/matrix.txt", 
                metavar="character",
                help="Fingerprint file name [default: %default]"),
    make_option(c("-m", "--max_clusters"), type="integer", default=6, metavar="integer",
                help="Maximum number of clusters [default: %default]"),
    make_option(c("-t", "--classif_type"), type="integer", default=4, metavar="integer",
                help="Type of classification [default: Complete links] (1: K-menoids; 2: K-means; 3: Ward; 4: Complete links; 5: Single links; 6: UPGMA; 7: WPGMA; 8: WPGMC; 9: UPGMC)"),
    make_option(c("-adv", "--advanced"), type="logical", action="store_true", 
                help="Activate advanced mode (print more outputs)"),
    make_option(c("-q", "--quiet"), type="logical", action="store_true",
                help="Activate quiet mode"),
    make_option(c("-V", "--verbose"), type="logical", action="store_true",
                help="Activate verbose mode"),
    make_option(c("-T", "--text"), type="logical", action="store_true",
                help="DO NOT print values on graph"),
    make_option(c("-n", "--nb_clusters"), type="integer", metavar="integer",
                help="Fix the number of clusters"),
    make_option(c("-r", "--removeDoublets"), type="logical", action="store_true", 
                help="Discard line containing the same information on all columns from analysis"),
    make_option(c("-b", "--bootstrap"), type="integer", default=500, metavar="integer",
                help="Number of bootstraps for Gap statistic (advanced mode)"),
    make_option(c("-D", "--distance"), type="integer", default=1, metavar="integer",
                help="Type of distance [default: Euclidian] (1: Euclidian, 2: Manhattan, 3: Jaccard, 4: Sokal & Michener, 5 Sorensen (Dice), 6: Ochiai)"),
    make_option(c("-H", "--header"), type="logical", action="store_true",
                help="Consider first row as header of columns"),
    make_option(c("-s", "--separator"), type="character", metavar="character", default="\t",
                help="Type of separator (default: tabulation)"),
    make_option(c("-N", "--nb_axis"), type="integer", default=2, metavar="integer",
                help="Number of axis for pca (default: 2)")
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
  
  checkMinCluster = function (o, def=""){
    if (opt[[o]] < 2){
      print_help(a)
      stop(paste("--",o ," must be upper or equal to 2",def,".\n",sep=""), call.=FALSE)
    }
  }
  checkMinCluster("max_clusters"," [by default: 6]")
  if(!is.null(opt$nb_clusters)) checkMinCluster("nb_clusters")
  
  if ((opt$classif_type < 1) || (opt$classif_type > 9)){
    print_help(a)
    stop("--classif_type must be comprise between 1 and 6 [by default: 2].\n", call.=FALSE)
  }
  
  if ((opt$nb_axis < 2) || (opt$nb_axis > 4)){
    print_help(a)
    stop("--nb_axis must be comprise between 2 and 4 [by default: 2].\n", call.=FALSE)
  }
  
  if ((opt$distance < 1) || (opt$distance > 6)){
    print_help(a)
    stop("--distance must be comprise between 1 and 6 [by default: 1 for Euclidian].\n", call.=FALSE)
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
  
  checkMaxCluster("max_clusters"," [by default: 6]")
  if(!is.null(opt$nb_clusters)) checkMaxCluster("nb_clusters")
}


#avoid doublets in row names
#r: row names vector
renameRownameDoublets = function(names.row){
  j=1
  for (i in 2:length(names.row)){
    if (names.row[i] == names.row[i-1]){
      j = j+1
      names.row[i] = paste(names.row[i], ".", j, sep="")
    }else{
      j = 1
    }
  }
  return (names.row)
}

#rename row and avoid doublets errors
renameRowname = function(d){
  names.row = as.character(d[,1])
  d=d[,-1]
  names.row = renameRownameDoublets(names.row)
  tryCatch({
    substr(names.row, 1, 25) -> rownames(d)
    return(d)
  }, warning = function(w) {
    names.row = renameRownameDoublets(substr(names.row, 1, 25))
    names.row -> rownames(d)
    return(d)
  }, error = function(e) {
    return(d)
  })
}

# Discard row from a reaction dataset that have the same conditions in each columns
#x: dataframe
discardRowCondDoublets = function(x){
  row_doublets <- list()
  j = 0
  for (i in 1:nrow(x)){
    #uniq remove doublets in a vector, so return 1 only if there is only 1
    if( (length(unique(as.integer(x[i,])))==1)){
      #print(row.names(x[i,]))
      j = j +1
      row_doublets[[j]] = i
    }
  }
  if(length(row_doublets)!=0){
    removed_reacs = row.names(x[unlist(row_doublets),])
    removed_conds = x[unlist(row_doublets), 1]
    removed = cbind(removed_reacs, removed_conds)
    colnames(removed) = c("condition", "")
    assign("removed", removed,.GlobalEnv)
    writeTsv("removed", v=F)
  }
  return (x[-unlist(row_doublets),])
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

#df: dataframe
#d: distance type
getDistance = function(df, d){
  dists=c("euclidian", "manhattan", 1, 2, 5, 7)
  if(d < 3) dist(df, method = dists[d])
  else dist.binary(df, method = dists[d])
}

#Inputs: x : a matrix
#filename of the saved file
#Prints the matrix, save the matrix
writeTsv = function(x, cl=T, v=T){
  #print on stdout
  if (isTRUE(v)) cat(paste("\n", gsub("_", " ", toupper(x)), ":\n", sep=""))
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
  if (isTRUE(v)){
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

plotBestClustering = function(sub_title, values, values_type, optimal_nb_clusters, interval = 1, min_x=2, best=NULL, val2=NULL){
  plotAxis(1, 2, max_clusters)
  
  if (interval >= 1) axisSeq=round(values)
  else axisSeq = c(0, max(values) +0.1)
  
  #case of plotting gap statistics
  if (min_x < 2) best_y=values[optimal_nb_clusters]
  #case of fusion levels
  else if (!is.null(val2)) best_y=values[optimal_nb_clusters -1]
  else best_y = max(values)
  
  #for non-elbow plots
  if (!is.null(val2)) best = round(max(val2),2)
  else if (is.null(best)) best = round(max(values),2)
  
  plotAxis(2, min(axisSeq), max(axisSeq), interval)
  title(main="Optimal number of clusters", line=2, cex.main=2)
  mtext(text=sub_title, font=3, cex=1.2, line = 0.5)
  abline(v=optimal_nb_clusters, col="red", lty=2, lwd=2)
  points(optimal_nb_clusters, best_y, pch=19, col="red", cex=2)
  
  if (!is.null(val2)) t_values = val2
  else t_values = values
  
  if(isTRUE(text)) text(y=values, x=min_x:max_clusters, labels=round(t_values,2), cex=1.2, pos=4, col="red")
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
# d: data
# dis: distance matrix
#Ouput: Hierarchical classification
getCAH = function(t, d, dis){
  if(t>2){
    if (t==8 | t==9 ) checkEuclidean(dis)
    #cah: classification hierarchic ascending
    cah = hclust(dis, method=getClassifType(t))
  #automaticly ordering by clusters
  return (reorder.hclust(cah, dis))
  }
}

# Selects best algo based on cophenetic calculation
# df: data (or distance matrix for hierarchic)
selectBestCAH = function (d, dis, v=F){
  temp = 0
  for (i in 3:9){
    cah = getCAH(d, i)
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
  methods = c("kmedoids","kmeans","ward.D2", "complete", "single", "average", "mcquitty", "median", "centroid")
  methods[t]
}

#Agglomerative coefficient ()
getCoefAggl = function(c)
  coef.hclust(c)

#Inputs: 
# t: number of type of classification
# d: data (or distance for pam)
# k: number of clusterting
#Ouput: Non-hierarchical classification
getCNH = function(t, d, dis, k){
  if (t==1) return (pam(dis, k, diss=T))
  else if (t==2){
    checkEuclidean(dis)
    return (kmeans(d, centers=k, nstart=100))
  }
}

getClassif = function(t, n, d, dis){
  if(t>2) getCAH(t, d, dis)
  else {
    list_cnh = list("method"=getClassifType(t))
    for (k in 2:(n+1)) list_cnh[[k]] = getCNH(t, d, dis, k)
    return(list_cnh)
  }
}

checkEuclidean = function(dis){
  if(attributes(dis)$method != "euclidean") 
    stop("Distance should be euclidean with this classification method.", call.=FALSE)
}

# Inputs: 
# t: number of type of classification
# k: number of clusters
# c: hierarchical classification
# d: data
# Output: partitionning contening k clusters
getClusters = function(k, c) {
  if (c$method == "kmedoids") c[[k]]$clustering
  else if(c$method == "kmeans" ) c[[k]]$cluster
  else cutree(c, k)
}

getClusterPerPart = function (n, c){
  cl = list()
  for (k in 2:n){
    cl[[k-1]] = getClusters(k, c)
  }
  return (cl)
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
writeClusters = function(d, cl, r=FALSE, v=FALSE){
  nb_cl = length(levels(as.factor(cl)))
  clusters = matrix(NA, length(cl), nb_cl)
  for (i in 1:nb_cl ){
    if (r == FALSE){
      #get metabolites from clusters and put into a column of the output matrix
      # from the begining of the column to the end of the vector of metabolites names
      clusters[c(1:length(cl[cl==i])),i] =  names(cl[cl==i])
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
  writeTsv("clusters", F, v=v)
}

############################################################
#          Cophenetic (dendrogram distance matrix)
############################################################

# Distance matrix between each leaf of the dendogramm
#Inputs:
# d : data
# cah : hierarchical classification
plotCohenetic=function(dis, cah){
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
# n: maximum number of clusters
# d: dataframe
# cl: list of clusters per partition
getRelativeBetweenPerPart = function(n, d, cl){
  d=as.matrix(d)
  between = rep(0, n-1)
  # total sum of square
  TSS = sum(scale(d, scale = FALSE)^2)
  for (i in 2:n) {
    cl_k = as.factor(cl[[i-1]])
    # apply(d, 2, mean) : centroids for each column
    # as.vector(table(cl) : size of each clusters
    # t : vector rotation for arithmetic with other row or column vectors
    between[i-1] = sum(t((t(getClusterCentroids(d,cl_k))-apply(d, 2, mean))^2) * as.vector(table(cl_k)))/TSS
  }
  return (100*between)
}

# tapply(data[,i], Cla, mean) :
# centroids of each clusters for a column i
# sapply(1:ncol(data), function(i) tapply(data[,i], Cla, mean)) :
# centroids of each clusters for each column
getClusterCentroids = function(d, cl){
  sapply(1:ncol(d), function(i) tapply(d[,i], cl, mean))
}

# Difference between each case of a vector
getBetweenDifferences = function(between){
  # apply produce a list, unlist convert in vector
  diff = unlist(sapply(1:length(between), function(i) between[i]-between[i-1]))
  return (as.vector(cbind(between[1], t(diff))))
  #-n-1 to remove the last NA value (pairwise comparison)
  #between[1] to get the difference with 1 cluster
}

getWithin = function(d, cl, k) {
  nk = length(cl[cl==k]) #number of individuals in the cluster
  d1 = scalecenter(d)
  return( nk * sum(getClusterCentroids(d1,cl)[k,]^2) / nrow(d))
}

# cl: list of clusters per partition
getRelativeWithinPerCluster = function(cls, d) {
  n = length(cls)
  within = matrix(NA, n-1, n)
  rownames(within) = seq(2, n)
  colnames(within) = paste("G", seq(1, n), sep="")
  for (k in 2:n){
    cl = cls[[k-1]]
    for (i in 1:length(table(cl)) ){
      within[k-1, i] = getWithin(d, cl, i)
    }
    within[k-1,] = within[k-1,] / sum(as.numeric(na.omit(within[k-1,])))
  }
  return (within)
}

# Between inertia differences between a partionning and the previous
plotBetweenDiff = function(between_diff) {
  if (isTRUE(verbose)) cat("\nBETWEEN DIFFERENCES:\n")
  optimal_nb_clusters = which.max(between_diff)+1
  savePdf("between_differences.pdf")
  plot(2:(length(between_diff)+1), between_diff, type="b", ylim=c(round(min(between_diff))-1,round(max(between_diff))+1), xlim=c(2,(length(between_diff)+2)), xlab="Nb. of clusters", ylab="Between-cluster variation (%)", col="grey", axes=F)
  plotBestClustering("Largest between differences method", between_diff, " variation with the previous partitionning (%)", optimal_nb_clusters)
  suprLog = dev.off()
}

plotFusionLevels = function(n, c) {
  if (isTRUE(verbose)) cat("\nFUSION LEVELS:\n")
  fusion = rev(c$height)
  diff = unlist(sapply(1:n, function(i) fusion[i-1]-fusion[i]))
  fusion = fusion[1:(n-1)]
  optimal_nb_clusters = which.max(diff)+1
  savePdf("fusion_levels.pdf")
  plot(2:n, fusion, type="b", ylim=c(round(min(fusion))-1,round(max(fusion))+1), xlim=c(2,n+1), xlab="Nb. of clusters", ylab="Cophenetic distance", col="grey", axes=F)
  plotBestClustering("Fusion level method", fusion, " gain with the previous fusion level", optimal_nb_clusters, val2=diff)
  suprLog = dev.off()
}

#x: vector of between inertia for k partitions
plotElbow = function(x) {
  if (isTRUE(verbose)) cat("\nELBOW:\n")
  n = length(between) +1
  within = c(100, 100 - x)
  ratio = within[1:(n-1)] / within[2:n]
  best = round(min(ratio),2)
  optimal_nb_clusters = which.min(ratio)
  savePdf("elbow.pdf")
  plot(1:n, within, type="b", ylim=c(min(within)-1,101), xlim=c(1,n+1), xlab="Nb. of clusters", ylab="Relative within inertia (%)", col="grey", axes=F)
  plotBestClustering("Elbow method", within, " Wk/Wk+1 ratio ", optimal_nb_clusters, 5, 1, best)
  suprLog = dev.off()
}

################################
#          Silhouette
################################

#Ouput: an ordered silhouette object
getSilhouette = function(d, cl_k, dis){
  sil = sortSilhouette(silhouette(cl_k, dis))
  rownames(sil) = row.names(d)[attr(sil,"iOrd")]
  return (sil)
}

getSilhouettePerPart =function(d, cl, dis){
  list_sil = list()
  for (k in 2:length(cl)) {
    list_sil[[k-1]] = getSilhouette(d, cl[[k-1]], dis)
  }
  return(list_sil)
}

# sils: list of silhouettes objects per partition
getMeanSilhouettePerPart = function(sils){
  unlist(sapply(1:length(sil), function(i) summary(sil[[i]])$avg.width))
}

# Plot the best average silhouette width for all clustering possible
# mean_sils: vector of silhouette average width
plotSilhouettePerPart = function(mean_silhouette){
  if (isTRUE(verbose)) cat("\nSILHOUETTE:\n")
  savePdf("average_silhouettes.pdf")
  optimal_nb_clusters = which.max(mean_silhouette)+1
  plot(2:(length(sil)+1), mean_silhouette, type="b", xlim=c(2,length(sil)+2), ylim=c(0,max(mean_silhouette)+0.1), col="grey", xlab="Nb. of clusters", ylab="Average silhouette width", axes=F)
  plotBestClustering("Silhouette method", mean_silhouette,"n average width", optimal_nb_clusters, 0.1)
  suprLog = dev.off()
  return (optimal_nb_clusters)
}

#sil_k: a silhouette object
plotSilhouette = function(sil_k){
  pdf("silhouette.pdf")
  setGraphicBasic()
  par(mar=c(4, 12, 3, 2))
  plot(sil_k, max.strlen=25, main=" ", sub= "", do.clus.stat=TRUE, xlab="Silhouette width", cex.names=0.8, col=colorClusters(sil_k[,1]), nmax.lab=100, do.n.k = FALSE, axes=F)
  mtext(paste("Average silhouette width:", round(summary(sil_k)$avg.width,3)), font=2, cex=1.5, line=1)
  plotAxis(1, 0, 1, 0.2)
  sil_scores = cbind(row.names(sil_k), sil_k[,1], sil_k[,3])
  #colnames(sil_scores) = c("Chemicals", "Cluster", "Silhouette score")
  assign("sil_scores", sil_scores, .GlobalEnv)
  writeTsv("sil_scores", v=F)
  suprLog = dev.off()
}

###################################
#          GAP STATISTICS
###################################

#B: nb of bootstrap
getGapPerPart = function(n, d, c, B=500){
  #FUN mus have only two args in this order and return a list with an object cluster
  gapFun = function(x, k) list(cluster = getClusters(k, c))
  clusGap(d, FUN=gapFun, K.max=n, verbose=F, B=B)
}

#g: gap object
getGapBest = function (g, M="Tibs2001SEmax"){
  with(g, maxSE(Tab[,"gap"], Tab[,"SE.sim"], method=M))
}

# Plot the gap statistics width for all clustering possible
#TODO: HERE
plotGapPerPart = function(n, d, c, B=500, v=T){
  if (isTRUE(v)) cat("\nGAP STATISTICS:\n")
  if(classif$method=="kmeans" & (n > 10 | nrow(d)>=100)) plural=c("few ", "s")
  else plural=c("","")
  if(classif$method=="kmeans" | nrow(d)>=100 ) cat(paste("It could take a ",plural[1], "minute",plural[2],"...\n",sep=""))
  gap = getGapPerPart(n, d, c, B)
  savePdf("gap_statistics.pdf")
  optimal_nb_clusters = getGapBest(gap)
  gap_k=round(gap$Tab,3)
  best = paste(gap_k[,"gap"][optimal_nb_clusters], ">",gap_k[,"gap"][optimal_nb_clusters+1],"-",gap_k[,"SE.sim"][optimal_nb_clusters])
  plot(gap, arrowArgs = list(col="gray", length=1/15, lwd=2, angle=90, code=3), type="b", xlim=c(1,n+1), ylim=c(0,max(gap$Tab[,"gap"])+0.1), col="grey", xlab="Nb. of clusters", ylab=expression(Gap[k]), main="",axes=F)
  plotBestClustering("Gap statistics method", gap$Tab[,"gap"]," gap value", optimal_nb_clusters, 0.1, 1, best)
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

printSummary = function(between, diff, sil, adv, gap=NULL){ 
  #TODO: no n = nrow(data)
  summary = cbind(between, diff, 100-between, sil)
  names = c("Between-inertia (%)", "Between-differences (%)", "Within-inertia (%)", "Silhouette index") 
  if(isTRUE(adv)) {
    summary = cbind(summary, gap$Tab[,"gap"][-1], gap$Tab[,"SE.sim"][-1])
    names=c(names, "Gap", "Gap SE")
  }
  rownames(summary) = seq(2, (length(between)+1))
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
# df: data frame
# d: a distance object
# s: an organised silhouette object
# c: CAH
# c: clusters from CAH
heatMap = function(df, d, s=NULL, c=NULL, cl=NULL, text=FALSE){
  
  if(!is.null(s)){
    order = attr(s,"iOrd")
    cl_sizes = summary(s)$clus.size
    title = "silhouette\'s scores"
    colors = colPers(length(levels(as.factor(s[,1]))))
  }else{
    order = c$order
    cl_sizes = getOrderedClusterSize(cl[order])
    title="dendrogram"
    colors = orderColors(c, cl)
  }

  matrix=as.matrix(d)
  matrix=matrix[order, order]
  rownames(matrix) <- rownames(df)[order] -> labels
  #if(tri == TRUE) matrix[!lower.tri(matrix)] = NA
  #image(1:ncol(matrix), 1:ncol(matrix), t(matrix), axes=F, xlab="", ylab="")

  options(warn = -1)
  if(nrow(df) > 200 ) png("heat_map.png", 2000, 2000)
  else pdf("heat_map.pdf")
  
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
plotDendrogram = function(t, k, c, d, n, cl){

  pdf("dendrogram.pdf")
  setGraphicBasic()
  par(mar=c(2,5,5,1))
  plot(c, hang=-1, ylim=c(0,max(c$height)), xlim=c(0,length(c$labels)), sub="", cex=0.8, font=3, ylab="Cophenetic distance", main="Dendrogram", axes=F)
  plotAxis(2, 0, max(c$height))
  abline(h=rev(c$height)[1:n], col="gray", lty=2, lwd=1)
  #projection of the clusters
  rect.hclust(c, k=as.numeric(k), border=orderColors(c, cl))
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

#nf: number of factorial axis
plotPca = function(t, k, cl, d, nf=2){
  pca = dudi.pca(d, scannf=F, nf=nf)
  if(nrow(d) > 200 ) {
    png("pca.png", 1000, 1000); cex = 1
  }else{
    pdf("pca.pdf"); cex = 0.6
  }
  par(mar=c(0,0,4.1,0))
  title = paste("Cumulated inertia:", round((pca$eig[nf-1]+pca$eig[nf])/sum(pca$eig),4)*100, "%")
  s.class(addaxes=F, cbind(pca$li[,nf-1] , pca$li[,nf]), ylim=c(min(pca$li[,nf])-3, max(pca$li[,nf])+3), xlim=c(min(pca$li[,nf-1])-3, max(pca$li[,nf-1])+3), csub=1.5, as.factor(cl), grid=F, col=colPers(k))
  mtext(title, font=2, cex=1.5, line=1)
  abline(h=0, v=0, lty=2, lwd=2, col="grey")
  text(x=pca$li[,nf-1], y=pca$li[,nf], labels=rownames(pca$li), col=colorClusters(cl), cex=cex)
  pca_coord = cbind(rownames(pca$li), pca$li[,1],pca$li[,2])
  #colnames(pca_coord) = c("Chemicals", "Axis 1", "Axis 2")
  assign("pca_coord", pca_coord, .GlobalEnv)
  writeTsv("pca_coord", v=F)
  par(fig=c(0.8,1,0.82,1),new=TRUE)
  if(isTRUE(advanced)) plotInertiaPca(pca)
  suprLog = dev.off()
}

# nf: number of inertia bar plot corresponding to factorial axis
plotInertiaPca = function (pca, nf=4){
  inertia = round(pca$eig/sum(pca$eig)*100, 1)
  par(mar=c(2, 0, 1, 1) + 0.1)
  plot(inertia, type="h", lwd=10, lend=1, xlim=c(0,nf+0.2),ylim=c(0,max(inertia+7)),col="grey75",font=2, axes=F, xlab="", ylab="")
  title(sub="Inertia (in %)", line=0, cex.sub=0.7, font.sub=3)
  text(1:nf,inertia[1:nf]+5, inertia[1:nf], cex = 0.6)
  par(new=TRUE); par(mar=c(0,0,0,0)) ; plot(0:1,0:1, axes=F, type="n")
  rect(0,0.1,0.9,0.9, border="grey65")
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
getCtrVar = function(t, k, cl, d) {
  
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
getPdis = function(t, k, cl, d) {
  
  #for each metabolite contribution (in column), sum the k clusters values
  return(apply(getCtrVar(t, k, cl, d), 2, sum))
}

# Inputs: 
# t: number of type of classification
# n: number max of clusters
# cls: list of clusters
# d: data
# index: pdis or rho2 calculation
getPdisPerPartition = function(t, n, cls, d){
  
  pdis_per_partition = matrix(NA, n-1, ncol(d))
  rownames(pdis_per_partition) = seq(2, n)
  colnames(pdis_per_partition) = colnames(d)
  
  for (k in 2:n){
    res = getPdis(t, k, cls[[k-1]], d)
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
librairies = c("cluster", "optparse", "gclus", "ade4", "scales", "beepr")
for (l in librairies){
  if (! (l %in% installed.packages()[,"Package"])) install.packages(l, repos = "http://cran.us.r-project.org", quiet = T)
  library(l, character.only = TRUE)
}

#Get arguments
args = getArgs()
checkArg(args)
opt = parse_args(args)

#Global variables settings
nb_clusters = opt$nb_clusters
max_clusters = opt$max_clusters
classif_type = opt$classif_type
bootstrap = opt$bootstrap
advanced = "advanced" %in% names(opt)
verbose = ( !("quiet" %in% names(opt)) | ("verbose" %in% names(opt)))
verboseNiv2 = ("verbose" %in% names(opt))
remove_doublets = ("removeDoublets" %in% names(opt))
text = !("text" %in% names(opt))
header = ("header" %in% names(opt))
if (!is.null(opt$workdir)) setwd(opt$workdir)

#Loading data "R_t03.2.csv"
data = read.table(opt$infile, header=header, sep=opt$separator, dec=".")
postChecking(args, data)
#rename row and avoid doublets errors
data = renameRowname(data)
#if(isTRUE(remove_doublets)) data = discardRowCondDoublets(data)

printProgress = function (v, val){
  if(isTRUE(v)) 
    cat(paste("\n[", format(Sys.time(), "%X"), "] ", val ," in progress...\n"), sep="")
}

#Perform classification
printProgress(verboseNiv2, "Distance calculation")
dis = getDistance(data, opt$distance)
printProgress(verboseNiv2, "Classification")
classif = getClassif(classif_type, max_clusters, data, dis)
printProgress(verboseNiv2, "Clustering")
list_clus = getClusterPerPart(max_clusters+1, classif)

#Indexes
if(classif_type>2){
  printProgress(verboseNiv2, "Cophenetic calculation")
  plotCohenetic(dis, classif)
  if(isTRUE(advanced) & isTRUE(verbose)) cat(paste("\nAGGLOMERATIVE COEFFICIENT: ", round(getCoefAggl(classif),3), "\n", sep=""))
  plotFusionLevels(max_clusters, classif)
}

#Inertia
between = getRelativeBetweenPerPart(max_clusters, data, list_clus)
diff = getBetweenDifferences(between)
plotElbow(between)

#Silhouette analysis
printProgress(verboseNiv2, "Silhouette calculation")
sil = getSilhouettePerPart(data, list_clus, dis)
mean_silhouette = getMeanSilhouettePerPart(sil)
optimal_nb_clusters = plotSilhouettePerPart(mean_silhouette)
if(!is.null(nb_clusters)) optimal_nb_clusters = nb_clusters
sil_k = sil[[optimal_nb_clusters-1]]
plotSilhouette(sil_k)

#cl_temp, because writeTsv(clusters) recreate a different object named clusters
clusters = list_clus[[optimal_nb_clusters-1]]
cl_temp = clusters
gap = NULL

#Advanced indexes
if (isTRUE(advanced)){

  gap = plotGapPerPart(max_clusters, data, classif, bootstrap)
  plotGapPerPart2(gap, max_clusters)
  
  contribution = 100 * getCtrVar(classif_type, optimal_nb_clusters, clusters, data)
  discriminant_power = 100 * getPdisPerPartition(classif_type, max_clusters, list_clus, data)
  within_k = getRelativeWithinPerCluster(list_clus, data)
  
  for (i in c("contribution", "discriminant_power", "within_k"))
    writeTsv(i, v=F)
}

#Plots
if(classif_type > 2) plotDendrogram(classif_type, optimal_nb_clusters, classif, data, max_clusters, clusters)
printProgress(verboseNiv2, "PCA")
plotPca(classif_type, optimal_nb_clusters, clusters, data, opt$nb_axis)
printProgress(verboseNiv2, "Heatmap calculation")
if(classif_type <= 2 || isTRUE(advanced)){
  heatMap(data, dis, sil_k, text=(nrow(data) < 100))
}else{
  heatMap(data, dis, c=classif, cl=clusters, text=(nrow(data) < 100))
}

#Final outputs
summary = printSummary(between, diff, mean_silhouette, advanced, gap)
writeTsv("summary", v=verbose)
writeClusters(data, clusters, TRUE, v=( (verbose) & (nrow(data) < 100) ) )
if (!isTRUE(verbose)) cat(paste("Optimal number of clusters:", optimal_nb_clusters,"\n"))
if (isTRUE(verboseNiv2)) beep("ping")

#errors
if (optimal_nb_clusters==max_clusters) message("\n[WARNING] The optimal number of clusters equals the maximum number of clusters. \nNo cluster structure has been found.")
if(min(table(cl_temp))==1) message("\n[WARNING] A cluster with an only singleton biased the silhouette score.")