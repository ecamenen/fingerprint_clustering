getArgs = function(){
  option_list = list(
    make_option(c("-i", "--infile"), type="character", metavar="character",
                help="Fingerprint file name"),
    make_option(c( "--output1"), type="character", default="average_silhouette", 
                metavar="character", help="Average silhouettes file name [default: %default]"),
    make_option(c( "--output2"), type="character", default="silhouette", 
                metavar="character", help="Silhouette file name [default: %default]"),
    make_option(c( "--output3"), type="character", default="dendrogram", 
                metavar="character", help="Dendrogram file name [default: %default]"),
    make_option(c( "--output4"), type="character", default="pca", 
                metavar="character", help="PCA file name [default: %default]"),
    make_option(c( "--output5"), type="character", default="heatmap", 
                metavar="character", help="Heatmap file name [default: %default]"),
    make_option(c( "--output6"), type="character", default="summary", 
                metavar="character", help="Summary file name [default: %default]"),
    make_option(c( "--output7"), type="character", default="clusters", 
                metavar="character", help="Clusters file name [default: %default]"),
    make_option(c( "--output8"), type="character", default="shepard_graph", 
                metavar="character", help="Shepard graph file name [default: %default]"),
    make_option(c( "--output9"), type="character", default="fusion_levels", 
                metavar="character", help="Fusion levels file name [default: %default]"),
    make_option(c("-m", "--maxClusters"), type="integer", default=6, metavar="integer",
                help="Maximum number of clusters [default: %default]"),
    make_option(c("-t", "--classifType"), type="integer", default=4, metavar="integer",
                help="Type of classification [default: Complete links] (1: K-menoids; 2: K-means; 3: Ward; 4: Complete links; 5: Single links; 6: UPGMA; 7: WPGMA; 8: WPGMC; 9: UPGMC)"),
    make_option(c("-a", "--advanced"), type="logical", action="store_true", 
                help="Activate advanced mode (print more outputs)"),
    make_option(c("-q", "--quiet"), type="logical", action="store_true",
                help="Activate quiet mode"),
    make_option(c("-v", "--verbose"), type="logical", action="store_true",
                help="Activate \"super-verbose\" mode"),
    make_option(c("-n", "--nbClusters"), type="integer", metavar="integer",
                help="Fix the number of clusters"),
    make_option(c("-d", "--distance"), type="integer", default=1, metavar="integer",
                help="Type of distance [default: Euclidian] (1: Euclidian, 2: Manhattan, 3: Jaccard, 4: Sokal & Michener, 5 Sorensen (Dice), 6: Ochiai)"),
    make_option(c("-H", "--header"), type="logical", action="store_true",
                help="Consider first row as header of columns"),
    make_option(c("-s", "--separator"), type="integer", metavar="integer", default=1,
                help="Type of separator [default: tabulation] (1: Tabulation, 2: Semicolon"),
    make_option(c("-N", "--nbAxis"), type="integer", default=2, metavar="integer",
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
  
  checkMinCluster = function (o, def=""){
    if (opt[[o]] < 3){
      stop(paste("--",o ," must be upper or equal to 3",def,".\n",sep=""), call.=FALSE)
    }
  }
  checkMinCluster("maxClusters"," [by default: 6]")
  if(!is.null(opt$nbClusters)) checkMinCluster("nbClusters")
  
  if ((opt$classifType < 1) || (opt$classifType > 9)){
    stop("--classifType must be comprise between 1 and 9 [by default: 2].\n", call.=FALSE)
  }
  
  if ((opt$nbAxis < 2) || (opt$nbAxis > 4)){
    stop("--nbAxis must be comprise between 2 and 4 [by default: 2].\n", call.=FALSE)
  }
  
  if ((opt$distance < 1) || (opt$distance > 6)){
    stop("--distance must be comprise between 1 and 6 [by default: 1 for Euclidian].\n", call.=FALSE)
  }
  
  
  checkFile = function (o){
    if(!file.exists(opt[[o]])){
      stop(paste("--", o, " name does not exist\n", sep=""), call.=FALSE)
    }
  }
  #if(!is.null(opt$workdir)) checkFile("workdir")
  if(!is.null(opt$infile)) checkFile("infile")
  else stop(paste("--infile is required\n", sep=""), call.=FALSE)
  
  return (opt)
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
  
  checkMaxCluster("maxClusters"," [by default: 6]")
  if(!is.null(opt$nbClusters)) checkMaxCluster("nbClusters")
}

printProgress = function (v, val){
  if(isTRUE(v)) 
    cat(paste("\n[", format(Sys.time(), "%X"), "] ", val ,"in progress...\n"), sep="")
}

getTimeElapsed = function(start_time){
  time = as.numeric(as.difftime(Sys.time()-start_time), units="secs")
  secs = time %% 60
  time = (time - secs) /60
  mins = time %% 60
  hours = time / 60
  time = paste(mins, "min ", round(secs), "s\n", sep="")
  if (hours >= 1)  time = paste(round(hours), "h ", time, sep="")
  cat(paste("\nTime to run the process : ", time, sep=""))
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
source("fingerprint_clustering.R")

#Get arguments
args = getArgs()
tryCatch({
  opt = checkArg(args)
}, error = function(e) {
  print_help(args)
  stop(e[[1]], call.=FALSE)
})

#Global variables settings
NB_CLUSTERS = opt$nbClusters
MAX_CLUSTERS = opt$maxClusters
CLASSIF_TYPE = opt$classifType
NB_BOOTSTRAP = 500 #should be comprise between 100 and 1000
ADVANCED = "advanced" %in% names(opt)
VERBOSE = ( !("quiet" %in% names(opt)) | ("verbose" %in% names(opt)))
VERBOSE_NIV2 = ("verbose" %in% names(opt))
NB_AXIS =opt$nbAxis
REMOVE_DOUBLETS = T #Discard line containing the same information on all columns from analysis
TEXT = T #print values on graph (for optimum partition and heatmap)
HEAD = ("header" %in% names(opt))
#if (!is.null(opt$workdir)) setwd(opt$workdir)
if (isTRUE(VERBOSE_NIV2)) start_time = Sys.time()
NB_ROW_MAX = 200 #max row to have pdf, otherwise, some plots are in png
DIM_PNG = 2000

if (opt$separator==1){ SEP="\t"
}else if (opt$separator==2){ SEP=";"
}else{
  print_help(args)
  stop(paste("--separator must be 1: Tabulation or 2: Semicolon."), call.=FALSE)
}

#Loading data
data = read.table(opt$infile, header=HEAD, sep=SEP, dec=".")
if(ncol(data)==1){
  print_help(args)
  stop(paste("Check for the --separator (by default, tabulation)."), call.=FALSE)
}
postChecking(args, data)
#rename row and avoid doublets errors
data = renameRowname(data)
#remove columns containing characters
data = data[, unlist(sapply(1:ncol(data), function(i) is.numeric(data[,i])))]
if(isTRUE(REMOVE_DOUBLETS)){
  printProgress(VERBOSE_NIV2, "Loading data")
  data = discardRowCondDoublets(data)
}
if ( (nrow(data) > 3500) & (CLASSIF_TYPE > 2) ) stop("With more than 3000 rows to analyse, --classType must be 1: K-medoids or 2: K-means", call.=FALSE)
if ( isSymmetric(as.matrix(data)) & !HEAD) colnames(data) = rownames(data)

#Perform classification
printProgress(VERBOSE_NIV2, "Distance calculation")
dis = getDistance(data, opt$distance)
if(CLASSIF_TYPE < 3) printProgress(VERBOSE_NIV2, "Classification")
classif = getClassif(CLASSIF_TYPE, MAX_CLUSTERS, data, dis)
list_clus = getClusterPerPart(MAX_CLUSTERS+1, classif)

#Indexes
if(CLASSIF_TYPE > 2){
  printProgress(VERBOSE_NIV2, "Cophenetic calculation")
  plotCohenetic(dis, classif)
  if(isTRUE(ADVANCED) & isTRUE(VERBOSE)) cat(paste("\nAGGLOMERATIVE COEFFICIENT: ", round(getCoefAggl(classif),3), "\n", sep=""))
  plotFusionLevels(MAX_CLUSTERS, classif)
}

#Inertia
printProgress(VERBOSE_NIV2, "Index calculation")
between = getRelativeBetweenPerPart(MAX_CLUSTERS, data, list_clus)
diff = getBetweenDifferences(between)

#Silhouette analysis
sil = getSilhouettePerPart(data, list_clus, dis)
mean_silhouette = getMeanSilhouettePerPart(sil)
optimal_nb_clusters = plotSilhouettePerPart(mean_silhouette)
if(!is.null(NB_CLUSTERS)) optimal_nb_clusters = NB_CLUSTERS
sil_k = sil[[optimal_nb_clusters-1]]
plotSilhouette(sil_k)

#cl_temp, because writeTsv(clusters) recreate a different object named clusters
clusters = list_clus[[optimal_nb_clusters-1]]
cl_temp = clusters
gap = NULL

#ADVANCED indexes
if (isTRUE(ADVANCED)){
  
  if (nrow(data) < 100){
    printProgress(VERBOSE_NIV2, "Gap statistics calculation")
    gap = plotGapPerPart(MAX_CLUSTERS, data, classif, NB_BOOTSTRAP, v=VERBOSE)
    plotGapPerPart2(gap, MAX_CLUSTERS)
  }
  
  plotElbow(between)
  #decomment to have contribution per variable to the inertia of each clusters and to each partionning
  #contribution = 100 * getCtrVar(CLASSIF_TYPE, optimal_nb_clusters, clusters, data)
  #discriminant_power = 100 * getPdisPerPartition(CLASSIF_TYPE, MAX_CLUSTERS, list_clus, data)
  within_k = getRelativeWithinPerCluster(list_clus, data)
  
  # for (i in c("contribution", "discriminant_power"))
  #   writeTsv(i, v=F)
  writeTsv("within_k", v=F)
}

#dendrogram
if(CLASSIF_TYPE > 2) plotDendrogram(CLASSIF_TYPE, optimal_nb_clusters, classif, data, MAX_CLUSTERS, clusters)

#pca
printProgress(VERBOSE_NIV2, "PCA")
pca = dudi.pca(data, scannf=F, nf=NB_AXIS)
for (i in 1:NB_AXIS)
  for (j in i:NB_AXIS)
    if(i != j) plotPca(pca, data, clusters, i, j)

#Heatmap
printProgress(VERBOSE_NIV2, "Heatmap calculation")
textHeatmap = isTRUE(  isTRUE(TEXT) & (nrow(data) < 100) )
if(CLASSIF_TYPE <= 2 || isTRUE(ADVANCED)){
  heatMap(data, dis, sil_k, text=textHeatmap)
}else{
  heatMap(data, dis, c=classif, cl=clusters, text=textHeatmap)
}

#Final outputs
summary = printSummary(between, diff, mean_silhouette, ADVANCED, gap)
writeTsv("summary", v=VERBOSE)
#decomment to have the mean of the variables for each clusters
#getClusterCentroids(data, clusters)
#decomment to have the mean of every variables (only if each column is a different condition of the same variable)
#apply(getClusterCentroids(data, clusters), 1, mean)
if (nrow(data) > 100){
  cat("\nCLUSTER SIZES:")
  # t<- do not print "clusters"
  table(t <- clusters)
} 
writeClusters(data, clusters, TRUE, v=( (VERBOSE) & (nrow(data) < 100) ) )
if (!isTRUE(VERBOSE)) cat(paste("Optimal number of clusters:", optimal_nb_clusters,"\n"))
if (isTRUE(VERBOSE_NIV2)) beep("ping")
if (isTRUE(VERBOSE_NIV2)) getTimeElapsed(start_time)

#errors
if (optimal_nb_clusters==MAX_CLUSTERS) message("\n[WARNING] The optimal number of clusters equals the maximum number of clusters. \nNo cluster structure has been found.")
if(min(table(cl_temp))==1) message("\n[WARNING] A cluster with an only singleton biased the silhouette score.")