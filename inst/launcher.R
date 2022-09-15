# Global variables settings
TEXT <- TRUE # print values on graph (for optimum partition and heatmap)
NB_ROW_MAX <- 200 # max row to have pdf, otherwise, some plots are in png
DIM_PNG <- 2000
VERBOSE_NIV2 <- FALSE
VERBOSE <- FALSE
MAX_CHAR_LEN <- 25 # maximum length of individual s names

library("optparse")
library("autoCluster")
source("R/fingerprint_clustering.R")
source("R/plot.R")

#' Load libraries
#'
#' Load libraries
#'
#' @param librairies A list of libraries
#'
#' @export
load_libraries <- function(librairies) {
  for (l in librairies) {
    suppressPackageStartupMessages(
      library(
        l,
        character.only = TRUE,
        warn.conflicts = FALSE,
        quietly = TRUE,
        verbose = FALSE
      )
    )
  }
}
load_libraries(pkgload::pkg_desc()$get_deps()[, 2])

getArgs <- function() {
    option_list <- list(
        make_option(
            c("-i", "--infile"),
            type = "character",
            metavar = "character",
            help = "Fingerprint file name"
        ),
        make_option(
            c("--output1"),
            type = "character",
            default = "average_silhouette.pdf",
            metavar = "character",
            help = "Average silhouettes file name [default: %default]"
        ),
        make_option(
            c("--output2"),
            type = "character",
            default = "silhouette.pdf",
            metavar = "character",
            help = "Silhouette file name [default: %default]"
        ),
        make_option(
            c("--output3"),
            type = "character",
            default = "pca.pdf",
            metavar = "character",
            help = "PCA file name [default: %default]"
        ),
        make_option(
            c("--output4"),
            type = "character",
            default = "heatmap.pdf",
            metavar = "character",
            help = "Heatmap file name [default: %default]"
        ),
        make_option(
            c("--output5"),
            type = "character",
            default = "summary.tsv",
            metavar = "character",
            help = "Summary file name [default: %default]"
        ),
        make_option(
            c("--output6"),
            type = "character",
            default = "clusters.tsv",
            metavar = "character",
            help = "Clusters file name [default: %default]"
        ),
        make_option(
            c("--output7"),
            type = "character",
            default = "dendrogram.tsv",
            metavar = "character",
            help = "Dendrogram file name [default: %default]"
        ),
        make_option(
            c("--output8"),
            type = "character",
            default = "shepard_graph.pdf",
            metavar = "character",
            help = "Shepard graph file name [default: %default]"
        ),
        make_option(
            c("--output9"),
            type = "character",
            default = "fusion_levels.pdf",
            metavar = "character",
            help = "Fusion levels file name [default: %default]"
        ),
        make_option(
            c("-m", "--maxClusters"),
            type = "integer",
            default = 6,
            metavar = "integer",
            help = "Maximum number of clusters [default: %default]"
        ),
        make_option(
            c("-t", "--classifType"),
            type = "integer",
            default = 3,
            metavar = "integer",
            help = "Type of classification [default: Complete links] (1: K-menoids; 2: K-means; 3: Ward; 4: Complete links; 5: Single links; 6: UPGMA; 7: WPGMA; 8: WPGMC; 9: UPGMC)"
        ),
        make_option(
            c("-a", "--advanced"),
            type = "logical",
            action = "store_true",
            help = "Activate advanced mode (print more outputs)"
        ),
        make_option(
            c("-q", "--quiet"),
            type = "logical",
            action = "store_true",
            help = "Activate quiet mode"
        ),
        make_option(
            c("-v", "--verbose"),
            type = "logical",
            action = "store_true",
            help = "Activate \"super-verbose\" mode"
        ),
        make_option(
            c("-n", "--nbClusters"),
            type = "integer",
            metavar = "integer",
            help = "Fix the number of clusters"
        ),
        make_option(
            c("-d", "--distance"),
            type = "integer",
            default = 1,
            metavar = "integer",
            help = "Type of distance [default: Euclidian] (1: Euclidian, 2: Manhattan, 3: Jaccard, 4: Sokal & Michener, 5 Sorensen (Dice), 6: Ochiai)"
        ),
        make_option(
            c("-H", "--header"),
            type = "logical",
            action = "store_false",
            help = "Consider first row as header of columns"
        ),
        make_option(
            c("-s", "--separator"),
            type = "integer",
            metavar = "integer",
            default = 1,
            help = "Type of separator [default: tabulation] (1: Tabulation, 2: Semicolon"
        ),
        make_option(
            c("-N", "--nbAxis"),
            type = "integer",
            default = 2,
            metavar = "integer",
            help = "Number of axis for pca (default: 2)"
        )
    )
    # if -h, avoid exit with error
    args <- commandArgs(trailingOnly = TRUE)
    if ("-h" %in% args) {
        q("no", status = 0, runLast = FALSE)
    }
    return(OptionParser(option_list = option_list))
}

# Check the arguments validity
# Inputs:
# a: arguments (optionParser object)
checkArg <- function(a) {
    opt <- parse_args(a)
    # o: one argument from the list of arguments
    # def: defaul message

    checkMinCluster <- function(o, n, def = "") {
        if (opt[[o]] < n) {
            stop(paste("--", o, " must be upper or equal to ", n, " ", def, ".\n", sep = ""), call. = FALSE)
        }
    }

    checkMinCluster("maxClusters", 3, " [by default: 6]")
    if (!is.null(opt$nbClusters)) {
        checkMinCluster("nbClusters", 2)
    }

    if ((opt$classifType < 1) || (opt$classifType > 9)) {
        stop("--classifType must be comprise between 1 and 9 [by default: 2].\n", call. = FALSE)
    }

    if ((opt$nbAxis < 2) || (opt$nbAxis > 4)) {
        stop("--nbAxis must be comprise between 2 and 4 [by default: 2].\n", call. = FALSE)
    }

    if ((opt$distance < 1) || (opt$distance > 6)) {
        stop("--distance must be comprise between 1 and 6 [by default: 1 for Euclidian].\n", call. = FALSE)
    }


    checkFile <- function(o) {
        if (!file.exists(opt[[o]])) {
            stop(paste("--", o, " name does not exist\n", sep = ""), call. = FALSE)
        }
    }
    # if(!is.null(opt$workdir)) checkFile("workdir")
    if (!is.null(opt$infile)) {
        checkFile("infile")
    } else {
        stop(paste("--infile is required\n", sep = ""), call. = FALSE)
    }

    return(opt)
}

# Checking clusters args after data loading
# Inputs:
# a: arguments (optionParser object)
# d: data
# o: one argument from the list of arguments
# def: defaul message
postChecking <- function(a, d) {
    opt <- parse_args(a)

    checkMaxCluster <- function(o, def = "") {
        if (opt[[o]] > nrow(d)) {
            # print_help(a)
            stop(
                paste(
                    "--",
                    o,
                    " must be lower or equal to the fingerprint",
                    def,
                    ".\n",
                    sep = ""
                ),
                call. = FALSE
            )
        }
    }

    checkMaxCluster("maxClusters", " [by default: 6]")
    if (!is.null(opt$nbClusters)) {
        checkMaxCluster("nbClusters")
    }
}


# Get arguments
args <- getArgs()
tryCatch(
    {
        opt <- checkArg(args)
    },
    error = function(e) {
        # print_help(args)
        stop(e[[1]], call. = FALSE)
    }
)


# Global variables settings
NB_CLUSTERS <- opt$nbClusters
MAX_CLUSTERS <- opt$maxClusters
CLASSIF_TYPE <- opt$classifType
ADVANCED <- "advanced" %in% names(opt)
NB_AXIS <- opt$nbAxis
# Discard line containing the same information on all columns from analysis
REMOVE_DOUBLETS <- TRUE
HEAD <- !("header" %in% names(opt))
# if (!is.null(opt$workdir)) setwd(opt$workdir)
if (isTRUE(VERBOSE_NIV2)) {
    start_time <- Sys.time()
}

if (opt$separator == 1) {
    SEP <- "\t"
} else if (opt$separator == 2) {
    SEP <- ";"
} else {
    # print_help(args)
    stop(paste("--separator must be 1: Tabulation or 2: Semicolon."), call. = FALSE)
}

# Loading data
data <- read.table(
    opt$infile,
    header = HEAD,
    sep = SEP,
    dec = ".",
    row.names = 1
)
data <- preProcessData(data, HEAD)
# data=data[,1:10]


# Perform classification
printProgress(VERBOSE_NIV2, "Distance calculation")
dis <- getDistance(data, opt$distance)
if (CLASSIF_TYPE < 3) {
    printProgress(VERBOSE_NIV2, "Classification")
}
classif <- getClassif(CLASSIF_TYPE, MAX_CLUSTERS, data, dis)
list_clus <- getClusterPerPart(MAX_CLUSTERS + 1, classif)
optimal_nb_clusters <- 2
clusters <- list_clus[[optimal_nb_clusters - 1]]

NB_BIOMARK <- 20
discr <- getDiscriminantVariables(
    CLASSIF_TYPE,
    optimal_nb_clusters,
    clusters,
    data,
    NB_BIOMARK
)
plotDiscriminantVariables(discr)

discriminant_power <- getPdisPerPartition(
    CLASSIF_TYPE,
    MAX_CLUSTERS,
    list_clus,
    data
) * 100

# Silhouette analysis
sil <- getSilhouettePerPart(data, list_clus, dis)
mean_silhouette <- getMeanSilhouettePerPart(sil)
plotSilhouettePerPart(mean_silhouette, sil, VERBOSE)
