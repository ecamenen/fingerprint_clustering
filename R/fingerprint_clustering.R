
################################
#            MAIN
################################

# Pseudo-random settings:
# milisec * PID
set.seed(as.numeric(format(Sys.time(), "%OS2")) * 100 * Sys.getpid())

# Global variables settings
TEXT <- TRUE # print values on graph (for optimum partition and heatmap)
NB_ROW_MAX <- 200 # max row to have pdf, otherwise, some plots are in png
MAX_CHAR_LEN <- 25 # maximum length of individual s names

# Loading data
# dec=".")
# if (isTRUE(VERBOSE_NIV2)) start_time = Sys.time()


# decomment to have the mean of the variables for each clusters
# getClusterCentroids(data, clusters)
# decomment to have the mean of every variables (only if each column is a different condition of the same variable)
# apply(getClusterCentroids(data, clusters), 1, mean)
# if (nrow(data) > 100){
#   cat("\nCLUSTER SIZES:")
#   # t<- do not print "clusters"
#   table(t <- clusters)
# }

# if (!isTRUE(VERBOSE)) cat(paste("Optimal number of clusters:", optimal_nb_clusters,"\n"))
# if (isTRUE(VERBOSE_NIV2)) getTimeElapsed(start_time)

printProgress <- function(v, val) {
    if (isTRUE(v)) {
        cat(paste("\n[", format(Sys.time(), "%X"), "] ", val, "in progress..."), sep = "")
    }
}

getTimeElapsed <- function(start_time) {
    time <- as.numeric(as.difftime(Sys.time() - start_time), units = "secs")
    secs <- time %% 60
    time <- (time - secs) / 60
    mins <- time %% 60
    hours <- time / 60
    time <- paste(mins, "min ", round(secs), "s\n", sep = "")

    if (hours >= 1) {
        time <- paste(round(hours), "h ", time, sep = "")
    }

    cat(paste("\nTime to run the process : ", time, sep = ""))
}

################################
#          Parsing
################################


# rename row and avoid doublets errors
preProcessData <- function(d, header = FALSE, verbose = FALSE) {
    if (ncol(d) == 1) {
        stop(paste("Check for the separator (by default, tabulation)."),
            call. = FALSE
        )
    }

    # Discard line containing the same information on all columns from analysis
    REMOVE_DOUBLETS <- (nrow(d) > NB_ROW_MAX)
    d <- renameRowname(d)
    # remove columns containing characters
    # if(nrow(d) > NB_ROW_MAX) verbose = TRUE
    # get only columns with numeric values
    d <- d[, unlist(sapply(1:ncol(d), function(i) is.numeric(d[, i])))]

    if (isTRUE(REMOVE_DOUBLETS)) {
        printProgress(verbose, "Loading data")
        d <- discardRowCondDoublets(d)
    }

    if (isSymmetric(as.matrix(d)) & !header) {
        colnames(d) <- rownames(d)
    }

    return(d)
}


# avoid doublets in row names
# r: row names vector
renameRownameDoublets <- function(names.row) {
    j <- 1

    for (i in 2:length(names.row)) {
        if (names.row[i] == names.row[i - 1]) {
            j <- j + 1
            names.row[i] <- paste(names.row[i], ".", j, sep = "")
        } else {
            j <- 1
        }
    }
    return(names.row)
}

# rename row and avoid doublets errors
renameRowname <- function(d) {
    # names.row = as.character(d[,1])
    names.row <- rownames(d)
    d <- d[, -1]
    names.row <- renameRownameDoublets(names.row)

    tryCatch(
        {
            substr(names.row, 1, MAX_CHAR_LEN) -> rownames(d)
            return(d)
        },
        warning = function(w) {
            names.row <- renameRownameDoublets(substr(names.row, 1, MAX_CHAR_LEN))
            names.row -> rownames(d)
            return(d)
        },
        error = function(e) {
            return(d)
        }
    )
}

# Discard row from a reaction dataset that have the same conditions in each columns
# x: dataframe
discardRowCondDoublets <- function(x) {
    row_doublets <- list()
    j <- 0

    for (i in 1:nrow(x)) {
        # uniq remove doublets in a vector, so return 1 only if there is only 1
        if ((length(unique(as.integer(x[i, ]))) == 1)) {
            # print(row.names(x[i,]))
            j <- j + 1
            row_doublets[[j]] <- i
        }
    }

    if (length(row_doublets) != 0) {
        removed_reacs <- row.names(x[unlist(row_doublets), ])
        removed_conds <- x[unlist(row_doublets), 1]
        removed <- cbind(removed_reacs, removed_conds)
        colnames(removed) <- c("condition", "")
        assign("removed", removed, .GlobalEnv)
        writeTsv("removed", v = FALSE)
    }

    if (length(row_doublets) > 0) {
        return(x[-unlist(row_doublets), ])
    } else {
        return(x)
    }
}

# rename row and avoid doublets errors
renameRowname <- function(d) {
    names.row <- as.character(d[, 1])
    d <- d[, -1]
    names.row <- renameRownameDoublets(names.row)
    tryCatch(
        {
            substr(names.row, 1, 25) -> rownames(d)
            return(d)
        },
        warning = function(w) {
            names.row <- renameRownameDoublets(substr(names.row, 1, 25))
            names.row -> rownames(d)
            return(d)
        },
        error = function(e) {
            return(d)
        }
    )
}

# Discard row from a reaction dataset that have the same conditions in each columns
# x: dataframe
discardRowCondDoublets <- function(x) {
    row_doublets <- list()
    j <- 0
    for (i in 1:nrow(x)) {
        # uniq remove doublets in a vector, so return 1 only if there is only 1
        if ((length(unique(as.integer(x[i, ]))) == 1)) {
            # print(row.names(x[i,]))
            j <- j + 1
            row_doublets[[j]] <- i
        }
    }
    if (length(row_doublets) != 0) {
        removed_reacs <- row.names(x[unlist(row_doublets), ])
        removed_conds <- x[unlist(row_doublets), 1]
        removed <- cbind(removed_reacs, removed_conds)
        colnames(removed) <- c("condition", "")
        assign("removed", removed, .GlobalEnv)
        writeTsv("removed", v = FALSE)
    }
    if (length(row_doublets) > 0) {
        return(x[-unlist(row_doublets), ])
    } else {
        return(x)
    }
}

# Inputs: x : a matrix
# filename of the saved file
# Prints the matrix, save the matrix
writeTsv <- function(x, f = NULL, cl = FALSE, v = TRUE) {
    # print on stdout
    if (isTRUE(v)) {
        cat(paste("\n", gsub("_", " ", toupper(x)), ":\n", sep = ""))
    }
    # disabling warning
    options(warn = -1)
    # get variable
    tab <- base::get(x)
    if (!isTRUE(cl)) {
        output <- as.matrix(rbind(c("", colnames(tab)), cbind(rownames(tab), tab)))
    } else {
        output <- tab
    }
    # discard empty rows
    output <- output[rowSums(is.na(output)) != ncol(output), ]
    # TODOD:
    # output = output[,colSums(is.na(output)) != nrow(output)]
    output[is.na(output)] <- ""
    colnames(output) <- rep("", ncol(output))
    rownames(output) <- rep("", nrow(output))
    if (isTRUE(v)) {
        if (!isTRUE(cl)) {
            printed <- round(apply(output[-1, -1], 2, as.numeric), 2)
            rownames(printed) <- rownames(tab)
            colnames(printed) <- colnames(tab)
        } else {
            printed <- output
        }
        print(printed, quote = FALSE)
    }
    if (is.null(f)) {
        f <- paste(x, ".tsv", sep = "")
    }
    write(t(output), f, ncolumns = ncol(output), sep = "\t")
    options(warn = 0)
}

################################
#          Graphic
################################

# Usage: colPers(x), x a number of colours in output
# Gradient of color
colPers <- colorRampPalette(c(
    rgb(0.6, 0.1, 0.5, 1),
    rgb(1, 0, 0, 1),
    rgb(0.9, 0.6, 0, 1),
    rgb(0.1, 0.6, 0.3, 1),
    rgb(0.1, 0.6, 0.5, 1),
    rgb(0, 0, 1, 1)
),
alpha = TRUE
)


setGraphic <- function() {
    setGraphicBasic()
    par(mar = c(5.1, 5.1, 5.1, 2.1))
}

setGraphicBasic <- function() {
    par(
        cex.lab = 1.5,
        font.lab = 3,
        font.axis = 3,
        cex.axis = 0.8,
        cex.main = 2,
        cex = 1,
        lwd = 3
    )
}

plotAxis <- function(side, min, max, interval = 1, lwd = 3) {
    axis(side, seq(min, max, interval), lwd = lwd)
}

plotBestClustering <- function(
    sub_title,
    values,
    values_type,
    optimal_nb_clusters,
    n,
    interval = 1,
    min_x = 2,
    best = NULL,
    val2 = NULL,
    verbose = FALSE) {

    plotAxis(1, 2, n)

    if (interval >= 1) {
        axisSeq <- round(values)
    } else {
        axisSeq <- c(0, max(values) + 0.1)
    }

    # case of plotting gap statistics
    if (min_x < 2) {
        best_y <- values[optimal_nb_clusters]
    } # case of fusion levels
    else if (!is.null(val2)) {
        best_y <- values[optimal_nb_clusters - 1]
    } else {
        best_y <- max(values)
    }

    # for non-elbow plots
    if (!is.null(val2)) {
        best <- round(max(val2), 2)
    } else if (is.null(best)) {
        best <- round(max(values), 2)
    }

    plotAxis(2, min(axisSeq), max(axisSeq), interval)
    title(
        main = "Optimal number of clusters",
        line = 2,
        cex.main = 2
    )
    mtext(
        text = sub_title,
        font = 3,
        cex = 1.2,
        line = 0.5
    )
    abline(
        v = optimal_nb_clusters,
        col = "red",
        lty = 2,
        lwd = 2
    )
    points(
        optimal_nb_clusters,
        best_y,
        pch = 19,
        col = "red",
        cex = 2
    )

    if (!is.null(val2)) {
        t_values <- val2
    } else {
        t_values <- values
    }

    if (isTRUE(TEXT)) {
        text(
            y = values,
            x = min_x:n,
            labels = round(t_values, 2),
            cex = 1.2,
            pos = 4,
            col = "red"
        )
    }
    if (isTRUE(verbose)) {
        cat(
            "Optimal number of clusters k = ",
            optimal_nb_clusters,
            "\n",
            "With a",
            values_type,
            " of ",
            best,
            "\n",
            sep = ""
        )
    }
}

# f: filename
savePdf <- function(f) {
    pdf(f)
    setGraphic()
}

################################
#          Statistics
################################

# Get the normalized distance between each points and the center
# Outputs:
# for each column, the mean=0 and the variance is the same
scalecenter <- function(d) {
    # output scale function: for each column, mean=0, sd=1
    return(scale(d) * sqrt(nrow(d) / (nrow(d) - 1)))
    # ponderation for sampling index (var use n-1)
    # without this constante, for advanced outputs, total (max_cluster=nrow(data)) will be different from 1
}

# df: dataframe
# d: distance type
getDistance <- function(df, d) {
    dists <- c("euclidian", "manhattan", 1, 2, 5, 7)

    if (d < 3) {
        dist(df, method = dists[d])
    } else {
        dist.binary(df, method = as.integer(dists[d]))
    }
}

# d distance object
checkEuclidean <- function(d) {
    if (attributes(d)$method != "euclidean") {
        stop("Distance should be euclidean with this classification method.", call. = FALSE)
    }
}

isSymmetric <- function(d) {
    if (nrow(d) == ncol(d)) {
        isReflexivity <- unique(d[cbind(1:nrow(d), 1:nrow(d))] == 0)

        if (length(isReflexivity) == 1 & isTRUE(isReflexivity)) {
            isCommutativity <- unique(d[lower.tri(d)] == t(d)[lower.tri(d)])

            if (length(isCommutativity) == 1 & isTRUE(isCommutativity)) {
                return(TRUE)
            }
        }
    }
    return(FALSE)
}

################################
#          Clustering
################################

# Inputs:
# t: number of type of classification
# df: data
# d: distance matrix
# Ouput: Hierarchical classification
getCAH <- function(t, df, d) {
    if (t > 2) {
        if (t == 8 | t == 9) {
            checkEuclidean(d)
        }

        # cah: classification hierarchic ascending
        cah <- hclust(d, method = getClassifType(t))
        # automaticly ordering by clusters
        return(reorder.hclust(cah, d))
    }
}

# Selects best algo based on cophenetic calculation
# df: data
# d: distance matrix
selectBestCAH <- function(df, d, v = FALSE) {
    temp <- 0
    for (i in 3:9) {
        cah <- getCAH(df, i)
        res <- cor(d, cophenetic(cah))
        if (isTRUE(v)) {
            cat(paste(getClassifType(i), ":", round(res, 3), "\n"))
        }
        if (res > temp) {
            temp <- res
            t <- i
        }
    }
    cat(paste("Selected:", getClassifType(t), "\n"))
    return(t)
}

# Inputs:
# t: number of type of classification
getClassifType <- function(t) {
    methods <- c(
        "kmedoids",
        "kmeans",
        "ward.D2",
        "complete",
        "single",
        "average",
        "mcquitty",
        "median",
        "centroid"
    )
    methods[t]
}

# Agglomerative coefficient ()
getCoefAggl <- function(c) {
    coef.hclust(c)
}

# Inputs:
# t: number of type of classification
# df: data
# d: distance for pam
# k: number of clusterting
# Ouput: Non-hierarchical classification
getCNH <- function(t, df, d, k) {
    if (t == 1) {
        return(pam(d, k, diss = TRUE))
    } else if (t == 2) {
        checkEuclidean(d)
        return(kmeans(df, centers = k, nstart = 100))
    }
}

getClassif <- function(t, n, df, d) {
    if (t > 2) {
        getCAH(t, df, d)
    } else {
        list_cnh <- list("method" = getClassifType(t))
        for (k in 2:(n + 1)) {
            list_cnh[[k]] <- getCNH(t, df, d, k)
        }
        return(list_cnh)
    }
}

# Inputs:
# t: number of type of classification
# k: number of clusters
# c: hierarchical classification
# d: data
# Output: partitionning contening k clusters
getClusters <- function(k, c) {
    if (c$method == "kmedoids") {
        c[[k]]$clustering
    } else if (c$method == "kmeans") {
        c[[k]]$cluster
    } else {
        cutree(c, k)
    }
}

getClusterPerPart <- function(n, c) {
    cl <- list()
    for (k in 2:n) {
        cl[[k - 1]] <- getClusters(k, c)
    }
    return(cl)
}

# Input:
# cl: clusters
colorClusters <- function(cl) {
    NB_CLUSTERS <- length(levels(as.factor(cl)))
    for (i in 1:NB_CLUSTERS) {
        cl[cl == i] <- colPers(NB_CLUSTERS)[i]
    }
    return(cl)
}

# Inputs:
# cl: clusters
# f : filename
# r: ordered alphabetically
writeClusters <- function(f, sil_k, pca, v = FALSE) {
    cluster <- cbind(sil_k[, 1], sil_k[, 3], pca$li[attr(sil_k, "iOrd"), c(1, 2)])
    colnames(cluster) <- c("Cluster", "Silhouette", "Axis1", "Axis2")
    assign("cluster", cluster, .GlobalEnv)
    writeTsv("cluster", f, cl = FALSE, v = v)
}

############################################################
#          Cophenetic (dendrogram distance matrix)
############################################################

# Distance matrix between each leaf of the dendogramm
# Inputs:
# d : distance matrix
# cah : hierarchical classification
plotCohenetic <- function(d, cah, is_png = FALSE, verbose = FALSE) {
    coph_matrix <- cophenetic(cah)
    cor_coph <- cor(d, coph_matrix)
    if (isTRUE(verbose)) {
        cat(
            paste(
                "\nCOPHENETIC:\nExplained variance (%):",
                round(cor_coph^2, 3),
                "\nCorrelation with the data:",
                round(cor_coph, 3),
                "\n"
            )
        )
    }

    # if(is_png) {
    #   #png(paste(opt$output8, ".png", sep=""), DIM_PNG/2, DIM_PNG/2)
    #   par(cex.lab=1.5*2, font.lab=3, font.axis=3, cex.axis=0.8*2, cex.main=2*2, cex=1, lwd=3*2)
    #   par(mar=c(5.1,5.1,5.1,2.1)+7)
    #   lwd=3*2
    #   line.lab = 5
    # }else{
    setGraphic()
    # savePdf(paste(opt$output8, ".pdf", sep=""))
    lwd <- 3
    line.lab <- 3
    # }

    plot(
        d,
        coph_matrix,
        pch = 19,
        col = alpha("red", 0.2),
        axes = FALSE,
        xlim = c(0, max(d)),
        xlab = "",
        ylab = "",
        ylim = c(0, max(coph_matrix)),
        asp = 1,
        main = paste("Cophenetic correlation: ", round(cor_coph, 3))
    )
    title(
        xlab = "Distance between metabolites",
        ylab = "Cophenetic distance",
        line = line.lab
    )
    plotAxis(2, 0, max(coph_matrix), lwd = lwd)
    plotAxis(1, 0, max(d), lwd = lwd)
    abline(0, 1, col = "grey", lty = 2, lwd = lwd)
    # suprLog = dev.off()
}

##############################################
#          Inertia
##############################################

# Relative inter-group inertia for each partitionning
# Inputs:
# n: maximum number of clusters
# d: dataframe
# cl: list of clusters per partition
getRelativeBetweenPerPart <- function(n, d, cl) {
    d <- as.matrix(d)
    between <- rep(0, n - 1)
    # total sum of square
    TSS <- sum(scale(d, scale = FALSE)^2)
    for (i in 2:n) {
        cl_k <- as.factor(cl[[i - 1]])
        # apply(d, 2, mean) : centroids for each column
        # as.vector(table(cl) : size of each clusters
        # t : vector rotation for arithmetic with other row or column vectors
        between[i - 1] <- sum(t((t(getClusterCentroids(d, cl_k)) - apply(d, 2, mean))^2) * as.vector(table(cl_k))) / TSS
    }
    return(100 * between)
}

# tapply(data[,i], Cla, mean) :
# centroids of each clusters for a column i
# sapply(1:ncol(data), function(i) tapply(data[,i], Cla, mean)) :
# centroids of each clusters for each column
getClusterCentroids <- function(d, cl) {
    sapply(1:ncol(d), function(i) tapply(d[, i], cl, mean))
}

# Difference between each case of a vector
getBetweenDifferences <- function(between) {
    # apply produce a list, unlist convert in vector
    diff <- unlist(sapply(1:length(between), function(i) between[i] - between[i - 1]))
    return(as.vector(cbind(between[1], t(diff))))
    #-n-1 to remove the last NA value (pairwise comparison)
    # between[1] to get the difference with 1 cluster
}

getWithin <- function(d, cl, k) {
    nk <- length(cl[cl == k]) # number of individuals in the cluster
    d1 <- scalecenter(d)
    return(nk * sum(getClusterCentroids(d1, cl)[k, ]^2) / nrow(d))
}

# cl: list of clusters per partition
getRelativeWithinPerCluster <- function(cls, d) {
    n <- length(cls)
    within <- matrix(NA, n - 1, n)
    rownames(within) <- seq(2, n)
    colnames(within) <- paste("G", seq(1, n), sep = "")
    for (k in 2:n) {
        cl <- cls[[k - 1]]
        for (i in 1:length(table(cl))) {
            within[k - 1, i] <- getWithin(d, cl, i)
        }
        within[k - 1, ] <- within[k - 1, ] / sum(as.numeric(na.omit(within[k - 1, ])))
    }
    return(within)
}

# Between inertia differences between a partionning and the previous
plotBetweenDiff <- function(between_diff, verbose = FALSE) {
    if (isTRUE(verbose)) {
        cat("\nBETWEEN DIFFERENCES:\n")
    }
    optimal_nb_clusters <- which.max(between_diff) + 1
    n <- length(between_diff)
    # savePdf("between_differences.pdf")
    setGraphic()
    plot(
        2:(n + 1),
        between_diff,
        type = "b",
        ylim = c(round(min(between_diff)) - 1, round(max(between_diff)) + 1),
        xlim = c(2, (n + 2)),
        xlab = "Nb. of clusters",
        ylab = "Between-cluster variation (%)",
        col = "grey",
        axes = FALSE
    )
    plotBestClustering(
        "Largest between differences method",
        between_diff,
        " variation with the previous partitionning (%)",
        optimal_nb_clusters,
        verbose = verbose,
        n = n
    )
    # suprLog = dev.off()
}

plotFusionLevels <- function(n, c, verbose = FALSE) {
    if (isTRUE(verbose)) {
        cat("\nFUSION LEVELS:\n")
    }
    fusion <- rev(c$height)
    diff <- unlist(sapply(1:n, function(i) fusion[i - 1] - fusion[i]))
    fusion <- fusion[1:(n - 1)]
    optimal_nb_clusters <- which.max(diff) + 1
    setGraphic()
    # savePdf(paste(opt$output9, ".pdf", sep=""))
    plot(
        2:n,
        fusion,
        type = "b",
        ylim = c(round(min(fusion)) - 1, round(max(fusion)) + 1),
        xlim = c(2, n + 1),
        xlab = "Nb. of clusters",
        ylab = "Cophenetic distance",
        col = "grey",
        axes = FALSE
    )
    plotBestClustering(
        "Fusion level method",
        fusion,
        " gain with the previous fusion level",
        optimal_nb_clusters,
        val2 = diff,
        verbose = verbose,
        n = n
    )
    # suprLog = dev.off()
}

# x: vector of between inertia for k partitions
plotElbow <- function(x, verbose = FALSE) {
    if (isTRUE(verbose)) {
        cat("\nELBOW:\n")
    }
    n <- length(x) + 1
    within <- c(100, 100 - x)
    ratio <- within[1:(n - 1)] / within[2:n]
    best <- round(min(ratio), 2)
    optimal_nb_clusters <- which.min(ratio)
    setGraphic()
    # savePdf("elbow.pdf")
    plot(
        1:n,
        within,
        type = "b",
        ylim = c(min(within) - 1, 101),
        xlim = c(1, n + 1),
        xlab = "Nb. of clusters",
        ylab = "Relative within inertia (%)",
        col = "grey",
        axes = FALSE
    )
    plotBestClustering(
        "Elbow method",
        within,
        " Wk/Wk+1 ratio ",
        optimal_nb_clusters,
        5,
        1,
        best,
        n = n,
        verbose = verbose
    )
    # suprLog = dev.off()
}

################################
#          Silhouette
################################

# Ouput: an ordered silhouette object
getSilhouette <- function(df, cl_k, d) {
    s <- sortSilhouette(silhouette(cl_k, d))
    rownames(s) <- row.names(df)[attr(s, "iOrd")]
    return(s)
}

getSilhouettePerPart <- function(df, cl, d) {
    list_sil <- list()
    for (k in 2:length(cl)) {
        list_sil[[k - 1]] <- getSilhouette(df, cl[[k - 1]], d)
    }
    return(list_sil)
}

# sils: list of silhouettes objects per partition
getMeanSilhouettePerPart <- function(sils) {
    unlist(sapply(1:length(sils), function(i) summary(sils[[i]])$avg.width))
}

# Plot the best average silhouette width for all clustering possible
# mean_sils: vector of silhouette average width
plotSilhouettePerPart <- function(mean_silhouette, sil = sil, verbose = FALSE) {
    if (isTRUE(verbose)) {
        cat("\nSILHOUETTE:\n")
    }
    setGraphic()
    # savePdf(opt$output1)
    optimal_nb_clusters <- which.max(mean_silhouette) + 1
    n <- length(sil)
    plot(
        2:(n + 1),
        mean_silhouette,
        type = "b",
        xlim = c(2, n + 2),
        ylim = c(0, max(mean_silhouette) + 0.1),
        col = "grey",
        xlab = "Nb. of clusters",
        ylab = "Average silhouette width",
        axes = FALSE
    )
    plotBestClustering(
        "Silhouette method",
        mean_silhouette,
        "n average width",
        optimal_nb_clusters,
        0.1,
        verbose = verbose,
        n = (n + 1)
    )
    # suprLog = dev.off()
    # return (optimal_nb_clusters)
}

# sil_k: a silhouette object
plotSilhouette <- function(sil_k) {
    # pdf(opt$output2)
    setGraphicBasic()
    par(mar = c(4, 12, 3, 2))
    plot(
        sil_k,
        max.strlen = MAX_CHAR_LEN,
        main = " ",
        sub = "",
        do.clus.stat = TRUE,
        xlab = "Silhouette width",
        cex.names = 0.8,
        col = colorClusters(sil_k[, 1]),
        nmax.lab = 100,
        do.n.k = FALSE,
        axes = FALSE
    )
    mtext(
        paste("Average silhouette width:", round(summary(sil_k)$avg.width, 3)),
        font = 2,
        cex = 1.5,
        line = 1
    )
    plotAxis(1, 0, 1, 0.2)
    # suprLog = dev.off()
}

###################################
#          GAP STATISTICS
###################################

# B: nb of NB_BOOTSTRAP
getGapPerPart <- function(n, d, c, B = 500, v = FALSE) {
    # FUN mus have only two args in this order and return a list with an object cluster

    if (isTRUE(v)) {
        cat("\nGAP STATISTICS:\n")
    }
    if (c$method == "kmeans" & (n > 10 | nrow(d) >= 100)) {
        plural <- c("few ", "s")
    } else {
        plural <- c("", "")
    }
    if (c$method == "kmeans" |
        nrow(d) >= 100) {
        cat(paste("It could take a ", plural[1], "minute", plural[2], "...\n", sep = ""))
    }

    gapFun <- function(x, k) list(cluster = getClusters(k, c))
    clusGap(d, FUNcluster = gapFun, K.max = n, verbose = FALSE, B = B)
}

# g: gap object
getGapBest <- function(g, M = "Tibs2001SEmax") {
    with(g, maxSE(Tab[, "gap"], Tab[, "SE.sim"], method = M))
}

# Plot the gap statistics width for all clustering possible
# TODO: HERE
plotGapPerPart <- function(g, n, verbose = FALSE) {
    setGraphic()
    # savePdf("gap_statistics.pdf")
    optimal_nb_clusters <- getGapBest(g)
    gap_k <- round(g$Tab, 3)
    best <- gap_k[, "gap"][optimal_nb_clusters]
    if (optimal_nb_clusters < n) {
        best <- paste(best, ">", gap_k[, "gap"][optimal_nb_clusters + 1], "-", gap_k[, "SE.sim"][optimal_nb_clusters + 1])
    }
    plot(
        g,
        arrowArgs = list(
            col = "gray",
            length = 1 / 15,
            lwd = 2,
            angle = 90,
            code = 3
        ),
        type = "b",
        xlim = c(1, n + 1),
        ylim = c(0, max(g$Tab[, "gap"]) + 0.1),
        col = "grey",
        xlab = "Nb. of clusters",
        ylab = expression(Gap[k]),
        main = "",
        axes = FALSE
    )
    plotBestClustering(
        "Gap statistics method",
        g$Tab[, "gap"],
        " gap value",
        optimal_nb_clusters,
        0.1,
        1,
        best,
        verbose = verbose,
        n = n
    )
    # cat(paste("With a corrected index, optimal number of clusters k =",getGapBest(gap,"firstSEmax"), "\n"))
    # suprLog = dev.off()
}

# Plot the gap between the two function: within and random within average
plotGapPerPart2 <- function(g, n) {
    setGraphic()
    # savePdf("log_w_diff.pdf")
    min_y <- round(min(g$Tab[, c(1, 2)]), 1)
    max_y <- round(max(g$Tab[, c(1, 2)]), 1)
    plot(
        0,
        0,
        xlim = c(1, n),
        ylim = c(min_y - 0.1, max_y + 0.1),
        type = "n",
        xlab = "Nb. of clusters",
        ylab = "log(within-inertia)",
        axes = FALSE
    )
    title(
        main = "Optimal number of clusters",
        line = 2,
        cex.main = 2
    )
    mtext(
        text = "Gap statistics method",
        font = 3,
        cex = 1.2,
        line = 0.5
    )
    optimal_nb_clusters <- getGapBest(g)
    abline(
        v = optimal_nb_clusters,
        col = "gray",
        lty = 2,
        lwd = 2
    )
    lines(seq(1:n), g$Tab[, 1], type = "b", col = "red")
    lines(seq(1:n), g$Tab[, 2], type = "b", col = "blue")
    plotAxis(1, 1, n)
    plotAxis(2, min_y, max_y, 0.1)
    legend(
        "topright",
        c("log(W)", "E.log(W)"),
        col = c("red", "blue"),
        lty = 1,
        box.lwd = -1,
        bg = "white"
    )
    # suprLog = dev.off()
}

printSummary <- function(between, diff, sil, adv, gap = NULL) {
    # TODO: no n = nrow(data)
    summary <- cbind(between, diff, 100 - between, sil)
    names <- c(
        "Between-inertia (%)",
        "Between-differences (%)",
        "Within-inertia (%)",
        "Silhouette index"
    )
    if (!is.null(gap)) {
        summary <- cbind(summary, gap$Tab[, "gap"][-1], gap$Tab[, "SE.sim"][-1])
        names <- c(names, "Gap", "Gap SE")
    }
    rownames(summary) <- seq(2, (length(between) + 1))
    colnames(summary) <- names
    return(summary)
}

################################
#          HEATMAP
################################

# Inputs:
# cl_size: vector of size for each clusters
plotRect <- function(cl_sizes, colors, lwd = 3) {
    # size of each clusters
    temp_size <- 0
    for (i in 1:length(cl_sizes)) {
        # y begin at the top, so sum(cl_sizes) must be substracted to y coord.
        # rect(xleft, ybottom, xright, ytop)
        # +0.5 because x, y coord are shifted to 0.5 comparativly to plotcolors functions
        rect(
            temp_size + 0.5,
            sum(cl_sizes) - temp_size - cl_sizes[i] + 0.5,
            cl_sizes[i] + temp_size + 0.5,
            sum(cl_sizes) - temp_size + 0.5,
            border = colors[i],
            lwd = lwd
        )
        # memorize the size of the cluster (for a bottom-right shift)
        temp_size <- temp_size + cl_sizes[i]
    }
}

# Outputs:
# lenght of clusters ordered by the clusters order
getOrderedClusterSize <- function(cl) {
    nb_cl <- length(levels(as.factor(cl)))
    size_cl <- rep(0, nb_cl)
    temp_cl <- rep(0, length(cl))
    j <- 0

    for (i in 1:length(cl)) {
        if (!cl[i] %in% temp_cl) {
            j <- j + 1
        }
        size_cl[j] <- size_cl[j] + 1
        temp_cl[i] <- cl[i]
    }
    return(size_cl)
}

# Inputs:
# df: data frame
# d: a distance object
# s: an organised silhouette object
# c: CAH
# c: clusters from CAH
heatMap <- function(df, d, s = NULL, c = NULL, cl = NULL, is_png = FALSE, verbose = FALSE) {
    printProgress(verbose, "Heatmap calculation")
    text <- isTRUE(isTRUE(TEXT) & (nrow(data) < 100))

    if (!is.null(s)) {
        order <- attr(s, "iOrd")
        cl_sizes <- summary(s)$clus.size
        title <- "silhouette\'s scores"
        colors <- colPers(length(levels(as.factor(s[, 1]))))
    } else {
        order <- c$order
        cl_sizes <- getOrderedClusterSize(cl[order])
        title <- "dendrogram"
        colors <- orderColors(c, cl)
    }

    matrix <- as.matrix(d)
    matrix <- matrix[order, order]
    rownames(matrix) <- rownames(df)[order] -> labels
    # if(tri == TRUE) matrix[!lower.tri(matrix)] = NA
    # image(1:ncol(matrix), 1:ncol(matrix), t(matrix), axes=FALSE, xlab="", ylab="")

    options(warn = -1)
    if (nrow(df) > NB_ROW_MAX) {
        labels <- order
    }
    # png(opt$output4,DIM_PNG, DIM_PNG)
    if (is_png) {
        cex.main <- 5
        cex.legend <- 3
        cex.lab <- 2
        y_top <- 12
        x_lab <- 0.6
        lwd.rect <- 6
    } else {
        # pdf(opt$output4)
        cex.main <- 1.5
        cex.legend <- 0.85
        cex.lab <- 0.7
        y_top <- 8
        x_lab <- 0.5
        lwd.rect <- 3
    }

    par(fig = c(0, 0.9, 0, 1), new = TRUE)
    par(mar = c(1, 8, y_top, 1))
    plotcolors(
        dmat.color(matrix, colors = heat.colors(1000), byrank = FALSE),
        ptype = "image",
        na.color = "red",
        rlabels = FALSE,
        clabels = FALSE,
        border.color = 0
    )
    mtext(
        paste("Distance matrix ordered by", title),
        3,
        line = 6,
        font = 4,
        cex = cex.main
    )
    text(
        -0.5,
        0:(ncol(matrix) - 1) + 1,
        rev(labels),
        xpd = NA,
        adj = 1,
        cex = 0.7
    )
    text(
        0.5:(ncol(matrix) - 0.5),
        ncol(matrix) + 1,
        substr(labels, 0, 20),
        xpd = NA,
        cex = 0.7,
        srt = 65,
        pos = 4
    )
    plotRect(cl_sizes, colors, lwd.rect)
    if (isTRUE(text)) {
        text(expand.grid(1:ncol(matrix), ncol(matrix):1), sprintf("%d", matrix), cex = 0.4)
    }

    par(fig = c(0.85, 1, 0.3, 0.8), new = TRUE)
    par(mar = c(5, 0, 4, 0) + 0.1)
    legend_image <- as.raster(matrix(heat.colors(1000), ncol = 1))
    plot(
        c(0, 1),
        c(0, 1),
        type = "n",
        axes = FALSE,
        xlab = "",
        ylab = "",
        main = ""
    )
    rasterImage(legend_image, 0.4, 0, 0.5, 1)
    mtext("   Distance",
        3,
        line = 0.5,
        cex = cex.legend,
        font = 2
    )
    text(
        x = x_lab,
        y = seq(0, 1, l = 3),
        labels = round(seq(max(matrix), 2, l = 3)),
        cex = cex.lab,
        pos = 4
    )

    options(warn = 0)
    if (verbose) {
        cat("done.\n")
    }
    # suprLog = dev.off()
}

################################
#          Dendrogram
################################

# Inputs:
# k: number of clusters
plotDendrogram <- function(t, k, c, d, n, cl) {
    if (nrow(d) > NB_ROW_MAX) {
        c$labels <- 1:nrow(d)
        cex <- 0.4
    } else {
        cex <- 0.8
    }
    # pdf(opt$output7)
    setGraphicBasic()
    par(mar = c(2, 5, 5, 1))
    plot(
        c,
        hang = -1,
        ylim = c(0, max(c$height)),
        xlim = c(0, length(c$labels)),
        sub = "",
        cex = cex,
        font = 3,
        ylab = "Cophenetic distance",
        main = "Dendrogram",
        axes = FALSE
    )
    plotAxis(2, 0, max(c$height))
    abline(h = rev(c$height)[1:n], col = "gray", lty = 2, lwd = 1)
    # projection of the clusters
    rect.hclust(c, k = as.numeric(k), border = orderColors(c, cl))
    # suprLog = dev.off()
}

# Get colors ordered for dendrogram
orderColors <- function(c, cl) {
    col_in <- colorClusters(cl)[c$order]
    j <- 1
    col_ordered <- rep(NA, length(table(cl)))
    col_ordered[1] <- col_in[1]
    for (i in 2:length(col_in)) {
        if (col_in[i] != col_in[i - 1]) {
            j <- j + 1
            col_ordered[j] <- col_in[i]
        }
    }
    # vector of color: one by cluster
    return(col_ordered)
}

################################
#            PCA
################################

# nf: number of factorial axis
plotPca <- function(pca, d, cl, axis1 = 1, axis2 = 2, advanced = FALSE, is_png = FALSE) {
    k <- length(levels(as.factor(cl)))

    if (nrow(d) > NB_ROW_MAX) {
        cpoint <- 0
        cstar <- 0
        cellipse <- 0
        clabel <- 0
        labels <- 1:nrow(d)
    } else {
        cpoint <- 0
        cstar <- 1
        cellipse <- 1
        clabel <- 0
        labels <- rownames(d)
    }

    if (is_png) {
        par(mar = c(0, 0, 18, 0), lwd = 4)
        cex <- 2
        cex.main <- 6
        lwd.line <- 8
        line.main <- 7
    } else {
        # pdf(opt$output3)
        par(mar = c(0, 0, 4.1, 0))
        cex <- 0.8
        cex.main <- 1.5
        lwd.line <- 2
        line.main <- 1
    }

    title <- paste("Cumulated inertia:", round((pca$eig[axis1] + pca$eig[axis2]) / sum(pca$eig), 4) * 100, "%")
    s.class(
        addaxes = FALSE,
        cbind(pca$li[, axis1], pca$li[, axis2]),
        ylim = c(min(pca$li[, axis2]) + min(pca$li[, axis2]) / 4, max(pca$li[, axis2]) + max(pca$li[, axis2]) / 2),
        xlim = c(min(pca$li[, axis1]), max(pca$li[, axis1])),
        csub = 1.5,
        as.factor(cl),
        grid = FALSE,
        col = colPers(k),
        clabel = clabel,
        cstar = cstar,
        cellipse = cellipse,
        cpoint = cpoint
    )
    mtext(title, font = 2, line = line.main, cex = cex.main)
    abline(h = 0, v = 0, lty = 2, lwd = lwd.line, col = "grey")
    text(
        x = pca$li[, axis1],
        y = pca$li[, axis2],
        labels = labels,
        col = colorClusters(cl),
        cex = cex
    )
    # colnames(pca_coord) = c("Chemicals", "Axis 1", "Axis 2")

    if (isTRUE(advanced)) {
        par(fig = c(0.8, 1, 0.82, 1), new = TRUE)
        plotInertiaPca(pca, d, pca$nf)
    }
    # suprLog = dev.off()
}

# nf: number of inertia bar plot corresponding to factorial axis
plotInertiaPca <- function(pca, d, nf = 4) {
    if (nrow(d) > NB_ROW_MAX) {
        r_lim <- c(8, 0, 4, 5)
        r_main_cex <- 2.7
        r_main_text <- 2.4
        lwd.hist <- 40
        line.hist <- 2
    } else {
        # r_lim = c(-0.2, 0.3, 1.1, 1.1);
        r_lim <- c(2, 0, 1, 1)
        r_main_cex <- 0.7
        r_main_text <- 0.6
        lwd.hist <- 10
        line.hist <- 0
    }

    inertia <- round(pca$eig / sum(pca$eig) * 100, 1)
    par(mar = c(r_lim[1], r_lim[2], r_lim[3], r_lim[4]) + 0.1)
    plot(
        inertia,
        type = "h",
        lwd = lwd.hist,
        lend = 1,
        xlim = c(0, nf + 0.2),
        ylim = c(0, max(inertia + 7)),
        col = "grey75",
        font = 2,
        axes = FALSE,
        xlab = "",
        ylab = ""
    )
    title(
        sub = " Inertia (in %)",
        line = line.hist,
        cex.sub = r_main_cex,
        font.sub = 3
    )
    text(1:nf, inertia[1:nf] + 5, inertia[1:nf], cex = r_main_text)
    par(new = TRUE)
    par(mar = c(0, 0, 0, 0))
    plot(0:1, 0:1, axes = FALSE, type = "n")
    rect(0, 0.1, 0.9, 0.9, border = "grey65")
}

#########################################
#            Variables contribution
#########################################

# For a given partition (cl) and each variables (dataset columns)
# pondered distance between the centroid of each clusters and the global centroid of the cloud
# Inputs:
# d: data
# cl: clusters object
getDistPerVariable <- function(d, cl) {
    # Distance between the centroid of each variables
    # ponderation by the sd of the variable (=total inertia per var)
    d <- scalecenter(d)
    nb_cl <- length(levels(as.factor(cl)))
    ctr <- matrix(0, nrow = nb_cl, ncol = ncol(d))

    for (i in 1:nrow(d)) {
        # get the group number for each row
        cli <- cl[i]
        # in the dataset, for a metabolite row, loop an each metadabolite column
        # values are affected the corresponding cluster row and metabolite column in ctr
        for (j in 1:ncol(d)) {
            ctr[cli, j] <- ctr[cli, j] + d[i, j]
        }
    }

    return(ctr)
}

# For a given partition, relative contributions of each metabolites to inertia of each clusters (CTR)
# The total of the clusters for each column corresponds to PDIS
# Inputs:
# t: number of type of classification
# k: number of clusters
# c: hierarchical classification
# d: data
getCtrVar <- function(t, k, cl, d) {
    # if NA values appear, scale 0/0 could produce NA values, NA could correspond to 0
    nb_cl <- length(levels(as.factor(cl)))
    ncol <- ncol(d)
    ctr <- getDistPerVariable(d, cl)

    rownames(ctr) <- paste("G", seq(1, k), sep = "")
    colnames(ctr) <- colnames(d)

    for (i in 1:nb_cl) {
        for (j in 1:ncol(d)) {
            ctr[i, j] <-
                ctr[i, j]^2 / (nrow(d) * length(cl[cl == i]))
        }
    }

    return(ctr)
}

getCtrVar2 <- function(t, k, cl, d, scale = TRUE) {
    ctr <- getCtrVar(t, k, cl, d)

    if (isTRUE(scale)) {
        ctr_part <- getPdis(t, k, cl, d)

        for (i in 1:nb_cl) {
            for (j in 1:ncol(d)) {
                if (scale) {
                    ctr[i, j] <- ctr[i, j] / ctr_part[j]
                }
            }
        }
    }

    return(ctr)
}

# Discriminant power (PDIS)
# Relative contributions of the metabolites to inertia of a partitionning (in %)
# Inputs:
# t: number of type of classification
# k: number of clusters
# c: hierarchical classification
# d: data
getPdis <- function(t, k, cl, d) {
    # for each metabolite contribution (in column), sum the k clusters values
    return(apply(getCtrVar(t, k, cl, d), 2, sum))
}

# Inputs:
# t: number of type of classification
# n: number max of clusters
# cls: list of clusters
# d: data
# index: pdis or rho2 calculation
getPdisPerPartition <- function(t, n, cls, d) {
    pdis_per_partition <- matrix(NA, n - 1, ncol(d))
    rownames(pdis_per_partition) <- seq(2, n)
    colnames(pdis_per_partition) <- colnames(d)

    for (k in 2:n) {
        res <- getPdis(t, k, cls[[k - 1]], d)
        for (i in 1:length(res)) {
            pdis_per_partition[k - 1, i] <- res[i]
        }
    }
    return(pdis_per_partition)
}
