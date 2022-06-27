#' The application server-side
#'
#' @param input,output,session Internal parameters for {shiny}.
#' @noRd
app_server <- function(input, output, session) {
    options(shiny.maxRequestSize = 30 * 1024^2)
    # Global variables settings
    NB_BOOTSTRAP <- 500 # should be comprise between 100 and 1000
    TEXT <- T # print values on graph (for optimum partition and heatmap)
    NB_ROW_MAX <- 200 # max row to have pdf, otherwise, some plots are in png
    DIM_PNG <- 2000
    VERBOSE_NIV2 <- F
    VERBOSE <- F
    MAX_CHAR_LEN <- 25 # maximum length of individual s names
    PNG <- F
    classif_methods <- list("K-menoids" = 1, "K-means" = 2, "Ward" = 3, "Complete links" = 4, "Single links" = 5, "UPGMA" = 6, "WPGMA" = 7, "WPGMC" = 8, "UPGMC" = 9)

tryCatch(
    {
        data <- loadData("matrix.txt")
    },
    warning = function(w) {
        data <- NULL
        warning = function(w) {
            data <- NULL
            message("Default file \"matrix.txt\" is not in the folder. Please, load another one.")
        },
        error = function(e) {
            data <- NULL
        }
    )

    loadData <- function(f, s = "\t", h = F) {
        # !file.exists(
        if (!is.null(f)) {
            if (grepl("xlsx?", f)) {
                data <- openxlsx::read.xlsx(f, rowNames = TRUE)
            } else {
                data <- read.table(f,
                    header = h,
                    sep = s,
                    dec = ".",
                    row.names = 1
                )
            }
            # data = preProcessData(data)
            # substr(rownames(data), 1, MAX_CHAR_LEN) -> rownames(data)
        }
        return(data)
    }

    refresh <- reactiveValues()
    refresh$classif_type <- ""
    getClassifValue <- function(key) unlist(classif_methods[key])

    ###################################
    #          SETTINGS
    ###################################

    setVariables <- reactive({
        refresh2 <- c(refresh$classif_type, input$max_biomark)

        assign(
            "data",
            refresh$data,
            .GlobalEnv
        )

        assign(
            "dis",
            refresh$dis,
            .GlobalEnv
        )

        assign(
            "MAX_CLUSTERS",
            refresh$max,
            .GlobalEnv
        )
        assign(
            "CLASSIF_TYPE",
            refresh$classif_type,
            .GlobalEnv
        )

        assign(
            "NB_CLUSTERS",
            input$nb_clusters,
            .GlobalEnv
        )
        assign(
            "ADVANCED",
            input$advanced,
            .GlobalEnv
        )

        assign(
            "AXIS1",
            input$axis1,
            .GlobalEnv
        )
        assign(
            "AXIS2",
            input$axis2,
            .GlobalEnv
        )
    })

    # Each shiny server functions run in local environment
    # With assign, variables are forced to be in global env
    setClassif <- reactive({
        setVariables()
        req(input$infile)

        if ((nrow(data) > 3000) & (input$classif_type > 2)) message("[WARNING] With more than 3000 rows to analyse, classification method must be K-medoids or K-means", call. = FALSE)

        # Perform classification
        printProgress(VERBOSE_NIV2, "Classification")

        assign(
            "classif",
            getClassif(CLASSIF_TYPE, MAX_CLUSTERS, data, dis),
            .GlobalEnv
        )
        if (VERBOSE_NIV2) cat("done.\n")
        assign(
            "list_clus",
            getClusterPerPart(MAX_CLUSTERS + 1, classif),
            .GlobalEnv
        )

        # inertia
        assign(
            "between",
            getRelativeBetweenPerPart(MAX_CLUSTERS, data, list_clus),
            .GlobalEnv
        )
        assign(
            "diff",
            getBetweenDifferences(between),
            .GlobalEnv
        )

        printProgress(VERBOSE_NIV2, "PCA")
        assign(
            "pca",
            dudi.pca(data, scannf = F, scale = input$scale, nf = 4),
            .GlobalEnv
        )
        if (VERBOSE_NIV2) cat("done.\n")

        assign(
            "sil",
            getSilhouettePerPart(data, list_clus, dis),
            .GlobalEnv
        )
        assign(
            "mean_silhouette",
            getMeanSilhouettePerPart(sil),
            .GlobalEnv
        )

        setClusters()

        setPrintFuncs()
    })

    setClusters <- reactive({
        refr <- c(refresh$classif_type, input$advanced, input$nb_clusters, input$scale, input$transpose)

        if (NB_CLUSTERS > 0) {
            assign(
                "optimal_nb_clusters",
                input$nb_clusters,
                .GlobalEnv
            )
        } else {
            assign(
                "optimal_nb_clusters",
                which.max(mean_silhouette) + 1,
                .GlobalEnv
            )
        }

        # cl_temp, because writeTsv(clusters) recreate a different object named clusters
        assign(
            "clusters",
            list_clus[[optimal_nb_clusters - 1]],
            .GlobalEnv
        )

        if (!input$advanced) {
            assign(
                "gap",
                NULL,
                .GlobalEnv
            )
        }

        assign(
            "sil_k",
            sil[[optimal_nb_clusters - 1]],
            .GlobalEnv
        )

        # errors
        if (optimal_nb_clusters == MAX_CLUSTERS) {
            message("\n[WARNING] The optimal number of clusters equals the maximum number of clusters. \nNo cluster structure has been found.")
        }
        if (min(table(clusters)) == 1) {
            message("\n[WARNING] A cluster with an only singleton biased the silhouette score.")
        }

        writeClusters("clusters.tsv", v = F)
    })

    setPrintFuncs <- function() {
        assign(
            "ADVANCED",
            input$advanced,
            .GlobalEnv
        )
        refr <- c(input$nb_clusters, input$max_biomark)

        ###### plot funcs #####
        assign(
            "plotPCA",
            function() plotPca(pca, data, clusters, AXIS1, AXIS2),
            .GlobalEnv
        )
        assign(
            "plotBest",
            function() plotSilhouettePerPart(mean_silhouette),
            .GlobalEnv
        )
        assign(
            "plotSil",
            function() plotSilhouette(sil_k),
            .GlobalEnv
        )
        assign(
            "plotDend",
            function() plotDendrogram(CLASSIF_TYPE, optimal_nb_clusters, classif, data, MAX_CLUSTERS, clusters),
            .GlobalEnv
        )

        if (CLASSIF_TYPE <= 2 | isTRUE(ADVANCED)) {
            assign(
                "plotHeatmap",
                function() heatMap(data, dis, sil_k),
                .GlobalEnv
            )
        } else {
            assign(
                "plotHeatmap",
                function() heatMap(data, dis, c = classif, cl = clusters),
                .GlobalEnv
            )
        }

        ##### advanced #####
        if (isTRUE(ADVANCED)) {
            assign(
                "plotFus",
                function() plotFusionLevels(MAX_CLUSTERS, classif),
                .GlobalEnv
            )
            assign(
                "plotCoph",
                function() {
                    printProgress(VERBOSE_NIV2, "Cophenetic calculation")
                    plotCohenetic(dis, classif)
                    if (VERBOSE_NIV2) cat("done.\n")
                },
                .GlobalEnv
            )
            assign(
                "plotGap",
                function() {
                    if (nrow(data) < (NB_ROW_MAX / 2)) {
                        plotGapPerPart(gap, MAX_CLUSTERS, v = F)
                        # plotGapPerPart2(gap, MAX_CLUSTERS)
                    } else {
                        message("\n[WARNING] Dataset too big to calculate a gap statistics.")
                    }
                },
                .GlobalEnv
            )
            assign(
                "plotGap2",
                function() {
                    if (nrow(data) < (NB_ROW_MAX / 2)) {
                        plotGapPerPart2(gap, MAX_CLUSTERS)
                    } else {
                        message("\n[WARNING] Dataset too big to calculate a gap statistics.")
                    }
                },
                .GlobalEnv
            )
            assign(
                "plotElb",
                function() plotElbow(between),
                .GlobalEnv
            )
            assign(
                "within_k",
                getRelativeWithinPerCluster(list_clus, data),
                .GlobalEnv
            )
        }

        ##### print table func #####

        assign(
            "summary",
            printSummary(between, diff, mean_silhouette, ADVANCED, gap),
            .GlobalEnv
        )
        assign(
            "ctr_part",
            100 * getPdisPerPartition(CLASSIF_TYPE, MAX_CLUSTERS, list_clus, data),
            .GlobalEnv
        )
        assign(
            "centroids",
            getDistPerVariable(data, clusters),
            .GlobalEnv
        )
        assign(
            "discr",
            getDiscriminantVariables(CLASSIF_TYPE, optimal_nb_clusters, clusters, data, input$max_biomark),
            .GlobalEnv
        )
        assign(
            "ctr_clus_plot",
            function() plotDiscriminantVariables(discr),
            .GlobalEnv
        )
        assign(
            "ctr_clus",
            100 * getCtrVar(CLASSIF_TYPE, optimal_nb_clusters, clusters, data),
            .GlobalEnv
        )
        writeTsv("discr", "discr_var.tsv", v = F)
    }

    # post-process for data
    # check that the maximum_number of clusters fixed is not greater than the number of row of the datafile
    checkMaxCluster <- function() {
        if (MAX_CLUSTERS > (nrow(data) - 1)) {
            message(paste("[WARNING] Max number of clusters must be lower (and not equal) to the number of line of the dataset (", nrow(data), ")", sep = ""))
            return(F)
        } else {
            return(T)
        }
    }

    # f: filename
    # func: plot function
    savePlot <- function(f, func) {
        f <- paste(f, ".tiff", sep = "")
        tiff(f, res = 300, width = 700, height = 800)
        func
        suprLog <- dev.off()
    }

    savePng <- function(f, func) {
        assign("PNG", T, .GlobalEnv)
        f <- paste(f, ".png", sep = "")
        png(f, DIM_PNG, DIM_PNG)
        func
        suprLog <- dev.off()
        assign("PNG", F, .GlobalEnv)
    }

    setData <- reactive({
        setClassifPar()
        tryCatch(
            {
                if (input$infile$size > 3000000) {
                    assign(
                        "VERBOSE_NIV2",
                        T,
                        .GlobalEnv
                    )
                } else {
                    assign(
                        "VERBOSE_NIV2",
                        F,
                        .GlobalEnv
                    )
                }

                refresh$data <- loadData(input$infile$datapath, input$sep, input$header)

                if (input$scale) {
                    refresh$data <- scale(refresh$data)
                } else {
                    refresh$data <- scale(refresh$data, scale = F)
                }

                if (input$transpose) {
                    refresh$data <- t(refresh$data)
                }

                if (VERBOSE_NIV2) {
                    cat("done.\n")
                }

                if (isSymmetric(as.matrix(refresh$data)) & !input$header) {
                    colnames(refresh$data) <- rownames(refresh$data)
                }
            },
            error = function(e) {
                message("[WARNING] Wrong separator, please selects another one.")
            }
        )
        setDistance()
    })

    setDistance <- reactive({
        printProgress(VERBOSE_NIV2, "Distance calculation")
        refresh$dis <- getDistance(refresh$data, as.integer(input$dist_type))
        if (VERBOSE_NIV2) cat("done.\n")

        setClassif()
    })

    setClassifPar <- reactive({
        refresh$max <- input$max_clusters
        refresh$classif_type <- as.integer(getClassifValue(input$classif_type))
    })

    ###################################
    #          EVENTS
    ###################################

    observeEvent(c(input$infile, input$header, input$sep, input$scale, input$transpose), {
        assign(
            "HEAD",
            input$header,
            .GlobalEnv
        )
        if (!is.null(input$infile)) {
            setData()
        }
    })

    observeEvent(input$dist_type, {
        if (!is.null(input$infile)) {
            setDistance()
        }
    })

    observeEvent(c(input$max_clusters, input$classif_type), {
        if (!is.null(input$infile)) {
            setClassifPar()
            setClassif()
        }
    })

    observeEvent(input$nb_clusters, {
        if (!is.null(input$infile) & input$nb_clusters > 1) {
            setClassif()
        }
    })

    observeEvent(c(input$max_biomark), {
        if (!is.null(input$infile)) {
            setPrintFuncs()
        }
    })

    # events for advanced mode
    observeEvent(input$advanced, {
        if (!is.null(input$infile)) {
            if (isTRUE(input$advanced)) {
                if (refresh$classif_type > 2) {
                    cat(paste("\nAGGLOMERATIVE COEFFICIENT: ", round(getCoefAggl(classif), 3), "\n", sep = ""))
                }

                if (nrow(data) < (NB_ROW_MAX / 2)) {
                    printProgress(VERBOSE_NIV2, "Gap statistics calculation")

                    assign(
                        "gap",
                        getGapPerPart(MAX_CLUSTERS, data, classif, NB_BOOTSTRAP),
                        .GlobalEnv
                    )
                    if (VERBOSE_NIV2) cat("done.\n")
                } else {
                    assign(
                        "gap",
                        NULL,
                        .GlobalEnv
                    )
                    message("\n[WARNING] Dataset too big to calculate a gap statistics.")
                }
                setPrintFuncs()
            } else {
                assign(
                    "gap",
                    NULL,
                    .GlobalEnv
                )
            }
        }
    })

    # hide either input options or tabs
    observe({

        # set an id in tabsetPanel (here "navbar") and for each tabs

        # default behaviour
        show(selector = "#navbar li a[data-value=dendr]")
        hide(selector = "#navbar li a[data-value=elbow]")
        hide(selector = "#navbar li a[data-value=gap]")
        hide(selector = "#navbar li a[data-value=gap2]")
        hide(selector = "#navbar li a[data-value=coph]")
        hide(selector = "#navbar li a[data-value=fusion]")
        hide(selector = "#navbar li a[data-value=within]")
        toggle(condition = input$advanced, id = "nb_clusters")
        toggle(condition = input$advanced, id = "axis1")
        toggle(condition = input$advanced, id = "axis2")

        # if(!is.null(input$infile)){ #catch condition when no data are loaded and adv is selected

        # responsive for a given condition
        toggle(condition = input$advanced, selector = "#navbar li a[data-value=elbow]")
        toggle(condition = input$advanced, selector = "#navbar li a[data-value=gap]")
        toggle(condition = input$advanced, selector = "#navbar li a[data-value=gap2]")
        toggle(condition = input$advanced, selector = "#navbar li a[data-value=within]")
        toggle(condition = (input$advanced & as.integer(getClassifValue(input$classif_type) > 2)), selector = "#navbar li a[data-value=fusion]")
        toggle(condition = (input$advanced & as.integer(getClassifValue(input$classif_type) > 2)), selector = "#navbar li a[data-value=coph]")
        toggle(condition = (as.integer(getClassifValue(input$classif_type) > 2)), selector = "#navbar li a[data-value=dendr]")

        # }
    })

    # function save_all
    observeEvent(input$save_all, {
        if (is.data.frame(data)) {
            setVariables()
            setPrintFuncs()
            writeTsv("summary", "summary.tsv", v = F)

            savePlot("best_clustering", plotBest())
            savePlot("silhouette", plotSil())
            savePlot("pca", plotPCA())
            savePlot("heatmap", plotHeatmap())

            if (as.integer(getClassifValue(input$classif_type) > 2)) {
                savePlot("dendrogram", plotDend())
            }

            if (isTRUE(input$advanced)) {
                if (as.integer(getClassifValue(input$classif_type) > 2)) {
                    savePlot("cohenetic", plotCoph())
                    savePlot("fusion", plotFus())
                }
                savePlot("elbow", plotElb())
                # if (nrow(data) < (NB_ROW_MAX/2)){
                savePlot("gap", plotGap())
                # }
                writeTsv("ctr_clus", "ctr_clus.tsv", v = F)
                writeTsv("ctr_part", "ctr_part.tsv", v = F)
                writeTsv("within_k", "within_k.tsv", v = F)
            }
        }
    })


    ###################################
    #          RENDER
    ###################################

    output$summary <- renderTable({
        tryCatch(
            {
                setVariables()
                if (checkMaxCluster()) {
                    observeEvent(input$summary_save, writeTsv("summary", "summary.tsv", v = F))
                    summary
                }
            },
            error = function(e) {
            }
        )
    })

    output$best_cluster <- renderPlot({
        tryCatch(
            {
                setVariables()
                if (checkMaxCluster()) {
                    observeEvent(input$best_save, savePlot("best_clustering", plotBest()))
                    plotBest()
                }
            },
            error = function(e) {
            }
        )
    })

    output$silhouette <- renderPlot({
        tryCatch(
            {
                setVariables()
                if (checkMaxCluster()) {
                    observeEvent(input$sil_save, savePlot("silhouette", plotSil()))
                    plotSil()
                }
            },
            error = function(e) {
            }
        )
    })

    output$pca <- renderPlot({
        tryCatch(
            {
                setVariables()
                if (checkMaxCluster()) {
                    assign(
                        "AXIS1",
                        input$axis1,
                        .GlobalEnv
                    )
                    assign(
                        "AXIS2",
                        input$axis2,
                        .GlobalEnv
                    )
                    observeEvent(input$pca_save, {
                        if (nrow(as.matrix(data)) > NB_ROW_MAX) {
                            savePng("pca", plotPCA())
                        } else {
                            savePlot("pca", plotPCA())
                        }
                    })
                    plotPCA()
                }
            },
            error = function(e) {
            }
        )
    })

    output$heatmap <- renderPlot({
        tryCatch(
            {
                setVariables()
                if (checkMaxCluster()) {
                    par(mar = c(0, 0, 0, 0))
                    plot(0:1, 0:1, axes = F, type = "n") # delete sil plot
                    observeEvent(input$heatmap_save, {
                        if (nrow(as.matrix(data)) > NB_ROW_MAX) {
                            savePng("heatmap", plotHeatmap())
                        } else {
                            savePlot("heatmap", plotHeatmap())
                        }
                    })
                    plotHeatmap()
                }
            },
            error = function(e) {
            }
        )
    })

    output$cophenetic <- renderPlot({
        tryCatch(
            {
                setVariables()
                if (isTRUE(input$advanced) & CLASSIF_TYPE > 2) {
                    if (checkMaxCluster()) {
                        observeEvent(input$coph_save, {
                            if (nrow(as.matrix(data)) > NB_ROW_MAX) {
                                savePng("cohenetic", plotCoph())
                            } else {
                                savePlot("cohenetic", plotCoph())
                            }
                        })
                        plotCoph()
                    }
                }
            },
            error = function(e) {
            }
        )
    })

    output$dendrogram <- renderPlot({
        tryCatch(
            {
                setVariables()
                if (CLASSIF_TYPE > 2) {
                    if (checkMaxCluster()) {
                        observeEvent(input$dendr_save, savePlot("dendrogram", plotDend()))
                        plotDend()
                    }
                }
            },
            error = function(e) {
            }
        )
    })

    output$fusion <- renderPlot({
        tryCatch(
            {
                setVariables()
                if (isTRUE(input$advanced) & CLASSIF_TYPE > 2) {
                    if (checkMaxCluster()) {
                        observeEvent(input$fusion_save, savePlot("fusion", plotFus()))
                        plotFus()
                    }
                }
            },
            error = function(e) {
            }
        )
    })

    output$gap <- renderPlot({
        tryCatch(
            {
                setVariables()
                if (isTRUE(input$advanced)) {
                    if (checkMaxCluster()) {
                        observeEvent(input$gap_save, savePlot("gap", plotGap()))
                        plotGap()
                    }
                }
            },
            error = function(e) {
            }
        )
    })

    output$gap2 <- renderPlot({
        tryCatch(
            {
                setVariables()
                if (isTRUE(input$advanced)) {
                    if (checkMaxCluster()) {
                        observeEvent(input$gap2_save, savePlot("gap", plotGap2()))
                        plotGap2()
                    }
                }
            },
            error = function(e) {
            }
        )
    })

    output$elbow <- renderPlot({
        tryCatch(
            {
                setVariables()
                if (isTRUE(input$advanced)) {
                    if (checkMaxCluster()) {
                        observeEvent(input$elbow_save, savePlot("elbow", plotElb()))
                        plotElb()
                    }
                }
            },
            error = function(e) {
            }
        )
    })

    output$within <- renderTable(
        {
            tryCatch(
                {
                    setVariables()
                    if (isTRUE(input$advanced)) {
                        if (checkMaxCluster()) {
                            observeEvent(input$within_save, writeTsv("within_k", "within_k.tsv", v = F))
                            within_k
                        }
                    }
                },
                error = function(e) {
                }
            )
        },
        rownames = T,
        hover = T,
        striped = T,
        digits = 2,
        width = "100cm",
        align = "c",
        na = "",
        size = 200
    )

    output$ctr_clus_plot <- renderPlot({
        tryCatch(
            {
                setVariables()
                if (checkMaxCluster()) {
                    observeEvent(input$ctr_clus_plot_save, savePlot("discr_var", ctr_clus_plot()))
                    ctr_clus_plot()
                }
            },
            error = function(e) {
            }
        )
    })

    output$ctr_clus <- renderTable({
        tryCatch(
            {
                setVariables()
                if (checkMaxCluster()) {
                    observeEvent(input$ctr_clus_save, writeTsv("ctr_clus", "ctr_clus.tsv", v = F))
                    ctr_clus
                }
            },
            error = function(e) {
            }
        )
    })

    output$ctr_part <- renderTable({
        tryCatch(
            {
                setVariables()
                if (checkMaxCluster()) {
                    observeEvent(input$ctr_part_save, writeTsv("ctr_part", "ctr_part.tsv", v = F))
                    ctr_part
                }
            },
            error = function(e) {
            }
        )
    })

    output$centroids <- renderTable({
        tryCatch(
            {
                setVariables()
                if (checkMaxCluster()) {
                    observeEvent(input$centroids_save, writeTsv("centroids_save", "centroids.tsv", v = F))
                    aggregate(data, list(clusters), mean)
                }
            },
            error = function(e) {
            }
        )
    })
}
