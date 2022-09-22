#' The application server-side
#'
#' @param input,output,session Internal parameters for {shiny}.
#' @noRd
app_server <- function(input, output, session) {
    options(shiny.maxRequestSize = 30 * 1024^2)
    # Global variables settings
    NB_BOOTSTRAP <- 500 # should be comprise between 100 and 1000
    TEXT <- TRUE # print values on graph (for optimum partition and heatmap)
    NB_ROW_MAX <- 200 # max row to have pdf, otherwise, some plots are in png
    DIM_PNG <- 2000
    MAX_CHAR_LEN <- 25 # maximum length of individual s names
    classif_methods <- list(
        "K-menoids" = 1,
        "K-means" = 2,
        "Ward" = 3,
        "Complete links" = 4,
        "Single links" = 5,
        "UPGMA" = 6,
        "WPGMA" = 7,
        "WPGMC" = 8,
        "UPGMC" = 9
    )
    vars <- reactiveValues(
        data = NULL,
        dis = NULL,
        max_clusters = NULL,
        classif_type = NULL,
        nb_clusters = NULL,
        advanced = NULL,
        axis1 = NULL,
        axis2 = NULL,
        classif = NULL,
        clusters = NULL,
        cl_k = NULL,
        between = NULL,
        sils = NULL,
        mean_sils = NULL,
        optimal_k = NULL,
        gap = NULL,
        sil_k = NULL,
        advanced = NULL,
        png = FALSE,
        verbose2 = FALSE,
        pca = FALSE
    )
    tryCatch(
        {
            vars$data <- loadData("inst/extdata/matrix.txt")
        },
        warning = function(w) {
            vars$data <- NULL
            message("Default file \"matrix.txt\" is not in the folder. Please, load another one.")
        },
        error = function(e) {
            vars$data <- NULL
        }
    )

    loadData <- function(f, s = "\t", h = FALSE) {
        # !file.exists(
        if (!is.null(f)) {
            if (grepl("xlsx?", f)) {
                data <- openxlsx::read.xlsx(f, rowNames = TRUE)
            } else {
                data <- read.table(
                    f,
                    header = h,
                    sep = s,
                    dec = ".",
                    row.names = 1
                )
            }
            # vars$data = preProcessData(vars$data, vars$header)
            # substr(rownames(vars$data), 1, MAX_CHAR_LEN) -> rownames(vars$data)
        }
        return(data)
    }

    refresh <- reactiveValues()
    refresh$classif_type <- ""
    getClassifValue <- function(key) {
        unlist(classif_methods[key])
    }

    ###################################
    #          SETTINGS
    ###################################

    setVariables <- reactive({
        refresh2 <- c(refresh$classif_type, input$max_biomark)

        vars$data <- refresh$data
        vars$dis<- refresh$dis
        vars$max_clusters <- refresh$max
        vars$classif_type <- refresh$classif_type
        vars$nb_clusters <- input$nb_clusters
        vars$advanced <- input$advanced
        vars$axis1 <- input$axis1
        vars$axis2 <- input$axis2
    })

    # Each shiny server functions run in local environment
    # With assign, variables are forced to be in global env
    setClassif <- reactive({
        setVariables()
        req(input$infile)

        if ((nrow(vars$data) > 3000) &
            (input$classif_type > 2)) {
            message(
                "[WARNING] With more than 3000 rows to analyse, classification method must be K-medoids or K-means",
                call. = FALSE
            )
        }
        # Perform classification
        printProgress(vars$verbose2, "Classification")

        vars$classif <- getClassif(vars$classif_type, vars$max_clusters, vars$data, vars$dis)
        if (vars$verbose2) {
            cat("done.\n")
        }
        vars$clusters <- getClusterPerPart(vars$max_clusters + 1, vars$classif)

        # inertia
        vars$between <- getRelativeBetweenPerPart(vars$max_clusters, vars$data, vars$clusters)
        vars$diff_between <- getBetweenDifferences(vars$between)

        printProgress(vars$verbose2, "PCA")
        vars$pca <- dudi.pca(
            vars$data,
            scannf = FALSE,
            scale = input$scale,
            nf = 4
        )
        if (vars$verbose2) {
            cat("done.\n")
        }

        vars$sils <- getSilhouettePerPart(vars$data, vars$clusters, vars$dis)
        vars$mean_sils <- getMeanSilhouettePerPart(vars$sils)

        setClusters()

        setPrintFuncs()
    })

    setClusters <- reactive({
        refr <- c(
            refresh$classif_type,
            input$advanced,
            input$nb_clusters,
            input$scale,
            input$transpose
        )

        if (vars$nb_clusters > 0) {
          vars$optimal_k <- input$nb_clusters
        } else {
          vars$optimal_k <- which.max(vars$mean_sils) + 1
        }

        # cl_temp, because writeTsv(clusters) recreate a different object named clusters
        vars$cl_k <- vars$clusters[[vars$optimal_k - 1]]

        if (!input$advanced) {
          vars$gap <- NULL
        }

        vars$sil_k <- vars$sils[[vars$optimal_k - 1]]

        # errors
        if (vars$optimal_k == vars$max_clusters) {
            message(
                "\n[WARNING] The optimal number of clusters equals the maximum number of clusters. \nNo cluster structure has been found."
            )
        }
        if (min(table(vars$cl_k)) == 1) {
            message("\n[WARNING] A cluster with an only singleton biased the silhouette score.")
        }

        writeClusters("clusters.tsv", vars$sil_k, vars$pca, v = FALSE)
    })

    setPrintFuncs <- function() {
        vars$advanced <- input$advanced
        refr <- c(input$nb_clusters, input$max_biomark)

        ###### plot funcs #####
        assign(
            "plotPCA",
            function() {
                plotPca(vars$pca, vars$data, vars$cl_k, vars$axis1, vars$axis2, vars$advanced, vars$png)
            },
            .GlobalEnv
        )
        assign(
            "plotBest",
            function() {
                plotSilhouettePerPart(vars$mean_sils, vars$sils)
            },
            .GlobalEnv
        )
        assign(
            "plotSil",
            function() {
                plotSilhouette(vars$sil_k)
            },
            .GlobalEnv
        )
        assign(
            "plotDend",
            function() {
                plotDendrogram(
                    vars$classif_type,
                    vars$optimal_k,
                    vars$classif,
                    vars$data,
                    vars$max_clusters,
                    vars$cl_k
                )
            },
            .GlobalEnv
        )

        if (vars$classif_type <= 2 | isTRUE(vars$advanced)) {
            assign(
                "plotHeatmap",
                function() {
                    heatMap(vars$data, vars$dis, vars$sil_k, is_png = vars$png, verbose = vars$verbose2)
                },
                .GlobalEnv
            )
        } else {
            assign(
                "plotHeatmap",
                function() {
                    heatMap(vars$data, vars$dis, c = vars$classif, cl = vars$cl_k, is_png = vars$png, verbose = vars$verbose2)
                },
                .GlobalEnv
            )
        }

        ##### advanced #####
        if (isTRUE(vars$advanced)) {
            assign(
                "plotFus",
                function() {
                    plotFusionLevels(vars$max_clusters, vars$classif)
                },
                .GlobalEnv
            )
            assign(
                "plotCoph",
                function() {
                    printProgress(vars$verbose2, "Cophenetic calculation")
                    plotCohenetic(vars$dis, vars$classif, vars$png)
                    if (vars$verbose2) {
                        cat("done.\n")
                    }
                },
                .GlobalEnv
            )
            assign(
                "plotGap",
                function() {
                    if (nrow(vars$data) < (NB_ROW_MAX / 2)) {
                        plotGapPerPart(vars$gap, vars$max_clusters, v = FALSE)
                        # plotGapPerPart2(gap, vars$max_clusters)
                    } else {
                        message("\n[WARNING] Dataset too big to calculate a gap statistics.")
                    }
                },
                .GlobalEnv
            )
            assign(
                "plotGap2",
                function() {
                    if (nrow(vars$data) < (NB_ROW_MAX / 2)) {
                        plotGapPerPart2(vars$gap, vars$max_clusters)
                    } else {
                        message("\n[WARNING] Dataset too big to calculate a gap statistics.")
                    }
                },
                .GlobalEnv
            )
            assign(
                "plotElb",
                function() {
                    plotElbow(vars$between)
                },
                .GlobalEnv
            )
            assign(
                "within_k",
                getRelativeWithinPerCluster(vars$clusters, vars$data),
                .GlobalEnv
            )
        }

        ##### print table func #####

        assign(
            "summary_table",
            printSummary(vars$between, vars$diff_between, vars$mean_sils, vars$advanced, vars$gap),
            .GlobalEnv
        )
        assign(
            "ctr_part",
            100 * getPdisPerPartition(vars$classif_type, vars$max_clusters, vars$clusters, vars$data),
            .GlobalEnv
        )
        assign(
            "centroids",
            getDistPerVariable(vars$data, vars$cl_k),
            .GlobalEnv
        )
        assign(
            "discr",
            getDiscriminantVariables(
                vars$classif_type,
                vars$optimal_k,
                vars$cl_k,
                vars$data,
                input$max_biomark
            ),
            .GlobalEnv
        )
        assign(
            "ctr_clus_plot",
            function() {
                plotDiscriminantVariables(discr)
            },
            .GlobalEnv
        )
        assign(
            "ctr_clus",
            100 * getCtrVar(vars$classif_type, vars$optimal_k, vars$cl_k, vars$data),
            .GlobalEnv
        )
        writeTsv(discr, "discr_var.tsv", v = FALSE)
    }

    # post-process for data
    # check that the maximum_number of clusters fixed is not greater than the number of row of the datafile
    checkMaxCluster <- function() {
        if (vars$max_clusters > (nrow(vars$data) - 1)) {
            message(
                paste(
                    "[WARNING] Max number of clusters must be lower (and not equal) to the number of line of the dataset (",
                    nrow(vars$data),
                    ")",
                    sep = ""
                )
            )
            return(FALSE)
        } else {
            return(TRUE)
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
        vars$png <- TRUE
        f <- paste(f, ".png", sep = "")
        png(f, DIM_PNG, DIM_PNG)
        func
        suprLog <- dev.off()
        vars$png <- FALSE
    }

    setData <- reactive({
        setClassifPar()
        tryCatch(
            {
                if (input$infile$size > 3000000) {
                    vars$verbose2 <- TRUE
                } else {
                    vars$verbose2 <- FALSE
                }

                refresh$data <- loadData(input$infile$datapath, input$sep, input$header)

                if (input$scale) {
                    refresh$data <- scale(refresh$data)
                } else {
                    refresh$data <- scale(refresh$data, scale = FALSE)
                }

                if (input$transpose) {
                    refresh$data <- t(refresh$data)
                }

                if (vars$verbose2) {
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
        printProgress(vars$verbose2, "Distance calculation")
        refresh$dis <- getDistance(refresh$data, as.integer(input$dist_type))
        if (vars$verbose2) {
            cat("done.\n")
        }

        setClassif()
    })

    setClassifPar <- reactive({
        refresh$max <- input$max_clusters
        refresh$classif_type <- as.integer(getClassifValue(input$classif_type))
    })

    ###################################
    #          EVENTS
    ###################################

    observeEvent(c(
        input$infile,
        input$header,
        input$sep,
        input$scale,
        input$transpose
    ), {
        vars$head <- input$header
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
    observeEvent(c(input$advanced, input$max_clusters), {
        if (!is.null(input$infile)) {
            if (isTRUE(input$advanced)) {
                if (refresh$classif_type > 2) {
                    cat(paste(
                        "\nAGGLOMERATIVE COEFFICIENT: ",
                        round(getCoefAggl(vars$classif), 3),
                        "\n",
                        sep = ""
                    ))
                }

                if (nrow(vars$data) < (NB_ROW_MAX / 2)) {
                    printProgress(vars$verbose2, "Gap statistics calculation")

                    vars$gap <- getGapPerPart(vars$max_clusters, vars$data, vars$classif, NB_BOOTSTRAP)
                    if (vars$verbose2) {
                        cat("done.\n")
                    }
                } else {
                    vars$gap <- NULL
                    message("\n[WARNING] Dataset too big to calculate a gap statistics.")
                }
                setPrintFuncs()
            } else {
              vars$gap <- NULL
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
        toggle(
            condition = (input$advanced &
                as.integer(
                    getClassifValue(input$classif_type) > 2
                )),
            selector = "#navbar li a[data-value=fusion]"
        )
        toggle(
            condition = (input$advanced &
                as.integer(
                    getClassifValue(input$classif_type) > 2
                )),
            selector = "#navbar li a[data-value=coph]"
        )
        toggle(condition = (as.integer(
            getClassifValue(input$classif_type) > 2
        )), selector = "#navbar li a[data-value=dendr]")

        # }
    })

    # function save_all
    observeEvent(input$save_all, {
        if (is.data.frame(vars$data)) {
            setVariables()
            setPrintFuncs()
            writeTsv(summary_table, "summary.tsv", v = FALSE)

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
                # if (nrow(vars$data) < (NB_ROW_MAX/2)){
                savePlot("gap", plotGap())
                # }
                writeTsv(ctr_clus, "ctr_clus.tsv", v = FALSE)
                writeTsv(ctr_part, "ctr_part.tsv", v = FALSE)
                writeTsv(within_k, "within_k.tsv", v = FALSE)
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
                    observeEvent(
                        input$summary_save,
                        writeTsv("summary_table", "summary.tsv", v = FALSE)
                    )
                    summary_table
                }
            },
            error = function(e) {
            }
        )
    })

    output$best_cluster <- renderPlot({
        # tryCatch(
        # {
        setVariables()
        if (checkMaxCluster()) {
            observeEvent(
                input$best_save,
                savePlot("best_clustering", plotBest())
            )
            plotBest()
        }
        #     },
        #     error = function(e) {
        #     }
        # )
    })

    output$silhouette <- renderPlot({
        # tryCatch(
        #     {
        setVariables()
        if (checkMaxCluster()) {
            observeEvent(input$sil_save, savePlot("silhouette", plotSil()))
            plotSil()
        }
        #     },
        #     error = function(e) {
        #     }
        # )
    })

    output$pca <- renderPlot({
        # tryCatch(
        #     {
        setVariables()
        if (checkMaxCluster()) {
          vars$axis1 <- input$axis1
          vars$axis2 <- input$axis2
            observeEvent(input$pca_save, {
                if (nrow(as.matrix(vars$data)) > NB_ROW_MAX) {
                    savePng("pca", plotPCA())
                } else {
                    savePlot("pca", plotPCA())
                }
            })
            plotPCA()
        }
        #     },
        #     error = function(e) {
        #     }
        # )
    })

    output$heatmap <- renderPlot({
        # tryCatch(
        #     {
        setVariables()
        if (checkMaxCluster()) {
            par(mar = c(0, 0, 0, 0))
            plot(0:1, 0:1, axes = FALSE, type = "n") # delete sil plot
            observeEvent(input$heatmap_save, {
                if (nrow(as.matrix(vars$data)) > NB_ROW_MAX) {
                    savePng("heatmap", plotHeatmap())
                } else {
                    savePlot("heatmap", plotHeatmap())
                }
            })
            plotHeatmap()
        }
        #     },
        #     error = function(e) {
        #     }
        # )
    })

    output$cophenetic <- renderPlot({
        # tryCatch(
        #     {
        setVariables()
        if (isTRUE(input$advanced) & vars$classif_type > 2) {
            if (checkMaxCluster()) {
                observeEvent(input$coph_save, {
                    if (nrow(as.matrix(vars$data)) > NB_ROW_MAX) {
                        savePng("cohenetic", plotCoph())
                    } else {
                        savePlot("cohenetic", plotCoph())
                    }
                })
                plotCoph()
            }
        }
        #     },
        #     error = function(e) {
        #     }
        # )
    })

    output$dendrogram <- renderPlot({
        # tryCatch(
        #     {
        setVariables()
        if (vars$classif_type > 2) {
            if (checkMaxCluster()) {
                observeEvent(
                    input$dendr_save,
                    savePlot("dendrogram", plotDend())
                )
                plotDend()
            }
        }
        #     },
        #     error = function(e) {
        #     }
        # )
    })

    output$fusion <- renderPlot({
        tryCatch(
            {
                setVariables()
                if (isTRUE(input$advanced) & vars$classif_type > 2) {
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
                            observeEvent(
                                input$within_save,
                                writeTsv("within_k", "within_k.tsv", v = FALSE)
                            )
                            within_k
                        }
                    }
                },
                error = function(e) {
                }
            )
        },
        rownames = TRUE,
        hover = TRUE,
        striped = TRUE,
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
                    observeEvent(
                        input$ctr_clus_plot_save,
                        savePlot("discr_var", ctr_clus_plot())
                    )
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
                    observeEvent(
                        input$ctr_clus_save,
                        writeTsv("ctr_clus", "ctr_clus.tsv", v = FALSE)
                    )
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
                    observeEvent(
                        input$ctr_part_save,
                        writeTsv("ctr_part", "ctr_part.tsv", v = FALSE)
                    )
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
                    observeEvent(
                        input$centroids_save,
                        writeTsv("centroids_save", "centroids.tsv", v = FALSE)
                    )
                    aggregate(vars$data, list(vars$cl_k), mean)
                }
            },
            error = function(e) {
            }
        )
    })
}
