#' The application User-Interface
#'
#' @param request Internal parameter for `{shiny}`.
#' @noRd
app_ui <- function(request) {
    classif_methods <- list("K-menoids" = 1, "K-means" = 2, "Ward" = 3, "Complete links" = 4, "Single links" = 5, "UPGMA" = 6, "WPGMA" = 7, "WPGMC" = 8, "UPGMC" = 9)
    jscode <- "shinyjs.refresh = function() { location.reload(); }"

    # min-max: minimum and maximum position number in the list
    getClassifKeys <- function(min, max) {
        # similar to for i in min:max
        unlist(sapply(min:max, function(i) names(classif_methods)[i]))
    }

    tagList(
        golem_add_external_resources(),
        fluidPage(
            title = "Automatic Clustering",
            h1(toupper("Automatic Clustering")),
            tags$p(
                tags$a(
                    href = "https://github.com/ecamenen/tcgaViz",
                    "Etienne CAMENEN,"
                ),
                "Maxime CHAZALVIEL, Fabien JOURDAN"
            ),
            h4(
                paste0(
                    "Performs unsupervised clustering and automatically",
                    "determine the best number of cluster."
                )
            ),
            sidebarLayout(
                sidebarPanel(

                    useShinyjs(),
                    extendShinyjs(text = jscode, functions = "refresh"),
                    fileInput(
                        inputId = "infile",
                        label = h5("Choose a file: ")
                    ),
                    checkboxInput("header",
                        "Consider first row as header",
                        value = TRUE
                    ),
                    radioButtons("sep",
                        "Separator",
                        choices = c(
                            Comma = ",",
                            Semicolon = ";",
                            Tab = "\t"
                        ),
                        selected = "\t"
                    ),
                    checkboxInput("scale",
                        "Scale",
                        value = TRUE
                    ),
                    checkboxInput("transpose",
                        "Transpose",
                        value = FALSE
                    ),
                    selectInput("dist_type",
                        h5("Distance method: "),
                        selected = 1,
                        choices = list(
                            `Distance` = c("Euclidean" = 1, "Manhattan" = 2),
                            `Similarity` = c("Jaccard" = 3, "Sokal & Michener" = 4, "Sorensen (Dice)" = 5, "Ochiai" = 6)
                        )
                    ),
                    selectInput("classif_type",
                        h5("Classification method: "),
                        selected = "Ward",
                        choices = list(
                            `Partitionning clustering` = getClassifKeys(1, 2),
                            `Hierarchical clustering` = getClassifKeys(3, 9)
                        )
                    ),
                    sliderInput("max_clusters",
                        h5("Maximum number of clusters: "),
                        min = 2, max = 100, value = 6
                    ),
                    sliderInput("max_biomark",
                        h5("Maximum number of biomarkers: "),
                        min = 10, max = 100, value = 10
                    ),
                    checkboxInput("advanced",
                        "Advanced mode",
                        value = F
                    ),
                    sliderInput("nb_clusters",
                        h5("Number of clusters: "),
                        min = 0, max = 10, value = 0
                    ),
                    sliderInput("axis1",
                        h5("PCA axis 1: "),
                        min = 1, max = 4, value = 1
                    ),
                    sliderInput("axis2",
                        h5("PCA axis 2: "),
                        min = 2, max = 4, value = 2
                    ),
                    actionButton(
                        "save_all",
                        "Save all"
                    )
                ),
                mainPanel(
                    tabsetPanel(
                        type = "tabs",
                        id = "navbar",
                        tabPanel(
                            "Summary",
                            tableOutput("summary"),
                            actionButton("summary_save", "Save")
                        ),
                        tabPanel(
                            "Best clustering",
                            plotOutput("best_cluster"),
                            actionButton("best_save", "Save")
                        ),
                        tabPanel(
                            "Silhouette",
                            plotOutput("silhouette", height = 600),
                            actionButton("sil_save", "Save")
                        ),
                        tabPanel(
                            "PCA",
                            plotOutput("pca"),
                            actionButton("pca_save", "Save")
                        ),
                        tabPanel(
                            "Heatmap",
                            plotOutput("heatmap"),
                            actionButton("heatmap_save", "Save")
                        ),
                        tabPanel("Cophenetic",
                            value = "coph",
                            plotOutput("cophenetic"),
                            actionButton("coph_save", "Save")
                        ),
                        tabPanel("Dendrogram",
                            value = "dendr",
                            plotOutput("dendrogram"),
                            actionButton("dendr_save", "Save")
                        ),
                        tabPanel("Elbow",
                            value = "elbow",
                            plotOutput("elbow"),
                            actionButton("elbow_save", "Save")
                        ),
                        tabPanel("Gap",
                            value = "gap",
                            plotOutput("gap"),
                            actionButton("gap_save", "Save")
                        ),
                        tabPanel("Gap2",
                            value = "gap2",
                            plotOutput("gap2"),
                            actionButton("gap2_save", "Save")
                        ),
                        tabPanel("Fusion",
                            value = "fusion",
                            plotOutput("fusion"),
                            actionButton("fusion_save", "Save")
                        ),
                        tabPanel("Within-inertia",
                            value = "within",
                            tableOutput("within"),
                            actionButton("within_save", "Save")
                        ),
                        tabPanel(
                            "Discriminant variables",
                            plotOutput("ctr_clus_plot"),
                            actionButton("ctr_clus_plot_save", "Save")
                        ),
                        tabPanel(
                            "Discriminant variables",
                            tableOutput("ctr_clus"),
                            actionButton("ctr_clus_save", "Save")
                        ),
                        tabPanel(
                            "Discriminant power",
                            tableOutput("ctr_part"),
                            actionButton("ctr_part_save", "Save")
                        ),
                        tabPanel(
                            "Centroids",
                            tableOutput("centroids"),
                            actionButton("centroids_save", "Save")
                        )
                    )
                )
            )
        )
    )
}

#' Add external Resources to the Application
#'
#' This function is internally used to add external
#' resources inside the Shiny application.
#'
#' @importFrom golem add_resource_path activate_js favicon bundle_resources
#' @noRd
golem_add_external_resources <- function() {
    add_resource_path(
        "www", app_sys("app/www")
    )

    tags$head(
        favicon(),
        bundle_resources(
            path = app_sys("app/www"),
            app_title = "autoCluster"
        )
    )
}
