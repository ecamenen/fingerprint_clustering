library("sys")
library("shinyjs")



classif_methods = list("K-menoids" = 1,  "K-means" = 2, "Ward"=3, "Complete links"=4, "Single links"=5, "UPGMA"=6, "WPGMA"=7, "WPGMC"=8, "UPGMC"=9)

#min-max: minimum and maximum position number in the list
getClassifKeys = function(min, max)
  #similar to for i in min:max
  unlist(sapply(min:max, function(i) names(classif_methods)[i]))

ui = fluidPage(
  useShinyjs(),
  
  includeCSS("www/animate.min.css"),
  # provides pulsating effect
  includeCSS("www/loading-content.css"),

  pageWithSidebar(
    headerPanel("Fingerprint clustering"), #titre appl
    
    #pannel de cote
    sidebarPanel( #elements de controle de l'apli
      
      fileInput(inputId="infile", #nom (ou ID) associee a cet element de controle,
                   #sera utilisee dans la partie 'server'
                   label=h5("Choose a file that containt MetExplore metabolite ids: ") #libellé associee a cet element 
                   #(apparait au dessus de l'input)
                ),
      selectInput("algoShortestPath",
                  h5("Shortest path algorithm : "),
                  selected = "Biological shortest path",
                  choices = list(`Biological shortest path` =  "ValidShortest",
                                 `Proximity shortest path` = "ShortestAsUndirected")),
      
      actionButton("conputeShortestPath",
                   "Conpute Shortest Path"),
      
      hr(style="height: 2px; color: #0e325e; background-color: #0e325e; width: 75%; border: none;"),
      #checkValues and checkNames do not works on selectInput but only on checkbox
      selectInput("classif_type",
                  h5("Classification method: "),
                  selected = "Complete links",
                  choices = list(`Partitionning clustering` =  getClassifKeys(1,2),
                                 `Hierarchical clustering` = getClassifKeys(3,9))),
      sliderInput("max_clusters", 
                  h5("Maximum number of clusters: "), 
                  min=2, max= 100, value=6),
      sliderInput("nb_clusters",
                  h5("Number of clusters: "),
                  min=0, max= 10, value=0),
      checkboxInput("advanced",
                    "Advanced mode",
                    value=T),
      actionButton("save_all",
                   "Save all")
    ),
    
    mainPanel(
      div(
        id = "no-content",
        class = "no-content",
        h2("Select fingerprint file...")
      ),
      div(
        id = "loading-content",
        class = "loading-content",
        h2(class = "animated infinite pulse", "Loading data...")
      ),
      tabsetPanel(
        id = "tabsetPanel",
        type = "tabs",
        
        tabPanel("Summary",
                 tableOutput("summary"),
                 actionButton("summary_save","Save")),
        tabPanel("Best clustering", 
                 plotOutput("best_cluster"),
                 actionButton("best_save","Save")),
        tabPanel("Silhouette",
                 plotOutput("silhouette"),
                 actionButton("sil_save","Save")),
        tabPanel("PCA",
                 plotOutput("pca"),
                 actionButton("pca_save","Save")),
        tabPanel("Heatmap",
                 plotOutput("heatmap"),
                 actionButton("heatmap_save","Save")),
        tabPanel("Cophenetic",
                 plotOutput("cophenetic"),
                 actionButton("coph_save","Save")),
        tabPanel("Dendrogram",
                 plotOutput("dendrogram"),
                 actionButton("dendr_save","Save")),
        tabPanel("Partition contributions",
                 tableOutput("ctr_part"),
                 actionButton("ctr_part_save","Save")),
        tabPanel("Cluster contributions",
                 tableOutput("ctr_clus"),
                 actionButton("ctr_clus_save","Save"))
      )
    )
  )
  
)
