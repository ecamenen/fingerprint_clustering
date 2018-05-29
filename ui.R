library(shiny)

classif_methods = list("K-menoids" = 1,  "K-means" = 2, "Ward"=3, "Complete links"=4, "Single links"=5, "UPGMA"=6, "WPGMA"=7, "WPGMC"=8, "UPGMC"=9)

#min-max: minimum and maximum position number in the list
getClassifKeys = function(min, max)
  #similar to for i in min:max
  unlist(sapply(min:max, function(i) names(classif_methods)[i]))

shinyUI( pageWithSidebar(
  headerPanel("Fingerprint clustering"), #titre appl
  
  #pannel de cote
  sidebarPanel( #elements de controle de l'apli
    
    textInput(inputId="infile", #nom (ou ID) associee a cet element de controle,
                 #sera utilisee dans la partie 'server'
                 label=h5("Name of the dataset: "), #libellé associee a cet element 
                 #(apparait au dessus de l'input)
                 value="matrix.txt"),
    
    #checkValues and checkNames do not works on selectInput but only on checkbox
    selectInput("classif_type",h5("Classification method: "),  selected = "Complete links",
                 choices = list(`Partitionning clustering` =  getClassifKeys(1,2), `Hierarchical clustering` = getClassifKeys(3,9))),
    
    sliderInput("max_clusters", h5("Maximum number of clusters: "), min=2, max= 100, value=6),
    sliderInput("nb_clusters", h5("Number of clusters: "), min=0, max= 10, value=0),
    checkboxInput("advanced", "Advanced mode", value=F)
  ),
  
  mainPanel(
    tabsetPanel(
      type = "tabs",
      tabPanel("Summary",tableOutput("summary")),
      tabPanel("Best clustering", plotOutput("best_cluster")),
      tabPanel("Silhouette", plotOutput("silhouette")),
      tabPanel("PCA", plotOutput("pca")),
      tabPanel("Heatmap", plotOutput("heatmap")),
      tabPanel("Cophenetic", plotOutput("cophenetic")),
      tabPanel("Dendrogram", plotOutput("dendrogram")),
      tabPanel("Contribution per cluster",tableOutput("ctr_clus")),
      tabPanel("Contribution per partition",tableOutput("ctr_part"))
    )
  )
  
))
