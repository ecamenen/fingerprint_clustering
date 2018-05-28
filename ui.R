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
                 label=h4("Name of the dataset: "), #libellé associee a cet element 
                 #(apparait au dessus de l'input)
                 value="matrix.txt"),
    
    #checkValues and checkNames do not works on selectInput but only on checkbox
    selectInput("classif_type",h4("Classification method: "),  selected = "Complete links",
                 choices = list(`Partitionning clustering` =  getClassifKeys(1,2), `Hierarchical clustering` = getClassifKeys(3,9))),
    
    sliderInput("max_clusters", h4("Maximum number of clusters allowed: "), min=2, max= 100, value=6),
    sliderInput("nb_clusters", h4("Number of clusters: "), min=2, max= 100, value=2)
  ),
  
  mainPanel(
    #titre donne a la partie presentant les sorties
    #h3("Summary: ")
    #le contenu sera defini dans la partie server
    tableOutput("summary"),
    plotOutput("fusion_levels"),
    plotOutput("silhouette"),
    plotOutput("pca"),
    plotOutput("heatmap"),
    plotOutput("cophenetic"),
    plotOutput("dendrogram")
  )
  
))
