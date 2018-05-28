library(shiny)
classif_methods = list("K-menoids" = 1,  "K-means" = 2, "Ward"=3, "Complete links"=4, "Single links"=5, "UPGMA"=6, "WPGMA"=7, "WPGMC"=8, "UPGMC"=9)

#min-max: minimum and maximum position number in the list
getClassifKeys = function(min, max)
  #similar to for i in min:max
  unlist(sapply(min:max, function(i) names(classif_methods)[i]))


shinyUI( pageWithSidebar(
  #titlePanel("Hello"),
  headerPanel("Fingerprint clustering"), #titre appl
  
  #pannel de cote
  sidebarPanel( #elements de controle de l'apli
    
    textInput(inputId="infile", #nom (ou ID) associee a cet element de controle,
                 #sera utilisee dans la partie 'server'
                 label=h4("Name of the dataset: "), #libellé associee a cet element 
                 #(apparait au dessus de l'input)
                 value="matrix.txt"),
    
    #checkValues and checkNames do not works on selectInput but only on checkbox
    selectInput(inputId="classif_type",
                 label=h4("Méthode de classification: "),
                 selected = "Complete links",
                 choices = list(`Partitionning clustering` = c(names(classif_methods)[1], names(classif_methods)[2]), 
                                `Hierarchical clustering` = getClassifKeys(3,9)))

    #names(classif_methods)[2]
    #unlist(classif_methods[2])
    
  ),
  
  mainPanel(
    #titre donne a la partie presentant les sorties
    #h3("Summary: "),
    
    #va permettre d'afficher une sortie de type 'texte'
    #le contenu sera defini dans la partie server avec comme nom de variable
    #'sortie de texte'
    tableOutput("summary"),
    plotOutput("fusion_levels")
    #en sortie, tu attend du texte
  )
  
))
