library(shiny)
shinyUI( pageWithSidebar(
  #titlePanel("Hello"),
  headerPanel("Fingerprint clustering"), #titre appl
  
  #pannel de cote
  sidebarPanel( #elements de controle de l'apli
    
    textInput(inputId="infile", #nom (ou ID) associee a cet element de controle,
                 #sera utilisee dans la partie 'server'
                 label="Name of the dataset: ", #libellé associee a cet element 
                 #(apparait au dessus de l'input)
                 value="matrix.txt")
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