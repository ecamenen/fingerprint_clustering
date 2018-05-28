library(shiny)
source("fingerprint_clustering.R")
classif_methods = list("K-menoids" = 1,  "K-means" = 2, "Ward"=3, "Complete links"=4, "Single links"=5, "UPGMA"=6, "WPGMA"=7, "WPGMC"=8, "UPGMC"=9)


shinyServer(function(input,output){ 

  loadData = function(infile){
    data = read.table(infile, header=F, sep="\t", dec=".", row.names=1)
    colnames(data) <- substr(rownames(data), 1, 25) -> rownames(data)
    postChecking(args, data)
    return (data)
  }
  
  performClassif = function(input){
    data = loadData(input$infile)
    return (getCAH(data, getClassifValue(input$classif_type)))
  }
  
  getClassifValue = function(key)  unlist(classif_methods[key])
  
  output$summary = renderTable({
    data = loadData(input$infile)
    classif = getCAH(data, getClassifValue(input$classif_type))
    printSummary(classif_type, max_cluster, classif, data)
  })
  
  output$fusion_levels = renderPlot({
    data = loadData(input$infile)
    classif = getCAH(data, getClassifValue(input$classif_type))
    plotFusionLevels(classif_type, max_cluster, classif, data)
  })
})
