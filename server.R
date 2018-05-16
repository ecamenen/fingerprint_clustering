library(shiny)
source("fingerprint_clustering.R")

shinyServer(function(input,output){ 
  
  loadData = function(infile){
    data = read.table(infile, header=F, sep="\t", dec=".", row.names=1)
    colnames(data) <- substr(rownames(data), 1, 25) -> rownames(data)
    postChecking(args, data)
    return (data)
  }
  
  performClassif = function(data){
    #Perform classification
    if(classif_type == 0) classif_type = selectBestCAH(data, F)
    return(getCAH(data, classif_type))
  }
  
  output$summary = renderTable({
    data = loadData(input$infile)
    classif = performClassif(data)
    printSummary(classif_type, max_cluster, classif, data)
  })
  
  output$fusion_levels = renderPlot({
    data = loadData(input$infile)
    classif = performClassif(data)
    plotFusionLevels(classif_type, max_cluster, classif, data)
  })
})
