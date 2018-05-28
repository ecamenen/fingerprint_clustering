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
  
  setVariables = function(input){
    max_cluster = input$max_clusters
    nb_clusters = input$nb_clusters
  }
  
  output$summary = renderTable({
    data = loadData(input$infile)
    classif = getCAH(data, getClassifValue(input$classif_type))
    printSummary(getClassifValue(input$classif_type), max_cluster, classif, data)
  })
  
  output$fusion_levels = renderPlot({
    
    data = loadData(input$infile)
    
    #setVariables(input)
    classif_type = getClassifValue(input$classif_type)
    max_cluster = input$max_clusters
    nb_clusters = input$nb_clusters
    classif = getCAH(data, classif_type)

    plotSilhouettePerPart(classif_type, max_cluster + 1, classif, data)
    #plotFusionLevels(getClassifValue(input$classif_type), max_cluster, classif, data)
  })
  
  output$silhouette = renderPlot({
    
    data = loadData(input$infile)
    
    classif_type = getClassifValue(input$classif_type)
    max_cluster = input$max_clusters
    nb_clusters = input$nb_clusters
    classif = getCAH(data, classif_type)
    
    optimal_nb_clusters = plotSilhouettePerPart(classif_type, max_cluster + 1, classif, data)
    if(nb_clusters > 1) optimal_nb_clusters = nb_clusters

    sil = getSilhouette(classif_type, optimal_nb_clusters, classif, data)
    plotSilhouette(sil)
  })
  
  output$pca = renderPlot({
    
    data = loadData(input$infile)
    
    classif_type = getClassifValue(input$classif_type)
    max_cluster = input$max_clusters
    nb_clusters = input$nb_clusters
    classif = getCAH(data, classif_type)
    optimal_nb_clusters = plotSilhouettePerPart(classif_type, max_cluster + 1, classif, data)
    if(nb_clusters > 1) optimal_nb_clusters = nb_clusters
    
    plotPca(classif_type, optimal_nb_clusters, classif, data)
  })
  
  output$heatmap = renderPlot({
    
    data = loadData(input$infile)
    
    classif_type = getClassifValue(input$classif_type)
    max_cluster = input$max_clusters
    nb_clusters = input$nb_clusters
    classif = getCAH(data, classif_type)
    optimal_nb_clusters = plotSilhouettePerPart(classif_type, max_cluster + 1, classif, data)
    if(nb_clusters > 1) optimal_nb_clusters = nb_clusters
    #delete sil plot
    par(mar=c(0,0,0,0)) ; plot(0:1,0:1, axes=F, type="n")
    
    if(classif_type <= 2 | isTRUE(advanced)){
      sil = getSilhouette(classif_type, optimal_nb_clusters, classif, data)
      heatMap(data, sil, text=T)
    }else{
      clusters = getClusters(classif_type, optimal_nb_clusters, classif, data)
      heatMap(data, c=classif, cl=clusters, text=T)
    }
  })
  
  output$cophenetic = renderPlot({

    data = loadData(input$infile)
    
    classif_type = getClassifValue(input$classif_type)
    max_cluster = input$max_clusters
    nb_clusters = input$nb_clusters
    classif = getCAH(data, classif_type)
    
    if(classif_type > 2){
      plotCohenetic(classif_type, data, classif)
    }
  })
  
  output$dendrogram = renderPlot({
    
    data = loadData(input$infile)
    
    classif_type = getClassifValue(input$classif_type)
    max_cluster = input$max_clusters
    nb_clusters = input$nb_clusters
    classif = getCAH(data, classif_type)
    optimal_nb_clusters = plotSilhouettePerPart(classif_type, max_cluster + 1, classif, data)
    if(nb_clusters > 1) optimal_nb_clusters = nb_clusters
    
    if(classif_type > 2){
      plotDendrogram(classif_type, optimal_nb_clusters, classif, data, advanced)
    }
  })
  
})
