source("fingerprint_clustering.R")
classif_methods <- list("K-menoids" = 1,  "K-means" = 2, "Ward"=3, "Complete links"=4, "Single links"=5, "UPGMA"=6, "WPGMA"=7, "WPGMC"=8, "UPGMC"=9)

loadData = function(f){
  d = read.table(f, header=F, sep="\t", dec=".", row.names=1)
  colnames(d) <- substr(rownames(d), 1, 25) -> rownames(d)
  #postChecking(args, d)
  return (d)
}

server = function(input, output, session){ 
  
  getClassifValue = function(key)  unlist(classif_methods[key])
  
  setVariables = function(input){
    data <<- loadData(input$infile)
    classif_type <<- getClassifValue(input$classif_type)
    classif <<- getCAH(data, classif_type)
    max_cluster <<- input$max_clusters
    nb_clusters <<- input$nb_clusters
    advanced <<- input$advanced
  }
  
  savePlot = function(filename, plot, adv=F){
        filename = paste(filename, ".pdf", sep="")
        if(isTRUE(adv)) savePdf(filename)
        else pdf(filename)
        plot
        suprLog = dev.off()
  }
  
  output$summary = renderTable({
    setVariables(input)
    assign("summary",printSummary(classif_type, max_cluster, classif, data),.GlobalEnv)
    observeEvent(input$summary_save, writeTsv("summary")); summary
  })
  
  output$best_cluster = renderPlot({
    
    setVariables(input)
    best = plotSilhouettePerPart(classif_type, max_cluster + 1, classif, data)
    #plotFusionLevels(getClassifValue(input$classif_type), max_cluster, classif, data)
    
    observeEvent(input$best_save, {
      savePdf("best_clustering.pdf")
      print(classif_type, max_cluster + 1, classif, data)
      plotSilhouettePerPart(classif_type, max_cluster + 1, classif, data)
      suprLog = dev.off()
      #savePlot("best_clustering","best", T)
      }); best
  })
  
  output$silhouette = renderPlot({
    
    setVariables(input)
    
    optimal_nb_clusters = plotSilhouettePerPart(classif_type, max_cluster + 1, classif, data)
    if(nb_clusters > 1) optimal_nb_clusters = nb_clusters
    
    sil = getSilhouette(classif_type, optimal_nb_clusters, classif, data)
    sil_plot = plotSilhouette(sil)
    observeEvent(input$sil_save, savePlot("silhouette",plotSilhouette(sil))); sil_plot
  })
  
  output$pca = renderPlot({
    
    setVariables(input)
    optimal_nb_clusters = plotSilhouettePerPart(classif_type, max_cluster + 1, classif, data)
    if(nb_clusters > 1) optimal_nb_clusters = nb_clusters
    
    plotPca(classif_type, optimal_nb_clusters, classif, data)
  })
  
  output$heatmap = renderPlot({
    
    setVariables(input)
    optimal_nb_clusters <<- plotSilhouettePerPart(classif_type, max_cluster + 1, classif, data)
    if(nb_clusters > 1) optimal_nb_clusters <<- nb_clusters
    #delete sil plot
    par(mar=c(0,0,0,0)) ; plot(0:1,0:1, axes=F, type="n")
    
    if(classif_type <= 2 | isTRUE(advanced)){
      heatMap(data, getSilhouette(classif_type, optimal_nb_clusters, classif, data), text=T)
    }else{
      heatMap(data, c=classif, cl=getClusters(classif_type, optimal_nb_clusters, classif, data), text=T)
    }
  })
  
  output$cophenetic = renderPlot({
    if(classif_type > 2){
      setVariables(input)
      plotCohenetic(classif_type, data, classif)
    }
  })
  
  output$dendrogram = renderPlot({
    if(classif_type > 2){
      setVariables(input)
      optimal_nb_clusters = plotSilhouettePerPart(classif_type, max_cluster + 1, classif, data)
      if(nb_clusters > 1) optimal_nb_clusters = nb_clusters
    
      plotDendrogram(classif_type, optimal_nb_clusters, classif, data, advanced)
    }
  })
  
  output$ctr_part = renderTable({
    if(isTRUE(input$advanced)){
      setVariables(input)
      
      100 * getPdisPerPartition(classif_type, max_cluster, classif, data)
    }
  }, rownames=T, hover=T, striped=T, digits=2, align="c")
  
  output$ctr_clus = renderTable({
    if(isTRUE(input$advanced)){
      setVariables(input)
      optimal_nb_clusters = plotSilhouettePerPart(classif_type, max_cluster + 1, classif, data)
      if(nb_clusters > 1) optimal_nb_clusters = nb_clusters
      
      100 * getCtrVar(classif_type, optimal_nb_clusters, classif, data)
    }
  }, rownames=T, hover=T, striped=T, digits=2, width="100cm", align="c", size=200)
  
}
