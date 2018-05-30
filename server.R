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
  
  #Each shiny server functions run in local environment
  #With assign, variables are forced to be in global env
  setVariables = function(input){
    assign("data", loadData(input$infile), .GlobalEnv)
    assign("classif_type", getClassifValue(input$classif_type), .GlobalEnv)
    assign("classif", getCAH(data, classif_type), .GlobalEnv)
    assign("max_cluster", input$max_clusters, .GlobalEnv)
    assign("nb_clusters", input$nb_clusters, .GlobalEnv)
    assign("advanced", input$advanced, .GlobalEnv)
    assign("verbose", F, .GlobalEnv)
    assign("optimal_nb_clusters", getOptimalNbClus(classif_type, max_cluster, classif, data, nb_clusters), .GlobalEnv)
  }
  
  #t: classif_type
  #m: max_cluster
  #c: classif
  #d: data,
  #n: nb_clusters
  getOptimalNbClus = function(t, m, c, d, n){
    if(n > 1) return (n)
    else return (plotSilhouettePerPart(t, m + 1, c, d))
  }

  # f: filename
  # func: plot function
  savePlot = function(f, func, adv=F){
      print ('ok2')
        f = paste(f, ".pdf", sep="")
        if(isTRUE(adv)) savePdf(f)
        else pdf(f)
        func
        suprLog = dev.off()
  }
  
  # e: event variable
  # f: filename
  # func: plot function
  plot2 = function (e, f, func){
    observeEvent(e, savePlot(f,func))
    func
  }
  
  output$summary = renderTable({
    setVariables(input)
    assign("summary", printSummary(classif_type, max_cluster, classif, data), .GlobalEnv)
    observeEvent(input$summary_save, writeTsv("summary")); summary
  })
  
  output$best_cluster = renderPlot({
    setVariables(input)
    #plotFusionLevels(getClassifValue(input$classif_type), max_cluster, classif, data)
    #TODO: pass an event into a nested function (actually not working)
    # plot2(input$sil_save,
    #       "best_clustering",
    #       plotSilhouettePerPart(classif_type, max_cluster + 1, classif, data))
    func = function() plotSilhouettePerPart(classif_type, max_cluster + 1, classif, data)
    observeEvent(input$best_save, savePlot("best_clustering", func())); func()
  })
  
  output$silhouette = renderPlot({
    
    setVariables(input)
    sil = getSilhouette(classif_type, optimal_nb_clusters, classif, data)
    
    func = function() plotSilhouette(sil)
    observeEvent(input$sil_save, savePlot("silhouette", func())); func()
    
  })
  
  output$pca = renderPlot({
    setVariables(input)
    func = function() plotPca(classif_type, optimal_nb_clusters, classif, data)
    observeEvent(input$pca_save, savePlot("pca", func())); func()
  })
  
  output$heatmap = renderPlot({
    
    setVariables(input)
    par(mar=c(0,0,0,0)) ; plot(0:1,0:1, axes=F, type="n") #delete sil plot
    
    if(classif_type <= 2 | isTRUE(advanced)){
      func = function() heatMap(data, getSilhouette(classif_type, optimal_nb_clusters, classif, data), text=T)
    }else{
      func = function() heatMap(data, c=classif, cl=getClusters(classif_type, optimal_nb_clusters, classif, data), text=T)
    }
    
    observeEvent(input$heatmap_save, savePlot("heatmap", func())); func()
  })
  
  output$cophenetic = renderPlot({
    if(classif_type > 2){
      setVariables(input)
      func = function() plotCohenetic(classif_type, data, classif)
      observeEvent(input$coph_save, savePlot("cohenetic", func())); func()
    }
  })
  
  output$dendrogram = renderPlot({
    if(classif_type > 2){
      setVariables(input)
      func = function() plotDendrogram(classif_type, optimal_nb_clusters, classif, data, advanced)
      observeEvent(input$dendr_save, savePlot("dendrogram", func())); func()
    }
  })
  
  output$ctr_part = renderTable({
    if(isTRUE(input$advanced)){
      setVariables(input)
      assign("summary", printSummary(classif_type, max_cluster, classif, data), .GlobalEnv)
      assign("ctr_part", 
             100 * getPdisPerPartition(classif_type, max_cluster, classif, data),
            .GlobalEnv)
      observeEvent(input$ctr_part_save, writeTsv("ctr_part")); ctr_part
    }
  }, rownames=T, hover=T, striped=T, digits=2, align="c")
  
  output$ctr_clus = renderTable({
    if(isTRUE(input$advanced)){
      setVariables(input)
      assign("ctr_clus", 
             100 * getCtrVar(classif_type, optimal_nb_clusters, classif, data),
            .GlobalEnv)
      observeEvent(input$ctr_clus_save, writeTsv("ctr_clus")); ctr_clus
    }
  }, rownames=T, hover=T, striped=T, digits=2, width="100cm", align="c", size=200)
  
}
