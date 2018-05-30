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
    assign("optimal_nb_clusters", 
           getOptimalNbClus(classif_type, max_cluster, classif, data, nb_clusters),
           .GlobalEnv)
  }
  
  setPrintFuncs = function(){

    ###### plot funcs #####
    assign("plotPCA", 
           function() plotPca(classif_type, optimal_nb_clusters, classif, data),
           .GlobalEnv)
    assign("plotBest", 
           function() plotSilhouettePerPart(classif_type, max_cluster + 1, classif, data),
           .GlobalEnv)
    assign("plotSil", 
           function() plotSilhouette(getSilhouette(classif_type, optimal_nb_clusters, classif, data)),
           .GlobalEnv)
    assign("plotCoph", 
           function() plotCohenetic(classif_type, data, classif),
           .GlobalEnv)
    assign("plotDend", 
           function() plotDendrogram(classif_type, optimal_nb_clusters, classif, data, advanced),
           .GlobalEnv)

    if(classif_type <= 2 | isTRUE(advanced)){
      assign("plotHeatmap", 
             function() heatMap(data, getSilhouette(classif_type, optimal_nb_clusters, classif, data), text=T),
             .GlobalEnv)
    }else{
      assign("plotHeatmap",
             function() heatMap(data, c=classif, cl=getClusters(classif_type, optimal_nb_clusters, classif, data), text=T),
             .GlobalEnv)
    }
    
    ##### print table func #####

    assign("summary",
           printSummary(classif_type, max_cluster, classif, data),
           .GlobalEnv)
    assign("ctr_part", 
           100 * getPdisPerPartition(classif_type, max_cluster, classif, data),
           .GlobalEnv)
    assign("ctr_clus", 
           100 * getCtrVar(classif_type, optimal_nb_clusters, classif, data),
           .GlobalEnv)
    options(warn = 0)
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
    #BUG: cannot get e
    observeEvent(e, savePlot(f,func))
    func
  }
  
  observeEvent(input$save_all, {
    setVariables(input)
    setPrintFuncs()
    writeTsv("summary")
    
    savePlot("best_clustering", plotBest())
    savePlot("silhouette", plotSil())
    savePlot("pca", plotPCA())
    savePlot("heatmap", plotHeatmap())
    savePlot("cohenetic", plotCoph())
    savePlot("dendrogram", plotDend())
    
    if(isTRUE(input$advanced)){
      writeTsv("ctr_clus")
      writeTsv("ctr_part")
    }
  })
  
  output$summary = renderTable({
    setVariables(input)
    setPrintFuncs()
    observeEvent(input$summary_save, writeTsv("summary")); summary
  })
  
  output$best_cluster = renderPlot({
    setVariables(input)
    setPrintFuncs()
    #plotFusionLevels(getClassifValue(input$classif_type), max_cluster, classif, data)
    #TODO: pass an event into a nested function (actually not working)
    # plot2(input$sil_save,
    #       "best_clustering",
    #       plotSilhouettePerPart(classif_type, max_cluster + 1, classif, data))
    observeEvent(input$best_save, savePlot("best_clustering", plotBest())); plotBest()
  })
  
  output$silhouette = renderPlot({
    setVariables(input)
    setPrintFuncs()
    observeEvent(input$sil_save, savePlot("silhouette", plotSil())); plotSil()
    
  })
  
  output$pca = renderPlot({
    setVariables(input)
    setPrintFuncs()
    observeEvent(input$pca_save, savePlot("pca", plotPCA())); plotPCA()
  })
  
  output$heatmap = renderPlot({
    setVariables(input)
    par(mar=c(0,0,0,0)) ; plot(0:1,0:1, axes=F, type="n") #delete sil plot
    setPrintFuncs()
    observeEvent(input$heatmap_save, savePlot("heatmap", plotHeatmap())); plotHeatmap()
    
  })
  
  output$cophenetic = renderPlot({
    if(classif_type > 2){
      setVariables(input)
      setPrintFuncs()
      observeEvent(input$coph_save, savePlot("cohenetic", plotCoph())); plotCoph()
    }
  })
  
  output$dendrogram = renderPlot({
    if(classif_type > 2){
      setVariables(input)
      setPrintFuncs()
      observeEvent(input$dendr_save, savePlot("dendrogram", plotDend())); plotDend()
    }
  })
  
  output$ctr_part = renderTable({
    if(isTRUE(input$advanced)){
      setVariables(input)
      setPrintFuncs()
      observeEvent(input$ctr_part_save, writeTsv("ctr_part")); ctr_part
    }
  }, rownames=T, hover=T, striped=T, digits=2, align="c")
  
  output$ctr_clus = renderTable({
    if(isTRUE(input$advanced)){
      setVariables(input)
      observeEvent(input$ctr_clus_save, writeTsv("ctr_clus")); ctr_clus
    }
  }, rownames=T, hover=T, striped=T, digits=2, width="100cm", align="c", size=200)
  
}
