closeAllConnections()
rm(list=ls(all=TRUE))
source("fingerprint_clustering.R")

classif_methods <- list("K-menoids" = 1,  "K-means" = 2, "Ward"=3, "Complete links"=4, "Single links"=5, "UPGMA"=6, "WPGMA"=7, "WPGMC"=8, "UPGMC"=9)
d=NULL

loadData = function(f, s="\t"){
  
  if(!is.null(f)){
    d = read.table(f,
       header=F,
       sep=s,
       dec=".",
       row.names=1)
    colnames(d) <- substr(rownames(d), 1, 100) -> rownames(d)
  }
  return (d)
}

assign("fingerprint",
       F,
       .GlobalEnv)

assign("fingerprintFile",
       "nofile",
       .GlobalEnv)

assign("algoShortestPath",
       "ValidShortest",
       .GlobalEnv)

server = function(input, output, session){
  
  donnees <- reactiveValues()
  donnees$data <- NULL
  
  hide("loading-content")
  hide("tabsetPanel")
  getClassifValue = function(key)  unlist(classif_methods[key])
  
  #Each shiny server functions run in local environment
  #With assign, variables are forced to be in global env
  setVariables = function(input){
    
    assign("fingerprint",
           input$fingerprint,
           .GlobalEnv)
    
    assign("fingerprintFile",
           input$fingerprintFile,
           .GlobalEnv)
    
    assign("fingerprintFileName",
           input$infile$name,
           .GlobalEnv)
    
    assign("algoShortestPath",
           "input$algoShortestPath",
           .GlobalEnv)
    
    assign("classif_type",
           getClassifValue(input$classif_type),
           .GlobalEnv)
    assign("max_cluster",
           input$max_clusters,
           .GlobalEnv)
    
    assign("nb_clusters",
           input$nb_clusters,
           .GlobalEnv)
    assign("advanced",
           input$advanced,
           .GlobalEnv)
    
    assign("verbose",
           F,
           .GlobalEnv)
    
    if(!is.null(donnees$data)){
      assign("classif",
             getCAH(donnees$data, classif_type),
             .GlobalEnv)
      
      assign("summary",
             printSummary(classif_type, input$max_clusters, classif, donnees$data),
             .GlobalEnv)
    }
  }

  setPrintFuncs = function(){
    message("setPrintFuncs")
    ###### plot funcs #####
    assign("plotPCA", 
           function() plotPca(classif_type, optimal_nb_clusters, classif, donnees$data),
           .GlobalEnv)
    assign("plotBest", 
           function() plotSilhouettePerPart(classif_type, max_cluster + 1, classif, donnees$data),
           .GlobalEnv)
    assign("plotSil", 
           function() plotSilhouette(getSilhouette(classif_type, optimal_nb_clusters, classif, donnees$data)),
           .GlobalEnv)
    assign("plotCoph", 
           function() plotCohenetic(classif_type, donnees$data, classif),
           .GlobalEnv)
    assign("plotDend", 
           function() plotDendrogram(classif_type, optimal_nb_clusters, classif, donnees$data, advanced),
           .GlobalEnv)
    
    if(classif_type <= 2 | isTRUE(advanced)){
      assign("plotHeatmap", 
             function() heatMap(donnees$data,
                                getSilhouette(classif_type, optimal_nb_clusters, classif, donnees$data),
                                text=T),
             .GlobalEnv)
    }else{
      assign("plotHeatmap",
             function() heatMap(donnees$data,
                                c=classif,
                                cl=getClusters(classif_type, optimal_nb_clusters, classif, donnees$data),
                                text=T),
             .GlobalEnv)
    }
    
    ##### print table func #####
    
    assign("summary",
           printSummary(classif_type, max_cluster, classif, donnees$data),
           .GlobalEnv)
    assign("ctr_part", 
           100 * getPdisPerPartition(classif_type, max_cluster, classif, donnees$data),
           .GlobalEnv)
    assign("ctr_clus", 
           100 * getCtrVar(classif_type, optimal_nb_clusters, classif, donnees$data),
           .GlobalEnv)
    
    clusters = getClusters(classif_type, optimal_nb_clusters, classif, donnees$data)
    writeClusters(clusters, T)
  }
  
  checkMaxCluster = function(){
    
    if(max_cluster > (nrow(donnees$data) - 1)){
      message(paste("[WARNING] Max number of clusters must be lower (and not equal) to the number of line of the dataset (", nrow(donnees$data), ")", sep=""))
      return (F)
    }else{
      assign("optimal_nb_clusters", 
             getOptimalNbClus(classif_type, max_cluster, classif, donnees$data, nb_clusters),
             .GlobalEnv)
      return(T)
    }
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
  
  #BUGED
  # e: event variable
  # f: filename
  # func: plot function
  plot2 = function (e, f, func){
    #BUG: cannot get e
    observeEvent(e, savePlot(f,func))
    func
  }
  
  computeDistance = function (){
    hide("no-content")
    show("loading-content") # make the loading pane appear
    
    f=gsub(".txt", "", fingerprintFileName)
    message(f)
    fn<-paste('matrix_', f,'_',algoShortestPath,'.txt', sep = "")
    fileName <- fn
    message(fileName)
    message(fingerprintFile)
    
    if(!file.exists(fileName))
    { 
      status <- exec_wait(paste('java -jar FingerprintSubnetwork-1.1.jar -network "recon2.v03_ext_noCompartment_noTransport.xml" -fingerprint ', fingerprintFile, ' -atommapping "recon2.v03_ext_noCompartment_noTransport_C-AAM-weights.tab" -matrixresult "', fileName, '" -reactionresult "reactionsPath.txt" -metabInfo "metabInfo.tsv" -algo ', algoShortestPath, sep = ""))
      
      while (!file.exists(fileName) & !file.access("fileName", mode = 2)) {
        Sys.sleep(2)
      }
    }  
    hide("loading-content")
    show("tabsetPanel")
    
    donnees$data <- loadData(fileName, "\t")
  }
  
  observeEvent(input$save_all, {
    if(!is.null(donnees$data)){
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
    }
  })
  
  observeEvent(input$conputeShortestPath, {
    if(!is.null(input$algoShortestPath) & !is.null(fingerprintFile)){
      computeDistance()
      Sys.sleep(2)
    }
  })
  
  observeEvent(input$infile, {
    if(!identical(fingerprintFile,input$infile$datapath) & !is.null(input$infile$datapath)){
      assign("fingerprint",
             T,
             .GlobalEnv)
      
      assign("fingerprintFile",
             input$infile$datapath,
             .GlobalEnv)
      
      assign("fingerprintFileName",
             input$infile$name,
             .GlobalEnv)
    }
  })
  
  observeEvent(input$algoShortestPath, {
    if(!identical(algoShortestPath,input$algoShortestPath) & !is.null(input$algoShortestPath) & !is.null(input$infile$datapath)){
      assign("algoShortestPath",
             input$algoShortestPath,
             .GlobalEnv)
    }
  })
  
  output$summary = renderTable({
    if(!is.null(donnees$data)){
      message('summ')
      setVariables(input)
      if(checkMaxCluster()){
        setPrintFuncs()
        observeEvent(input$summary_save, writeTsv("summary")); summary
      }
    }
  })
  
  output$best_cluster = renderPlot({
    if(!is.null(donnees$data)){
      setVariables(input)
      if(checkMaxCluster()){
        setPrintFuncs()
        observeEvent(input$best_save, savePlot("best_clustering", plotBest()))
        plotBest()
      }
    }
  })
  
  output$silhouette = renderPlot({
    if(!is.null(donnees$data)){
      setVariables(input)
      if(checkMaxCluster()){
        setPrintFuncs()
        observeEvent(input$sil_save, savePlot("silhouette", plotSil()))
        plotSil()
      }
    }
  })
  
  output$pca = renderPlot({
    if(!is.null(donnees$data)){
      setVariables(input)
      if(checkMaxCluster()){
        setPrintFuncs()
        observeEvent(input$pca_save, savePlot("pca", plotPCA()))
        plotPCA()
      }
    }
  })
  
  output$heatmap = renderPlot({
    if(!is.null(donnees$data)){
      setVariables(input)
      if(checkMaxCluster()){
        par(mar=c(0,0,0,0)) ; plot(0:1,0:1, axes=F, type="n") #delete sil plot
        setPrintFuncs()
        observeEvent(input$heatmap_save, savePlot("heatmap", plotHeatmap()))
        plotHeatmap()
      }
    }
  })
  
  
  output$cophenetic = renderPlot({
    if( !is.null(donnees$data) & classif_type > 2){
      setVariables(input)
      if(checkMaxCluster()){
        setPrintFuncs()
        observeEvent(input$coph_save, savePlot("cohenetic", plotCoph()))
        plotCoph()
      }
    }
  })
  
  output$dendrogram = renderPlot({
    if(!is.null(donnees$data) & classif_type > 2){
      setVariables(input)
      if(checkMaxCluster()){
        setPrintFuncs()
        observeEvent(input$dendr_save, savePlot("dendrogram", plotDend()))
        plotDend()
      }
    }
  })
  
  output$ctr_part = renderTable({
    if(!is.null(donnees$data) & isTRUE(input$advanced)){
      setVariables(input)
      if(checkMaxCluster()){
        setPrintFuncs()
        observeEvent(input$ctr_part_save, writeTsv("ctr_part"))
        ctr_part
      }
    }
  }, rownames=T, hover=T, striped=T, digits=2, align="c")
  
  output$ctr_clus = renderTable({
    if(!is.null(donnees$data) & isTRUE(input$advanced)){
      setVariables(input)
      if(checkMaxCluster()){
        setPrintFuncs()
        observeEvent(input$ctr_clus_save, writeTsv("ctr_clus"))
        ctr_clus
      }
    }
  }, rownames=T, hover=T, striped=T, digits=2, width="100cm", align="c", size=200)
}