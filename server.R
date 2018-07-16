rm(list=ls())
source("fingerprint_clustering.R")
classif_methods <- list("K-menoids" = 1,  "K-means" = 2, "Ward"=3, "Complete links"=4, "Single links"=5, "UPGMA"=6, "WPGMA"=7, "WPGMC"=8, "UPGMC"=9)
library("shinyjs")

tryCatch({
  data <- loadData("matrix.txt")
  }, warning = function(w) {
    data=NULL
    message("Default file \"matrix.txt\" is not in the folder. Please, load another one.")
  }, error = function(e) {
  })

loadData = function(f, s="\t",h=F){
  #!file.exists(
  if(!is.null(f)){
    data = read.table(f,
                 header=h,
                 sep=s,
                 dec=".",
                 row.names=1)
    colnames(data) <- substr(rownames(data), 1, 25) -> rownames(data)
    data = preProcessData(data)
  }
  return (data)
}

server = function(input, output, session){
  
  getClassifValue = function(key)  unlist(classif_methods[key])
  
  #Each shiny server functions run in local environment
  #With assign, variables are forced to be in global env
  setVariables = function(input){
    assign("data", 
           loadData(input$infile$datapath, input$sep, input$header),
           .GlobalEnv)
    assign("CLASSIF_TYPE",
           as.integer(getClassifValue(input$classif_type)),
           .GlobalEnv)
    assign("MAX_CLUSTERS",
           input$max_clusters,
           .GlobalEnv)
    assign("NB_CLUSTERS",
           input$nb_clusters,
           .GlobalEnv)
    assign("ADVANCED",
           input$advanced,
           .GlobalEnv)
    assign("VERBOSE",
           F,
           .GlobalEnv)
    assign("AXIS1",
           input$axis1,
           .GlobalEnv)
    assign("AXIS2",
           input$axis2,
           .GlobalEnv)
    
    #Perform classification
    printProgress(VERBOSE_NIV2, "Distance calculation")
    assign("dis",
           getDistance(data, as.integer(input$dist_type)),
           .GlobalEnv)
    if(CLASSIF_TYPE < 3) printProgress(VERBOSE_NIV2, "Classification")
    assign("classif",
           getClassif(CLASSIF_TYPE, MAX_CLUSTERS, data, dis),
           .GlobalEnv)
    assign("list_clus",
           getClusterPerPart(MAX_CLUSTERS+1, classif),
           .GlobalEnv)
    
    #inertia
    printProgress(VERBOSE_NIV2, "Index calculation")
    assign("between",
           getRelativeBetweenPerPart(MAX_CLUSTERS, data, list_clus),
           .GlobalEnv)
    assign("diff",
           getBetweenDifferences(between),
           .GlobalEnv)
  }
  
  setPrintFuncs = function(){

    ###### plot funcs #####
    assign("plotPCA", 
           function() {
              printProgress(VERBOSE_NIV2, "PCA")
              assign("pca",
                    dudi.pca(data, scannf=F, nf=max(AXIS1,AXIS2)),
                    .GlobalEnv)
                  plotPca(pca, data, clusters, AXIS1, AXIS2)
                  # if(isTRUE(ADVANCED)){
                  #   par(fig=c(0.8,1,0.82,1), new=T)
                  #   plotInertiaPca(pca, d, pca$nf)
                  # }
             },
           .GlobalEnv)
    assign("plotBest", 
           function() plotSilhouettePerPart(mean_silhouette),
           .GlobalEnv)
    assign("plotSil", 
           function() {
             assign("sil_k",
                  sil[[optimal_nb_clusters-1]],
                  .GlobalEnv)
             plotSilhouette(sil_k)
           },
           .GlobalEnv)
    assign("plotCoph", 
           function()  {
             printProgress(VERBOSE_NIV2, "Cophenetic calculation")
             plotCohenetic(dis, classif)
             },
           .GlobalEnv)
    assign("plotDend", 
           function() plotDendrogram(CLASSIF_TYPE, optimal_nb_clusters, classif, data, MAX_CLUSTERS, clusters),
           .GlobalEnv)

    if(CLASSIF_TYPE <= 2 | isTRUE(ADVANCED)){
      assign("plotHeatmap", 
             function() heatMap(data, dis, sil_k),
             .GlobalEnv)
    }else{
      assign("plotHeatmap",
             function() heatMap(data, dis, c=classif, cl=clusters),
             .GlobalEnv)
      
    }
    
    ##### print table func #####

    assign("summary",
           printSummary(between, diff, mean_silhouette, ADVANCED, gap),
           .GlobalEnv)
    assign("ctr_part", 
           100 * getPdisPerPartition(CLASSIF_TYPE, MAX_CLUSTERS, list_clus, data),
           .GlobalEnv)
    assign("ctr_clus", 
           100 * getCtrVar(CLASSIF_TYPE, optimal_nb_clusters, clusters, data),
           .GlobalEnv)
  }
  
  checkMaxCluster = function(){
    assign("sil", 
           getSilhouettePerPart(data, list_clus, dis),
           .GlobalEnv)
    assign("mean_silhouette", 
           getMeanSilhouettePerPart(sil),
           .GlobalEnv)
    
    if(MAX_CLUSTERS > (nrow(data) - 1)){
      message(paste("[WARNING] Max number of clusters must be lower (and not equal) to the number of line of the dataset (", nrow(data), ")", sep=""))
      return (F)
    }else{
      if(NB_CLUSTERS > 0) 
        assign("optimal_nb_clusters",
               NB_CLUSTERS,
               .GlobalEnv)
      else assign("optimal_nb_clusters", 
             which.max(mean_silhouette)+1,
             .GlobalEnv)
      
      #cl_temp, because writeTsv(clusters) recreate a different object named clusters
      assign("clusters", 
             list_clus[[optimal_nb_clusters-1]],
             .GlobalEnv)
      assign("cl_temp", 
             clusters,
             .GlobalEnv)
      #writeClusters(data, clusters, "clusters.tsv", TRUE, v=F )
      assign("gap", 
             NULL,
             .GlobalEnv)
      
      #errors
      if (optimal_nb_clusters==MAX_CLUSTERS) message("\n[WARNING] The optimal number of clusters equals the maximum number of clusters. \nNo cluster structure has been found.")
      if(min(table(cl_temp))==1) message("\n[WARNING] A cluster with an only singleton biased the silhouette score.")
      
      return(T)
    }
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
  
  observeEvent(input$advanced, {
    if(isTRUE(input$advanced))  
      cat(paste("\nAGGLOMERATIVE COEFFICIENT: ", round(getCoefAggl(classif),3), "\n", sep=""))
  })
  
  observe({
    #set an id in tabsetPanel (here "navbar") and for each tabs
    
    #default behaviour
    show(selector = "#navbar li a[data-value=coph]")
    show(selector = "#navbar li a[data-value=dendr]")
    hide(selector = "#navbar li a[data-value=ctr_part]")
    hide(selector = "#navbar li a[data-value=ctr_clus]")
    
    #responsive for a given condition
    toggle(condition = input$advanced, selector = "#navbar li a[data-value=ctr_part]")
    toggle(condition = input$advanced, selector = "#navbar li a[data-value=ctr_clus]")
    toggle(condition = (as.integer(getClassifValue(input$CLASSIF_TYPE) > 2)), selector = "#navbar li a[data-value=coph]")
    toggle(condition = (as.integer(getClassifValue(input$CLASSIF_TYPE) > 2)), selector = "#navbar li a[data-value=dendr]")
  })
  
  observeEvent(input$save_all, {
    if(is.data.frame(data)){
      setVariables(input)
      setPrintFuncs()
      writeTsv("summary", "summary.tsv", v=F)
      
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
  
  output$summary = renderTable({
    tryCatch({
      setVariables(input)
      if(checkMaxCluster()){
        setPrintFuncs()
        observeEvent(input$summary_save, writeTsv("summary", "summary.tsv", v=F)); summary
      }
    }, error = function(e) {
    })
  })
  
  output$best_cluster = renderPlot({

    tryCatch({
      setVariables(input)

      if(checkMaxCluster()){
        setPrintFuncs()
        #plotFusionLevels(getClassifValue(input$CLASSIF_TYPE), MAX_CLUSTERS, classif, data)
        #TODO: pass an event into a nested function (actually not working)
        # plot2(input$sil_save,
        #       "best_clustering",
        #       plotSilhouettePerPart(CLASSIF_TYPE, MAX_CLUSTERS + 1, classif, data))
        observeEvent(input$best_save, savePlot("best_clustering", plotBest()))
        plotBest()

      }
    }, error = function(e) {
    })
  })
  
  output$silhouette = renderPlot({
    tryCatch({
      setVariables(input)
      if(checkMaxCluster()){
        setPrintFuncs()
        observeEvent(input$sil_save, savePlot("silhouette", plotSil()))
        plotSil()
      }
    }, error = function(e) {
    })
  })
  
  output$pca = renderPlot({
    tryCatch({
      setVariables(input)
      if(checkMaxCluster()){
        setPrintFuncs()
        observeEvent(input$pca_save, savePlot("pca", plotPCA()))
        plotPCA()
      }
    }, error = function(e) {
    })
  })
  
  output$heatmap = renderPlot({
    tryCatch({
      setVariables(input)
      if(checkMaxCluster()){
        par(mar=c(0,0,0,0)) ; plot(0:1,0:1, axes=F, type="n") #delete sil plot
        setPrintFuncs()
        observeEvent(input$heatmap_save, savePlot("heatmap", plotHeatmap()))
        plotHeatmap()
      }
    }, error = function(e) {
    })
  })
  
  output$cophenetic = renderPlot({
    tryCatch({
      if( CLASSIF_TYPE > 2){
        setVariables(input)
        if(checkMaxCluster()){
          setPrintFuncs()
          observeEvent(input$coph_save, savePlot("cohenetic", plotCoph()))
          plotCoph()
        }
      }
    }, error = function(e) {
    })
  })
  
  output$dendrogram = renderPlot({
    tryCatch({
      if( CLASSIF_TYPE > 2){
        setVariables(input)
        if(checkMaxCluster()){
          setPrintFuncs()
          observeEvent(input$dendr_save, savePlot("dendrogram", plotDend()))
          plotDend()
        }
      }
    }, error = function(e) {
    })
  })
  
  output$ctr_part = renderTable({
    tryCatch({
      if(isTRUE(input$advanced)){
        setVariables(input)
        if(checkMaxCluster()){
          setPrintFuncs()
          observeEvent(input$ctr_part_save, writeTsv("ctr_part"))
          ctr_part
        }
      }
    }, error = function(e) {
    })
  }, rownames=T, hover=T, striped=T, digits=2, align="c")
  
  output$ctr_clus = renderTable({
    tryCatch({
      if(isTRUE(input$advanced)){
        setVariables(input)
        if(checkMaxCluster()){
          setPrintFuncs()
          observeEvent(input$ctr_clus_save, writeTsv("ctr_clus"))
          ctr_clus
        }
      }
    }, error = function(e) {
    })
  }, rownames=T, hover=T, striped=T, digits=2, width="100cm", align="c", size=200)
  
}
