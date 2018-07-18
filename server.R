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
    colnames(data) <- substr(rownames(data), 1, MAX_CHAR_LEN) -> rownames(data)
    data = preProcessData(data)
  }
  return (data)
}

server = function(input, output, session){

  getClassifValue = function(key)  unlist(classif_methods[key])
  
  ###################################
  #          SETTINGS
  ###################################
  
  #Each shiny server functions run in local environment
  #With assign, variables are forced to be in global env
  setVariables = function(input){

    assign("MAX_CLUSTERS",
           input$max_clusters,
           .GlobalEnv)
    assign("CLASSIF_TYPE",
           as.integer(getClassifValue(input$classif_type)),
           .GlobalEnv)
    
    assign("NB_CLUSTERS",
           input$nb_clusters,
           .GlobalEnv)
    assign("ADVANCED",
           input$advanced,
           .GlobalEnv)

    
    #Perform classification
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
    
    printProgress(VERBOSE_NIV2, "PCA") 
    assign("pca", 
           dudi.pca(data, scannf=F, nf=4), 
           .GlobalEnv)
  }
  
  setPrintFuncs = function(){
    
    assign("sil_k",
           sil[[optimal_nb_clusters-1]],
           .GlobalEnv)

    ###### plot funcs #####
    assign("plotPCA", 
           function() plotPca(pca, data, clusters, AXIS1, AXIS2),
           .GlobalEnv)
    assign("plotBest", 
           function() plotSilhouettePerPart(mean_silhouette),
           .GlobalEnv)
    assign("plotSil", 
           function() plotSilhouette(sil_k),
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
    
    ##### advanced #####
    if(isTRUE(ADVANCED)){
      
      if (nrow(data) < (NB_ROW_MAX/2)){
        printProgress(VERBOSE_NIV2, "Gap statistics calculation")
        assign("gap",
               getGapPerPart(MAX_CLUSTERS, data, classif, NB_BOOTSTRAP),
               .GlobalEnv)
      }
      
      assign("plotFus",
             function() plotFusionLevels(MAX_CLUSTERS, classif),
             .GlobalEnv)
      assign("plotCoph", 
             function()  {
               printProgress(VERBOSE_NIV2, "Cophenetic calculation")
               plotCohenetic(dis, classif)
             },
             .GlobalEnv)
      assign("plotGap",
             function() {
               if (nrow(data) < (NB_ROW_MAX/2)){
                 plotGapPerPart(gap, MAX_CLUSTERS, v=F)
                 #plotGapPerPart2(gap, MAX_CLUSTERS)
               }
              },
             .GlobalEnv)
      assign("plotElb",
             function() plotElbow(between),
             .GlobalEnv)
      assign("within_k",
             getRelativeWithinPerCluster(list_clus, data),
             .GlobalEnv)
    }
    
    ##### print table func #####

    assign("summary", 
           printSummary(between, diff, mean_silhouette, ADVANCED, gap),
           .GlobalEnv)
    # assign("ctr_part", 
    #        100 * getPdisPerPartition(CLASSIF_TYPE, MAX_CLUSTERS, list_clus, data),
    #        .GlobalEnv)
    # assign("ctr_clus", 
    #        100 * getCtrVar(CLASSIF_TYPE, optimal_nb_clusters, clusters, data),
    #        .GlobalEnv)
    
    writeClusters("clusters.tsv", v=F)
  }
  
  # post-process for data
  # check that the maximum_number of clusters fixed is not greater than the number of row of the datafile
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
      assign("gap", 
             NULL,
             .GlobalEnv)
      
      #errors
      if (optimal_nb_clusters==MAX_CLUSTERS) message("\n[WARNING] The optimal number of clusters equals the maximum number of clusters. \nNo cluster structure has been found.")
      if(min(table(clusters))==1) message("\n[WARNING] A cluster with an only singleton biased the silhouette score.")
      
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

  setDataNDistance = reactive({
    assign("data",
           loadData(input$infile$datapath, input$sep, input$header),
           .GlobalEnv)
    printProgress(VERBOSE_NIV2, "Distance calculation")
    assign("dis",
           getDistance(data, as.integer(input$dist_type)),
           .GlobalEnv)
    setClassifPar()
  })
  
  setClassifPar = reactive({
    assign("MAX_CLUSTERS",
           input$max_clusters,
           .GlobalEnv)
    assign("CLASSIF_TYPE",
           as.integer(getClassifValue(input$classif_type)),
           .GlobalEnv)
    print('ok')
  })
  
  ###################################
  #          EVENTS
  ###################################
  
  observeEvent(input$infile, {
    assign("HEAD", 
           FALSE,
           .GlobalEnv)
    setDataNDistance()
  })
  
  # observeEvent(c(input$max_clusters, input$classif_type), {
  #     print('ok')
  #     setClassifPar()
  # })
  
  observeEvent(input$header, {
    assign("HEAD",
           input$header,
           .GlobalEnv)
    if(!is.null(input$infile)){
      setDataNDistance()
    }
  })

  #events for advanced mode
  observeEvent(input$advanced, {
    if(!is.null(input$infile)){
      #cat(paste("\nAGGLOMERATIVE COEFFICIENT: ", round(getCoefAggl(classif),3), "\n", sep=""))
    }
  })

  # hide either input options or tabs
  observe({
    
    #set an id in tabsetPanel (here "navbar") and for each tabs
      
      #default behaviour
      show(selector = "#navbar li a[data-value=dendr]")
      hide(selector = "#navbar li a[data-value=elbow]")
      hide(selector = "#navbar li a[data-value=gap]")
      hide(selector = "#navbar li a[data-value=coph]")
      hide(selector = "#navbar li a[data-value=fusion]")
      hide(selector = "#navbar li a[data-value=within]")
      toggle(condition = input$advanced, id="nb_clusters")
      toggle(condition = input$advanced, id="axis1")
      toggle(condition = input$advanced, id="axis2")
      
      if(!is.null(input$infile)){ #catch condition when no data are loaded and adv is selected
      
        #responsive for a given condition
        toggle(condition = input$advanced, selector = "#navbar li a[data-value=elbow]")
        toggle(condition = input$advanced, selector = "#navbar li a[data-value=gap]")
        toggle(condition = input$advanced, selector = "#navbar li a[data-value=within]")
        toggle(condition = ( input$advanced & as.integer(getClassifValue(input$classif_type) > 2) ), selector = "#navbar li a[data-value=fusion]")
        toggle(condition = ( input$advanced & as.integer(getClassifValue(input$classif_type) > 2) ), selector = "#navbar li a[data-value=coph]")
        toggle(condition = (as.integer(getClassifValue(input$classif_type) > 2)), selector = "#navbar li a[data-value=dendr]")
        
      }
  })
  
  observeEvent(input$refresh, {
      js$refresh()
  })
  
  #function save_all
  observeEvent(input$save_all, {
    if(is.data.frame(data)){
      setVariables(input)
      setPrintFuncs()
      writeTsv("summary", "summary.tsv", v=F)
      
      savePlot("best_clustering", plotBest())
      savePlot("silhouette", plotSil())
      savePlot("pca", plotPCA())
      savePlot("heatmap", plotHeatmap())
      
      if(as.integer(getClassifValue(input$classif_type) > 2)){
        savePlot("dendrogram", plotDend())
      }
      
      if(isTRUE(input$advanced)){
        if(as.integer(getClassifValue(input$classif_type) > 2)){
          savePlot("cohenetic", plotCoph())
          savePlot("fusion", plotFus())
        }
        savePlot("elbow", plotElb())
        #if (nrow(data) < (NB_ROW_MAX/2)){
          savePlot("gap", plotGap())
        #}
        # writeTsv("ctr_clus", "ctr_clus.tsv", v=F)
        # writeTsv("ctr_part", "ctr_part.tsv", v=F)
        writeTsv("within_k", "within_k.tsv", v=F)
      }
    }
  })
  
  ###################################
  #          RENDER
  ###################################
  
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
        assign("AXIS1",
               input$axis1,
               .GlobalEnv)
        assign("AXIS2",
               input$axis2,
               .GlobalEnv)
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
      if(isTRUE(input$advanced) & CLASSIF_TYPE > 2){
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
  
  output$fusion = renderPlot({
    tryCatch({
      if(isTRUE(input$advanced) & CLASSIF_TYPE > 2){
        setVariables(input)
        if(checkMaxCluster()){
          setPrintFuncs()
          observeEvent(input$fusion_save, savePlot("fusion", plotFus()))
          plotFus()
        }
      }
    }, error = function(e) {
    })
  })
  
  output$gap = renderPlot({
    tryCatch({
      if(isTRUE(input$advanced)){
        setVariables(input)
        if(checkMaxCluster()){
          setPrintFuncs()
          observeEvent(input$gap_save, savePlot("gap", plotGap()))
          plotGap()
        }
      }
    }, error = function(e) {
    })
  })
  
  output$elbow = renderPlot({
    tryCatch({
      if(isTRUE(input$advanced)){
        setVariables(input)
        if(checkMaxCluster()){
          setPrintFuncs()
          observeEvent(input$elbow_save, savePlot("elbow", plotElb()))
          plotElb()
          #plotFus()
          #print("Gap statistics calculation in progress...")
          #plotGap()
        }
      }
    }, error = function(e) {
    })
  })
  
  output$within = renderTable({
    tryCatch({
      if(isTRUE(input$advanced)){
        setVariables(input)
        if(checkMaxCluster()){
          setPrintFuncs()
          observeEvent(input$within_save, writeTsv("within_k", "within_k.tsv", v=F) )
          within_k
        }
      }
    }, error = function(e) {
    })
  }, rownames=T, hover=T, striped=T, digits=2, width="100cm", align="c", na="", size=200)
  
}
