library(shiny)
source("fingerprint_clustering.R")

shinyServer(function(input,output){ #les elements de la partie server vont etre
  #definies ici
  
  
  #
    
  output$silhouette = renderTable({
    #fonction reactive (render) shiny pour permettre lien dynamique et afficher le texte
    data = read.table(input$infile, header=F, sep="\t", dec=".", row.names=1)
    colnames(data) <- substr(rownames(data), 1, 25) -> rownames(data)
    postChecking(args, data)
    
    #Perform classification
    if(classif_type == 0) classif_type = selectBestCAH(data, F)
    classif = getCAH(data, classif_type)

    if(classif_type>2) plotCohenetic(classif_type, data, classif)
    plotFusionLevels(classif_type, max_cluster, classif, data)
    
    #Silhouette analysis
    optimal_nb_clusters = plotSilhouettePerPart(classif_type, max_cluster + 1, classif, data)
    if(!is.null(nb_clusters)) optimal_nb_clusters = nb_clusters
    sil = getSilhouette(classif_type, optimal_nb_clusters, classif, data)
    plotSilhouette(sil)
    summary = printSummary(classif_type, max_cluster, classif, data)
    
    #name2=getClassifType(classif_type)
    #paste0(hist(rnorm(100), main=name2))
    #paste0(plotCohenetic(classif_type, data, classif))
    summary
    })
})
