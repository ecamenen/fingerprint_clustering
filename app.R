rm(list=ls())
#Pseudo-random settings: 
#milisec * PID
set.seed(as.numeric(format(Sys.time(), "%OS2"))*100 * Sys.getpid())

source("fingerprint_clustering.R")
source("server.R")
source("ui")

shinyApp(server = server, ui = ui)