rm(list=ls())
#Pseudo-random settings: 
#milisec * PID
set.seed(as.numeric(format(Sys.time(), "%OS2"))*100 * Sys.getpid())

#Loading librairies
librairies = c("shiny","cluster", "optparse", "gclus", "ade4", "scales")
for (l in librairies){
  if (! (l %in% installed.packages()[,"Package"])) install.packages(l, repos = "http://cran.us.r-project.org", quiet = T)
  library(l, character.only = TRUE)
}

source("server.R")
source("ui")

shinyApp(server = server, ui = ui)