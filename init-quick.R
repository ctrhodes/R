#Remove packages and objects
pkg = names(sessionInfo()$otherPkgs)
pkgs = paste('package:', pkg, sep = "")
lapply(pkgs, detach, character.only = TRUE, unload = TRUE)
rm(list=ls(all=TRUE))

require(tidyverse)
require(plotly)

Sys.setenv("plotly_username"="crplotly")
Sys.setenv("plotly_api_key"="c6t0tpuq7j")

set.seed(123)

#sanity check to make sure everything groovy
search() #find attached datasets, objects, etc
######

setwd("~/R/working")

list.files()



