#Remove packages and objects
pkg = names(sessionInfo()$otherPkgs)
pkgs = paste('package:', pkg, sep = "")
lapply(pkgs, detach, character.only = TRUE, unload = TRUE)
rm(list=ls(all=TRUE))

require(tidyverse)
require(plotly)

Sys.setenv("plotly_username"="USERNAME")
Sys.setenv("plotly_api_key"="KEY")

set.seed(123)

#make sure everything groovy
search() #find attached datasets, objects, etc
######

setwd("~/R/working")

list.files()



