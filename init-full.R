#Remove packages and objects
pkg = names(sessionInfo()$otherPkgs)
pkgs = paste('package:', pkg, sep = "")
lapply(pkgs, detach, character.only = TRUE, unload = TRUE)
rm(list=ls(all=TRUE))

require(mosaic)
require(tidyverse)
require(plotly)

Sys.setenv("plotly_username"="USERNAME")
Sys.setenv("plotly_api_key"="KEY")

set.seed(123)
#sanity check to make sure everything groovy
search() #find attached datasets, objects, etc
#searchpaths() #find objects and paths
#sessionInfo() #check everythin else

######

#function for importing entire excel spreadsheet:
#convert xls workbook to list, each xls sheet to list element as tibble
read_excel_allsheets <- function(filename) {
  sheets <- readxl::excel_sheets(filename)
  x <- lapply(sheets, function(X) 
    readxl::read_excel(filename, col_names = TRUE, sheet = X))
  names(x) <- sheets
  x
}

#My boss likes everything in excel spreadsheet,
#so I use a function to aid export to single excel file
#populate list
#each list element will be a new tab in excel spreadsheet
#"name" list is for adding names while writing excel file
excel=list(dataFrame1[,c("col1","col2","col3")])
name=c("dataFrame1")
#append new sheets to list
excel[[length(excel)+1]] <- dataFrame2[,c("col1","col2","col3")]
name=c(name, "dataFrame2")

#alternatively, if you have large list (as seen in loops, etc),
#initialize an empty list to avoid memory fragmentation
excel <- vector("list", 2)
#then populate the list
excel[[1]] <- dataFrame1[,c("col1","col2","col3")]
excel[[2]] <- dataFrame2[,c("col1","col2","col3")]
name=c("dataFrame1", "dataFrame2")

#writing excel file
save2.xlsx <- function (file, x, name)
{
  require(xlsx, quietly = TRUE)
  objects <- x
  names(objects)=name
  nobjects <- length(objects)
  for (i in 1:nobjects) {
    if (i == 1)
      write.xlsx(objects[[i]], file, sheetName = names(objects)[i], row.names=FALSE)
    else write.xlsx(objects[[i]], file, sheetName = names(objects)[i], row.names=FALSE,
                    append = TRUE)
  }
  print(paste("Workbook", file, "has", nobjects, "worksheets."))
}

save2.xlsx("Spreadsheet.xlsx", excel, name)

######
#end seession initialization
######
###########
###########

setwd("~/R/working")

list.files()



