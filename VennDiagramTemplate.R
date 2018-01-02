require(Vennerable)

#Two Venn diagram functions to choose from:

#1) textVenn = venn with labels, counts, etc
#this function is great for exploratory visualization

# ... passes arbitrary number of arguments
# vname = passes name of each variable
textVenn=function(..., vname=deparse(substitute(cat(...)))) {
  #pull all variable names into list to construct file name
  x=list(...)
  valName=sub("cat\\((.*)\\)", "\\1", vname)
  noSep=gsub("(, )", "_", valName)
  fname=paste0(noSep, "_text.pdf")
  
  #make labels for venn construction
  labels=as.vector(unlist(strsplit(valName,", ")),mode="list")
  #make venn
  vennObj=Venn(x, SetNames=c(labels))
  #export venn object to file, pass fname as file name
  pdf(file = fname, width = 6, height = 6)
  plot(vennObj)
  dev.off()
}

#2) blankVenn = same as textVenn(), but with no labels, numbers, etc
#this function is useful for creating images for publications, as text often
#does not scale well when creating manuscript figures.

blankVenn=function(..., vname=deparse(substitute(cat(...)))) {
  #pull all variable names into list to construct file name
  x=list(...)
  valName=sub("cat\\((.*)\\)", "\\1", vname)
  noSep=gsub("(, )", "_", valName)
  fname=paste0(noSep, ".pdf")
  
  #make venn
  vennObj=Venn(x)
  #export venn object to file, pass fname as file name
  pdf(file = fname, width = 6, height = 6)
  plot(vennObj, show = list(SetLabels = FALSE, FaceText="", Faces = TRUE))
  dev.off()
}

#Alternatively, merge textVenn() and blankVenn() functions, so can run both with single call!
allVenns=function(...){
  blankVenn(...)
  textVenn(...)
}



#####
#####
#####
##### Example data

setA = c( "THE1", "DOG", "RAN", "CROSS", "THE2", "ROAD" )
setB = c( "WHY", "DID", "THE1", "CHICKEN", "CROSS", "THE2", "ROAD" )
setC = c( "I", "CHASED", "MY", "DOG", "CROSS", "THE1", "ROAD" )

dframeA = data.frame( setA, key=1:length( setA ) )
dframeB = data.frame( setB, key=1:length( setB ) )
dframeC = data.frame( setC, key=1:length( setC ) )




######
######
######
###### Let Draw some Venn Diagrams!

#Build Venns and export relevant gene lists
setwd( "~/R/working" )

#export simple gene list
CommonToAB = intersect( setA, setB )
write.table(CommonToAB, "CommonToAB.txt", sep = '\t')

#Venn with 2 character vectors
allVenns( setA, setB )

#Venn using columns from 2 dataframes
allVenns( dframeA$setA, dframeB$setB )

#Venn with 3 character vectors
allVenns( setA, setB, setC )

#Venn using columns from 3 dataframes
allVenns( dframeA$setA, dframeB$setB, dframeC$setC )
