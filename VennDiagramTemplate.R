#clear memory, with exception of attached objects (i.e. attach())
rm(list=ls(all=TRUE))
require(Vennerable)

########
########
#start venn function
#textVenn = venn with labels, counts, etc
# ... passes arbitrary number of arguments
# vname = passes name of each variable/set
textVenn=function(..., vname=deparse(substitute(cat(...)))) {
  #pull all variable names into list to construct file name
  x=list(...)
  valName=sub("cat\\((.*)\\)", "\\1", vname)
  noSep=gsub("(, )", "_", valName)
  fname=paste0(noSep, "_names.pdf")
  
  #make labels for venn construction
  labels=as.vector(unlist(strsplit(valName,", ")),mode="list")
  #make venn
  vennObj=Venn(x, SetNames=c(labels))
  #export venn object to file, pass fname as file name
  pdf(file = fname, width = 6, height = 6)
  plot(vennObj)
  dev.off()
}

#blankVenn = same as textVenn, but with no labels, numbers, etc
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

#merge functions, so can run both with single call
allVenns=function(...){
  blankVenn(...)
  textVenn(...)
}

#end venn function
######
######

#example data manipulation

setA = c( "THE1", "DOG", "RAN", "CROSS", "THE2", "ROAD" )
setB = c( "WHY", "DID", "THE1", "CHICKEN", "CROSS", "THE2", "ROAD" )
setC = c( "I", "CHASED", "MY", "DOG", "CROSS", "THE1", "ROAD" )

dframeA = data.frame( setA, key=1:length( setA ) )
dframeB = data.frame( setB, key=1:length( setB ) )
dframeC = data.frame( setC, key=1:length( setC ) )


######
#set operation for 2 variables
onlySetA = setdiff( setA, setB )

CommonToAB = intersect( setA, setB )

onlySetB = setdiff( setB, setA )

#examples of set operation for 3 variables,
#gets confusing doing all possible combinations - there's probably
#a better way to do this...
onlySetA_3var = intersect( setdiff(setA, setB), setdiff(setA, setC) )

CommonToABC = intersect( setA, intersect(setB, setC) )


#alternatively, set operations using dataframe instead of char vectors
intersect(dframeA$setA, dframeB$setB)

#end examples
######


######
#Build Venns and export relevant gene lists
setwd( "~/R/working" )

#basic gene list
write.table(CommonToAB, "CommonToAB.txt", sep = '\t')

#Venn with 2 character vectors
allVenns( setA, setB )

#Venn with columns from 2 dataframes
allVenns( dframeA$setA, dframeB$setB )

#Venn with 3 character vectors
allVenns( setA, setB, setC )

#Venn with columns from 2 dataframes
allVenns( dframeA$setA, dframeB$setB, dframeC$setC )
