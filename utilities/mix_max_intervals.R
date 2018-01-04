#Example coordinates are tuple
# with "START" coordinate
#end with "END" coordinate
#coordinates in example are full-open
#if need half-open as used in UCSC genome browser, subtract 1 from "START" integer

#example = c("START", "END")
yourCoords = c(1, 2)
yourCoords

refGenome = c(3, 4)
refGenome

#find minimal distance between 2 genomic elements
#similar to BEDOPS/ BEDTOOLS closest feature tool
min(abs(outer(yourCoords, refGenome, "-")))

#find max distance from "START" of 1 element to "END" of second
max(abs(outer(yourCoords, refGenome, "-")))
