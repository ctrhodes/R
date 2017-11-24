# source("https://bioconductor.org/biocLite.R")
# biocLite("DESeq2")
# biocLite("EDASeq")
# biocLite("TCGAbiolinks")
# biocLite("SummarizedExperiment")

pkg = names(sessionInfo()$otherPkgs)
pkgs = paste('package:', pkg, sep = "")
lapply(pkgs, detach, character.only = TRUE, unload = TRUE)

rm(list = ls(all=TRUE))

require(data.table)
require(tidyverse)
require(TCGAbiolinks)
require(SummarizedExperiment)
require(DESeq2)
require(EDASeq)
require(glmnet)
require(caret)
require(e1071)

require(data.table)



#MASS causes dplyr::select to crash
#either invoke MASS options on the fly
#or call dplyr::select...
#already have a lot of dplyr code, so invoke MASS on the fly...
#require(MASS)

if (!require("gplots")) {
  install.packages("gplots", dependencies = TRUE)
  library(gplots)
}
if (!require("RColorBrewer")) {
  install.packages
  library(RColorBrewer)
}

setwd("C:\\Users\\Chris\\Documents\\TCGA\\working")

# controls primary only and idhWT, possibly subset into groups
#full case IDs for GDCquery
TCGAbiolinks:::getGDCprojects()$project_id

TCGAbiolinks:::getProjectSummary("TCGA-GBM")

all_query <- GDCquery(project = "TCGA-GBM",
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification",
                  workflow.type = "HTSeq - FPKM-UQ"
                  )

tail(getResults(all_query))
dim(getResults(all_query))

tcga_query = getResults(all_query)
head(tcga_query)


tcga_query = tcga_query %>% 
  mutate(prefix = gsub("(TCGA-[0-9]{2,4}-[A-Z0-9]{2,4}).*", "\\1", cases))
head(tcga_query)
class(tcga_query$prefix)

caseTab <- read_csv("C:\\Users\\Chris\\Documents\\TCGA\\input\\PatientIDs_hg38-manual-noNA.csv")
head(caseTab)

combs = inner_join(tcga_query, caseTab, by = c("prefix" = "cases")) %>% 
  arrange(cases)
head(combs)
sort(combs$cases)

primary_38 = combs %>% 
  filter(tissue.definition == "Primary solid Tumor") %>% 
  select(cases, SVZ, CTX, group)
head(primary_38)
dim(primary_38)

recurrent_38 = combs %>% 
  filter(tissue.definition == "Recurrent Solid Tumor") %>% 
  select(cases, SVZ, CTX, group)
head(recurrent_38)
dim(recurrent_38)


contTab <- read_csv("C:\\Users\\Chris\\Documents\\TCGA\\input\\TGCA_Controls.csv")

head(contTab)
contTab = contTab %>% 
  select(cases, SVZ, CTX, group)
head(contTab)

primarySamples = rbind(primary_38, contTab)
head(primarySamples)
tail(primarySamples)
dim(primarySamples)

recurrentSamples = rbind(recurrent_38, contTab)
head(recurrentSamples)
tail(recurrentSamples)
dim(recurrentSamples)


primary19_liftedfrom38_query <- GDCquery(project = "TCGA-GBM", 
                  data.category = "Gene expression",
                  data.type = "Gene expression quantification",
                  experimental.strategy = "RNA-Seq",
                  platform = "Illumina HiSeq",
                  file.type = "results",
                  barcode = primarySamples$cases, 
                  legacy = TRUE)

tail(getResults(primary19_liftedfrom38_query))
dim(getResults(primary19_liftedfrom38_query))

hg19_results = getResults(primary19_liftedfrom38_query)
hg19_results$cases

primarySamples = primarySamples %>%
  filter(grepl(paste(hg19_results$cases, collapse="|"), cases))
dim(primarySamples)

primary19_query <- GDCquery(project = "TCGA-GBM",
                            data.category = "Gene expression",
                            data.type = "Gene expression quantification",
                            experimental.strategy = "RNA-Seq",
                            platform = "Illumina HiSeq",
                            file.type = "results",
                            barcode = hg19_results$cases,
                            legacy = TRUE)

# primary_query <- GDCquery(project = "TCGA-GBM",
#                       data.category = "Transcriptome Profiling",
#                       data.type = "Gene Expression Quantification",
#                       workflow.type = "HTSeq - Counts",
#                       barcode = hg19_results$cases
# )

# Download a list of barcodes with platform IlluminaHiSeq_RNASeqV2
GDCdownload(primary19_query)


# Prepare expression matrix with geneID in the rows and samples (barcode) in the columns
# rsem.genes.results as values
prim_exp19_SE <- GDCprepare(primary19_query)

prim_samplesTable = data.table(as.data.frame(colData(prim_exp19_SE)), 
                          rownames = FALSE)
head(prim_samplesTable)
dim(prim_samplesTable)

prim_samplesTable = cbind(prim_samplesTable, primarySamples$group)
names(prim_samplesTable)[names(prim_samplesTable) == 'V2'] <- 'group'
tail(prim_samplesTable)
dim(prim_samplesTable)

prim_Matrix <- assay(prim_exp19_SE,"raw_count")
head(prim_Matrix)
dim(prim_Matrix)

head(prim_samplesTable)
#prim_samplesTable[60:80,]
which(prim_samplesTable$barcode == "TCGA-02-0047-01A-01R-1849-01")

healthy = c(
  "TCGA-06-0675-11A-32R-A36H-07",
  "TCGA-06-0678-11A-32R-A36H-07",
  "TCGA-06-0680-11A-32R-A36H-07",
  "TCGA-06-0681-11A-41R-A36H-07",
  "TCGA-06-AABW-11A-31R-A36H-07"
)
healthy

for (i in healthy){
  prim_samplesTable[prim_samplesTable$barcode == i, "subtype_IDH.status"] <- "WT"
}
tail(prim_samplesTable$subtype_IDH.status)

#Another way using the dplyr:
#df %>% mutate(Name = ifelse(State == "WI" & Name == "John_Smith", "John_Smith1", Name))


prim_idhWT = prim_samplesTable %>% 
  filter(subtype_IDH.status != "Mutant") %>% 
  select(-rownames)

head(prim_idhWT)
tail(prim_idhWT)
dim(prim_idhWT)
prim_idhWT[1:5,1:5]

dim(prim_Matrix)
prim_Matrix[1:5,1:5]
prim_Matrix = subset(prim_Matrix, select = prim_idhWT$barcode)
prim_Matrix[1:5,1:5]



prim_1 = prim_idhWT %>% 
  filter(group %in% c(1,0))
head(prim_1)
tail(prim_1)
dim(prim_1)

prim_Matrix[1:5,1:5]
prim_Matrix_1 = subset(prim_Matrix, select = prim_1$barcode)
dim(prim_Matrix_1)

prim_2 = prim_idhWT %>% 
  filter(group %in% c(2,0))
tail(prim_2)
dim(prim_2)

prim_Matrix_2 = subset(prim_Matrix, select = prim_2$barcode)
dim(prim_Matrix_2)

prim_3 = prim_idhWT %>% 
  filter(group %in% c(3,0))
head(prim_3)
dim(prim_3)

prim_Matrix_3 = subset(prim_Matrix, select = prim_3$barcode)
dim(prim_Matrix_3)

prim_4 = prim_idhWT %>% 
  filter(group %in% c(4,0))
tail(prim_4)
dim(prim_4)

prim_Matrix_4 = subset(prim_Matrix, select = prim_4$barcode)
dim(prim_Matrix_4)

# # For gene expression if you need to see a boxplot correlation and AAIC plot to define outliers you can run
#dataPrep <- TCGAanalyze_Preprocessing(object = prim_Matrix, datatype = "HTSeq - Counts") 

# normalization of genes
dataNorm <- TCGAanalyze_Normalization(tabDF = prim_Matrix, geneInfo =  geneInfo)
dataNorm_1 <- TCGAanalyze_Normalization(tabDF = prim_Matrix_1, geneInfo =  geneInfo)
dataNorm_2 <- TCGAanalyze_Normalization(tabDF = prim_Matrix_2, geneInfo =  geneInfo)
dataNorm_3 <- TCGAanalyze_Normalization(tabDF = prim_Matrix_3, geneInfo =  geneInfo)
dataNorm_4 <- TCGAanalyze_Normalization(tabDF = prim_Matrix_4, geneInfo =  geneInfo)


boxplot(dataNorm, outline = FALSE)
boxplot(dataNorm_1, outline = FALSE)
boxplot(dataNorm_2, outline = FALSE)
boxplot(dataNorm_3, outline = FALSE)
boxplot(dataNorm_4, outline = FALSE)

# quantile filter of genes
dataFilt <- TCGAanalyze_Filtering(tabDF = dataNorm, method = "quantile", qnt.cut =  0.25)
dataFilt_1 <- TCGAanalyze_Filtering(tabDF = dataNorm_1, method = "quantile", qnt.cut =  0.25)
dataFilt_2 <- TCGAanalyze_Filtering(tabDF = dataNorm_2, method = "quantile", qnt.cut =  0.25)
dataFilt_3 <- TCGAanalyze_Filtering(tabDF = dataNorm_3, method = "quantile", qnt.cut =  0.25)
dataFilt_4 <- TCGAanalyze_Filtering(tabDF = dataNorm_4, method = "quantile", qnt.cut =  0.25)


# selection of normal samples "NT"
samplesNT <- TCGAquery_SampleTypes(barcode = colnames(dataFilt), typesample = c("NT"))
samplesNT
samplesNT_1 <- TCGAquery_SampleTypes(barcode = colnames(dataFilt_1), typesample = c("NT"))
samplesNT_1
samplesNT_2 <- TCGAquery_SampleTypes(barcode = colnames(dataFilt_2), typesample = c("NT"))
samplesNT_2
samplesNT_3 <- TCGAquery_SampleTypes(barcode = colnames(dataFilt_3), typesample = c("NT"))
samplesNT_3
samplesNT_4 <- TCGAquery_SampleTypes(barcode = colnames(dataFilt_4), typesample = c("NT"))
samplesNT_4

# selection of tumor samples "TP"
samplesTP <- TCGAquery_SampleTypes(barcode = colnames(dataFilt), typesample = c("TP"))
samplesTP
samplesTP_1 <- TCGAquery_SampleTypes(barcode = colnames(dataFilt_1), typesample = c("TP"))
samplesTP_1
samplesTP_2 <- TCGAquery_SampleTypes(barcode = colnames(dataFilt_2), typesample = c("TP"))
samplesTP_2
samplesTP_3 <- TCGAquery_SampleTypes(barcode = colnames(dataFilt_3), typesample = c("TP"))
samplesTP_3
samplesTP_4 <- TCGAquery_SampleTypes(barcode = colnames(dataFilt_4), typesample = c("TP"))
samplesTP_4



# Diff.expr.analysis (DEA)
dataDEGs <- TCGAanalyze_DEA(mat1 = dataFilt[,samplesNT],
                            mat2 = dataFilt[,samplesTP],
                            Cond1type = "Normal",
                            Cond2type = "Tumor",
                            fdr.cut = 0.01 ,
                            logFC.cut = 1,
                            method = "glmLRT")

head(dataDEGs)

dataDEGs_1 <- TCGAanalyze_DEA(mat1 = dataFilt_1[,samplesNT_1],
                            mat2 = dataFilt_1[,samplesTP_1],
                            Cond1type = "Normal",
                            Cond2type = "Tumor",
                            fdr.cut = 0.01 ,
                            logFC.cut = 1,
                            method = "glmLRT")

head(dataDEGs_1)

dataDEGs_2 <- TCGAanalyze_DEA(mat1 = dataFilt_2[,samplesNT_2],
                              mat2 = dataFilt_2[,samplesTP_2],
                              Cond1type = "Normal",
                              Cond2type = "Tumor",
                              fdr.cut = 0.01 ,
                              logFC.cut = 1,
                              method = "glmLRT")

head(dataDEGs_2)

dataDEGs_3 <- TCGAanalyze_DEA(mat1 = dataFilt_3[,samplesNT_3],
                              mat2 = dataFilt_3[,samplesTP_3],
                              Cond1type = "Normal",
                              Cond2type = "Tumor",
                              fdr.cut = 0.01 ,
                              logFC.cut = 1,
                              method = "glmLRT")

head(dataDEGs_3)

dataDEGs_4 <- TCGAanalyze_DEA(mat1 = dataFilt_4[,samplesNT_4],
                              mat2 = dataFilt_4[,samplesTP_4],
                              Cond1type = "Normal",
                              Cond2type = "Tumor",
                              fdr.cut = 0.01 ,
                              logFC.cut = 1,
                              method = "glmLRT")

head(dataDEGs_4)

dataDEGsFiltLevel <- TCGAanalyze_LevelTab(dataDEGs,"Tumor","Normal",
                                          dataFilt[,samplesTP], dataFilt[,samplesNT])
head(dataDEGsFiltLevel)

dataDEGsFiltLevel_1 <- TCGAanalyze_LevelTab(dataDEGs_1,"Tumor","Normal",
                                          dataFilt[,samplesTP_1], dataFilt[,samplesNT_1])
head(dataDEGsFiltLevel_1)

dataDEGsFiltLevel_2 <- TCGAanalyze_LevelTab(dataDEGs_2,"Tumor","Normal",
                                          dataFilt[,samplesTP_2], dataFilt[,samplesNT_2])
head(dataDEGsFiltLevel_2)

dataDEGsFiltLevel_3 <- TCGAanalyze_LevelTab(dataDEGs_3,"Tumor","Normal",
                                          dataFilt[,samplesTP_3], dataFilt[,samplesNT_3])
head(dataDEGsFiltLevel_3)

dataDEGsFiltLevel_4 <- TCGAanalyze_LevelTab(dataDEGs_4,"Tumor","Normal",
                                          dataFilt[,samplesTP_4], dataFilt[,samplesNT_4])
head(dataDEGsFiltLevel_4)

getwd()
setwd("C:/Users/Chris/Documents/TCGA/output")

write.csv(dataDEGsFiltLevel, file = "primaryTumor_hg19_DEGs.csv", row.names = FALSE)
write.csv(dataDEGsFiltLevel_1, file = "group_1_hg19_DEGs.csv", row.names = FALSE)
write.csv(dataDEGsFiltLevel_2, file = "group_2_hg19_DEGs.csv", row.names = FALSE)
write.csv(dataDEGsFiltLevel_3, file = "group_3_hg19_DEGs.csv", row.names = FALSE)
write.csv(dataDEGsFiltLevel_4, file = "group_4_hg19_DEGs.csv", row.names = FALSE)

common = intersect(
  intersect(
    intersect(dataDEGsFiltLevel_1$mRNA, dataDEGsFiltLevel_2$mRNA), 
    dataDEGsFiltLevel_3$mRNA), 
  dataDEGsFiltLevel_4$mRNA)
length(common)

uniq_1 = setdiff(setdiff(setdiff(setdiff(dataDEGsFiltLevel_1$mRNA, common), 
                         dataDEGsFiltLevel_2$mRNA),
                 dataDEGsFiltLevel_3$mRNA),
                 dataDEGsFiltLevel_4$mRNA)
uniq_1

uniq_2 = setdiff(setdiff(setdiff(setdiff(dataDEGsFiltLevel_2$mRNA, common), 
                                 dataDEGsFiltLevel_3$mRNA),
                         dataDEGsFiltLevel_4$mRNA),
                 dataDEGsFiltLevel_1$mRNA)
uniq_2

uniq_3 = setdiff(setdiff(setdiff(setdiff(dataDEGsFiltLevel_3$mRNA, common), 
                                 dataDEGsFiltLevel_4$mRNA),
                         dataDEGsFiltLevel_1$mRNA),
                 dataDEGsFiltLevel_2$mRNA)
uniq_3

uniq_4 = setdiff(setdiff(setdiff(setdiff(dataDEGsFiltLevel_4$mRNA, common), 
                                 dataDEGsFiltLevel_1$mRNA),
                         dataDEGsFiltLevel_2$mRNA),
                 dataDEGsFiltLevel_3$mRNA)
uniq_4

write.csv(common, file = "common_to_groups_1-4_hg19_DEGs.csv", row.names = FALSE)
write.csv(uniq_1, file = "distinct_group_1_hg19_DEGs.csv", row.names = FALSE)
write.csv(uniq_2, file = "distinct_group_2_hg19_DEGs.csv", row.names = FALSE)
write.csv(uniq_3, file = "distinct_group_3_hg19_DEGs.csv", row.names = FALSE)
write.csv(uniq_4, file = "distinct_group_4_hg19_DEGs.csv", row.names = FALSE)
