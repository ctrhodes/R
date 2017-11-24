# source("https://bioconductor.org/biocLite.R")
# biocLite("EDASeq")
# biocLite("DESeq2")
# biocLite("TCGAbiolinks")
# biocLite("SummarizedExperiment")
# biocLite("biomaRt")

rm(list = ls(all=TRUE))
require(TCGAbiolinks)
require(SummarizedExperiment)
require(DESeq2)
#require(biomaRt)
require(data.table)
require(tidyverse)

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

# gbm_hgnc_ensembl <- read_csv("~/R-scripts/TCGA/gbm_hgnc_ensembl.csv")
# gbm_hgnc_ensembl = gbm_hgnc_ensembl %>% 
#   select(-gene_id)
# head(gbm_hgnc_ensembl)

setwd("/home/chris/R-scripts/TCGA")
# PatientIDs_with_controls.csv has new RNA and mri
# caseTab <- read_csv("~/R-scripts/TCGA/PatientIDs_with_controls.csv")
# only primary samples
# need to filter out idh
caseTab <- read_csv("/home/chris/R-scripts/TCGA/input/TGCA_samples_Primary_Controls.csv")
# controls primary only and idhWT, possibly subset into groups
#caseTab <- read_csv("~/R-scripts/TCGA/TCGA_noIDH_TP_NT.csv")
head(caseTab)
samples = caseTab$cases

listSamples = samples
listSamples

query <- GDCquery(project = "TCGA-GBM", 
                  data.category = "Gene expression",
                  data.type = "Gene expression quantification",
                  experimental.strategy = "RNA-Seq",
                  platform = "Illumina HiSeq",
                  file.type = "results",
                  barcode = listSamples, 
                  legacy = TRUE)

# Download a list of barcodes with platform IlluminaHiSeq_RNASeqV2
GDCdownload(query)

# Prepare expression matrix with geneID in the rows and samples (barcode) in the columns
# rsem.genes.results as values
GBMRnaseqSE <- GDCprepare(query)

samplesTable = data.table(as.data.frame(colData(GBMRnaseqSE)), 
                          rownames = FALSE)

GBMMatrix <- assay(GBMRnaseqSE,"raw_count")

head(samplesTable)
idhWT = samplesTable %>% 
  filter(subtype_IDH.status != "Mutant")
head(idhWT)
idhWT$barcode

# # For gene expression if you need to see a boxplot correlation and AAIC plot to define outliers you can run
# BRCARnaseq_CorOutliers <- TCGAanalyze_Preprocessing(GBMMatrix)

# normalization of genes
dataNorm <- TCGAanalyze_Normalization(tabDF = GBMMatrix, geneInfo =  geneInfo)

# quantile filter of genes
dataFilt <- TCGAanalyze_Filtering(tabDF = dataNorm,
                                  method = "quantile", 
                                  qnt.cut =  0.25)

# selection of normal samples "NT"
samplesNT <- TCGAquery_SampleTypes(barcode = colnames(dataFilt),
                                   typesample = c("NT"))
samplesNT

# selection of tumor samples "TP"
samplesTP <- TCGAquery_SampleTypes(barcode = colnames(dataFilt), 
                                   typesample = c("TP"))
samplesTP
class(samplesTP)
samplesTP = intersect(samplesTP, idhWT$barcode)


# Diff.expr.analysis (DEA)
dataDEGs <- TCGAanalyze_DEA(mat1 = dataFilt[,samplesNT],
                            mat2 = dataFilt[,samplesTP],
                            Cond1type = "Normal",
                            Cond2type = "Tumor",
                            fdr.cut = 0.01 ,
                            logFC.cut = 1,
                            method = "glmLRT")

dataDEGsFiltLevel <- TCGAanalyze_LevelTab(dataDEGs,"Tumor","Normal",
                                          dataFilt[,samplesTP],dataFilt[,samplesNT])
dataDEGsFiltLevel

getwd()
write.csv(dataDEGsFiltLevel, file = "dataDEGsFiltLevel.csv", row.names = FALSE)



TGCA_samples_Primary_Controls = c(samplesTP, samplesNT)
TGCA_samples_Primary_Controls = as.data.frame(TGCA_samples_Primary_Controls)
colnames(TGCA_samples_Primary_Controls) <- "cases"
head(TGCA_samples_Primary_Controls)

write.csv(TGCA_samples_Primary_Controls, file = "TGCA_samples_Primary_Controls.csv", row.names = FALSE)




pca <- TCGAvisualize_PCA(dataFilt[,c(samplesNT, samplesTP)],dataDEGsFiltLevel, ntopgenes = 200, samplesNT, samplesTP)
pca$x



#################################
#################################
query <- GDCquery(project = c("TCGA-GBM"),
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification", 
                  workflow.type = "HTSeq - FPKM-UQ",
                  legacy = FALSE
)

head(getResults(query))
dim(getResults(query))

cases_expr = getResults(query, cols = c("data_type","cases")) %>% 
  arrange(cases) %>% 
  mutate(prefix = substr(cases, 1, 12)) %>% 
  distinct(prefix, .keep_all = TRUE)
head(cases_expr)

TCIA_TCGA_GBM_seriesInstanceUids1502381652704 <- read_csv("~/R-scripts/TCGA/TCIA-TCGA-GBM-seriesInstanceUids1502381652704.csv")
TCIA_TCGA_GBM_seriesInstanceUids1502381652704
tcia = TCIA_TCGA_GBM_seriesInstanceUids1502381652704 %>% 
  select(c(`Patient Id`, `Series Description`)) %>% 
  arrange(`Patient Id`) %>% 
  distinct(`Patient Id`, .keep_all = TRUE)
head(tcia)

case_list = inner_join(cases_expr, tcia, by=c("prefix" = "Patient Id") ) %>% 
  select(cases)
head(case_list)

write.csv(case_list, file = "PatientIDs.csv")

grades <- read_csv("~/R-scripts/TCGA/PatientIDs_with_controls.csv")
gradeDF = grades %>% 
  filter(!is.na(SVZ)) %>%
  select(prefix, group) 
# mutate(group = gsub("2", "1", group)) %>% 
# mutate(group = gsub("4", "3", group))
gradeDF
id = gradeDF$prefix
id

query2 <- GDCquery(project = c("TCGA-GBM"),
                   data.category = "Transcriptome Profiling",
                   data.type = "Gene Expression Quantification",
                   workflow.type = "HTSeq - FPKM-UQ",
                   barcode = c(case_list$cases)
)

GDCdownload(query2)
data <- GDCprepare(query2)


head(data)

sample = data.table(as.data.frame(colData(data)), 
                     rownames = FALSE)
head(samples)
idhWT = samples %>% 
  filter(subtype_IDH.status != "Mutant")
head(idhWT)

features = rowRanges(data)
head(features)

htseq = assay(data, "HTSeq - FPKM-UQ")
head(htseq)

#for heatmap only! altering colnames will make heatmap more readable, but likely
#unusable for additional summarizedExperiment calls

expr = htseq
colnames(expr) = substr(colnames(expr), 1, 12)
#expr = expr[ , -which(colnames(expr) %in% c("TCGA-14-0871"))]
colnames(expr)
head(expr)
#ensembls = rownames(expr)
#head(ensembls)

#keep expression idh-WT data
idhWT_use <- colnames(expr)[colnames(expr) %in% unique(as.character(idhWT$subtype_patient))]
expr <- expr[, idhWT_use]

#keep expression tcia data
names_use <- colnames(expr)[colnames(expr) %in% id]
expr <- expr[, names_use]
dim(expr)
head(expr)

expr = as.data.frame(expr) %>% 
  tibble::rownames_to_column()
names(expr)[1] = "ensembl_gene_id"
head(expr)

expr_gbm = inner_join(expr, gbm_hgnc_ensembl, by = c("ensembl_gene_id", "ensembl_gene_id"))
head(expr_gbm)

expr_gbm = as.data.frame(expr_gbm) %>% 
  distinct(ensembl_gene_id, .keep_all = TRUE) %>% 
  tibble::column_to_rownames(var = "ensembl_gene_id")
expr_gbm = as.matrix(expr_gbm)
dim(expr_gbm)

expr_gbm[expr_gbm == 0] <- NA
expr_gbm = expr_gbm[complete.cases(expr_gbm),]
is.na(expr_gbm)
dim(expr_gbm)
head(expr_gbm)
reg_gbm = t(expr_gbm)
reg_gbm[1:5,1:5]

variance = apply(expr_gbm, 1, var)
expr_gbm = cbind(expr_gbm, variance)
head(expr_gbm)
dim(expr_gbm)

# for pca using only tcga-gbm (no gbm1 or gbm2 in vivo case) and fpkm-uq,
# select genes with: 1.0 > variance > 0.90 (above 90th percentile)
# for comparing tcga-gbm and gbm1, gbm2 using fpkm (not fpkm-uq),
# use variance: 0.95 > variance > 0.4 (possibly 0.9)
expr2 = as.data.frame(expr_gbm) %>% 
  tibble::rownames_to_column() %>% 
  filter(variance <= quantile(variance, c(1.0)) & 
           variance > quantile(variance, c(.90))) %>% 
  tibble::column_to_rownames() %>%
  t()
head(expr2[1:5,1:5])
dim(expr2)
expr2 = expr2[!rownames(expr2) %in% "variance", ]
dim(expr2)

expr2[1:5, 1:5]
tail(expr2[,1:5])

expr2 = log2(expr2)
which(apply(expr2, 2, var)==0)


group = subset(gradeDF, prefix %in% names_use)
head(group)
dim(group)

gbm12 <- data.frame(
  prefix = c("GBM-1", "GBM-2"),
  group = c(1, 2)
)
head(gbm12)

group = rbind(group, gbm12)
head(group)
tail(group)

pr.out=prcomp(expr2, center = TRUE, scale. = TRUE)
head(pr.out$x[,1:2])
pr_indiv = tbl_df(pr.out$x[,1:2]) %>% 
  mutate(GBM = group$group)
head(pr_indiv)

pr_med = tbl_df(pr.out$x[,1:2]) %>% 
  mutate(GBM = group$group) %>% 
  group_by(GBM) %>% 
  summarize(PC1 = median(PC1), PC2 = median(PC2)) %>% 
  mutate(GBM = GBM + 4) %>% 
  select(PC1, PC2, GBM)
head(pr_med)

pr_full = rbind(pr_indiv, pr_med)
head(pr_full)
tail(pr_full)
dim(pr_full)
pr_full_sub = pr_full[58:63,]
pr_full_sub

qplot(PC1, PC2, data = data.frame(pr_full_sub), 
      xlab = "PC1" , ylab = "PC2", colour = as.factor(GBM))

dim(pr.out$x)
pr.tcga = pr.out$x[1:57,]
group.tcga = group$group[1:57]
pr.tcga.group = cbind(pr.tcga, group.tcga)
qplot(PC1, PC2, data = data.frame(pr.tcga.group), 
      xlab = "PC1" , ylab = "PC2", colour = as.factor(group.tcga))



d <- dist(expr2) # euclidean distances between the rows
fit <- cmdscale(d,eig=TRUE, k=2) # k is the number of dim
#fit <- isoMDS(d, k=2)
fit # view results

# plot solution 
x <- fit$points[,1]
y <- fit$points[,2]

qplot(x, y, xlab = "Coor1" , ylab = "Coor2", colour = as.factor(group$group))


#multiple models of log regression - just like titanic

