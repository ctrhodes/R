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


primary_query <- GDCquery(project = "TCGA-GBM",
                      data.category = "Transcriptome Profiling",
                      data.type = "Gene Expression Quantification",
                      workflow.type = "HTSeq - FPKM-UQ",
                      barcode = primarySamples$cases
)

tail(getResults(primary_query))
dim(getResults(primary_query))


# Download a list of barcodes with platform IlluminaHiSeq_RNASeqV2
GDCdownload(primary_query)


# Prepare expression matrix with geneID in the rows and samples (barcode) in the columns
# rsem.genes.results as values
prim_exp38_SE <- GDCprepare(primary_query)

prim_samplesTable = data.table(as.data.frame(colData(prim_exp38_SE)), 
                          rownames = FALSE)
head(prim_samplesTable)
dim(prim_samplesTable)

prim_samplesTable = cbind(prim_samplesTable, primarySamples$group)
names(prim_samplesTable)[names(prim_samplesTable) == 'V2'] <- 'group'
tail(prim_samplesTable)
dim(prim_samplesTable)

prim_Matrix <- assay(prim_exp38_SE,"HTSeq - FPKM-UQ")
head(prim_Matrix)
dim(prim_Matrix)

tail(prim_samplesTable)

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
  mutate(age_decades = findInterval(.$age_at_diagnosis/365, 
                                c(0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100))) %>% 
  select(-rownames)
tail(prim_idhWT)
dim(prim_idhWT)

# TCGAanalyze_survival(prim_idhWT,
#                      "group",
#                      main = "TCGA Set\n GBM",height = 10, width=10)
# 
# x = matrix(c(0,1,2,3,4,5,6,0,8), nrow = 3)
# x = as.data.frame(x)
# x
# x[x==0] = 1
# x

prim_mat_subset = as.data.frame(prim_Matrix) %>% 
  select(prim_idhWT$barcode)
#prim_mat_subset[prim_mat_subset==0] = 0.01

# try below to truly remove zeros
# prim_mat_subset[prim_mat_subset==0] <- NA
# prim_mat_subset = prim_mat_subset[complete.cases(prim_mat_subset), ]

dim(prim_mat_subset)
prim_mat_subset[1:5,1:5]

shapiro.test(as.numeric(prim_mat_subset[100,]))
qqnorm(prim_mat_subset[100,])
qqline(prim_mat_subset[100,])


prim_preSplits = t(prim_mat_subset)
dim(prim_preSplits)
prim_preSplits[57:62,51717:51722]


#####
##### to assess LASSO without "healthy" controls, drop last 5 rows here
#####

healthy
prim_preSplits = prim_preSplits[!rownames(prim_preSplits) %in% healthy, ]
prim_preSplits[51:57,51717:51722]

prim_idhWT$group
hlth_grp = which(prim_idhWT$group == 0)
prim_idhWT_grp = prim_idhWT$group[-hlth_grp]
prim_idhWT_grp

####
#### remove columns with no variance or remove later with caret preProcess()
####
# nzv <- nearZeroVar(prim_preSplits, saveMetrics= TRUE)
# nzv[nzv$nzv,][1:10,]
# 
# dim(prim_preSplits)
# nzv <- nearZeroVar(prim_preSplits)
# filteredDescr <- prim_preSplits[, -nzv]
# dim(filteredDescr)


set.seed(0)
# without controls
trainIndex <- createDataPartition(as.factor(prim_idhWT_grp), p = .75, list = FALSE)
head(trainIndex)

grp_Train <- prim_idhWT_grp[trainIndex]
grp_Train
length(grp_Train)
grp_Test  <- prim_idhWT_grp[ -trainIndex]
grp_Test
length(grp_Test)

# #for inclusions with controls
# trainIndex <- createDataPartition(as.factor(prim_idhWT$group), p = .75, list = FALSE)
# head(trainIndex)
# 
# grp_Train <- prim_idhWT$group[trainIndex]
# grp_Train
# length(grp_Train)
# grp_Test  <- prim_idhWT$group[ -trainIndex]
# grp_Test
# length(grp_Test)

prim_Train = prim_preSplits[ trainIndex,]
prim_Train[1:5,1:5]
class(prim_Train)
dim(prim_Train)
prim_Test = prim_preSplits[-trainIndex,]
prim_Test[1:5,1:5]
class(prim_Test)
dim(prim_Test)


# only select one transformation to get data normal. Either log2 OR BoxCox


preProcValues <- preProcess(prim_Train, method = c("center", "scale", "BoxCox", "nzv"))
preProcValues

train_BC <- predict(preProcValues, prim_Train)
train_BC[1:5,1:5]
hist(as.numeric(train_BC[,100]))
shapiro.test(as.numeric(train_BC[,100]))
qqnorm(train_BC[,100])
qqline(train_BC[,100])
which(is.na(train_BC), arr.ind=TRUE)

test_BC <- predict(preProcValues, prim_Test)
test_BC[1:5,1:5]
hist(as.numeric(test_BC[,100]))
shapiro.test(as.numeric(test_BC[,100]))
qqnorm(test_BC[,100])
qqline(test_BC[,100])
which(is.na(test_BC), arr.ind=TRUE)


is.infinite.data.frame <- function(obj){
  sapply(obj,FUN = function(x) all(is.infinite(x)))
}
is.infinite.data.frame(trainBC)

#W close to 1 and p > 0.05 is normal distribution
#W lower value and p <= 0.05 is non-normal distribution

train5_x = rbind(train_BC, train_BC, train_BC, train_BC, train_BC)
class(train5_x)
train5_y = c(grp_Train, grp_Train, grp_Train, grp_Train, grp_Train)
class(train5_y)

test5_x = rbind(test_BC, test_BC, test_BC, test_BC, test_BC)
class(test5_x)
test5_y = c(grp_Test, grp_Test, grp_Test, grp_Test, grp_Test)
class(test5_y)

dim(train5_x)
length(train5_y)

dim(test5_x)
length(test5_y)


#alpha = 0: ridge, alpha = 1: lasso
cvfit = cv.glmnet(x = train5_x, y = as.matrix(train5_y),
                  nfolds = 5,
                  alpha = 1,
                  family="multinomial", type.measure = "class")

summary(cvfit)
plot(cvfit)

cvfit$lambda.min

head(coef(cvfit, s = "lambda.min")[[1]])
#coerse dgCMatrix to DF
coef_df = as.data.frame(as.matrix(coef(cvfit, s = "lambda.min")[[1]]))

head(coef_df)

nz_coefs = as.vector(which(coef_df$`1` != 0))
nz_coefs
row.names(coef_df)[nz_coefs]
coef_df[nz_coefs,]


cvpred = predict(cvfit, newx = test5_x, s = "lambda.min",
                 type = "class")

head(cvpred)

table(as.numeric(cvpred), test5_y)
