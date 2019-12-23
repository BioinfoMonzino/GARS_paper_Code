rm(list = ls())
cat("\014")
graphics.off()

library(DaMiRseq)
library(RColorBrewer)
library(limma)
library(ggplot2)
library(GARS)
library(pheatmap)
library(caret)
library(randomForest)
library(e1071)
library(SummarizedExperiment)

### import raw data
#downloaded from MetaboLight: MTBLS24 (Zacharias et al 2103)
dataraw <- read.delim("dataset.txt")
classes <- read.delim("metadata.txt")

### filter sample (AKI2 and AKI3)
# 70 nonAKI
# 26 AKI-1

dataraw <- dataraw[,-c(28,33,99:106)]
classes <- droplevels(classes[-c(28,33,99:106),])
dataMatr.AKI <- as.matrix(t(scale(dataraw))))

## Generate learning set and independent test set
sample.per.class <- 20
set.seed(12345)
rnd.ind.sample.nonAKI <- sample(70,size = sample.per.class)
rnd.ind.sample.AKI <- sample(26,size = sample.per.class) + 70

rnd.ind.all <- c(rnd.ind.sample.nonAKI, rnd.ind.sample.AKI)

dataset <- dataMatr.AKI[rnd.ind.all,]
class.df <- classes[rnd.ind.all,]
testSet <- dataMatr.AKI[-rnd.ind.all,]
class.TS <- classes[-rnd.ind.all,]

# save data
save.image("./NormalizedDatasets.RData")
