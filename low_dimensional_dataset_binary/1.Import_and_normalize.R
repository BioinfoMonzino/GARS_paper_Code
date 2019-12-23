rm(list = ls())
cat("\014")
graphics.off()

library(DaMiRseq)
library(GARS)
library(RColorBrewer)
library(limma)
library(ggplot2)
library(MLSeq)

# load cervical dataset from MLSeq package
filepath <- system.file("extdata/cervical.txt", package = "MLSeq")
cervical <- read.table(filepath, header=TRUE)

# replace "wild-card" characters with "underscore"
rownames(cervical) <- gsub("*", "x", rownames(cervical), fixed = TRUE)
rownames(cervical) <- gsub("-", "_", rownames(cervical), fixed = TRUE)

# create the "class" vector
class_vector <- data.frame(gsub('[0-9]+', '', colnames(cervical)))
colnames(class_vector) <- "class"
rownames(class_vector) <- colnames(cervical)

# create a Summarized Experiment object
SE_obj <- DaMiR.makeSE(cervical, class_vector)

# filter and normalize the dataset by VST
datanorm <- DaMiR.normalization(SE_obj, th.cv = 100)
classes <- as.data.frame(colData(datanorm))
dataMatr <- as.matrix(t(assay(datanorm)))

## Generate learning set and independent test set
sample.per.class <- 25
set.seed(123)
rnd.ind.sample.N <- sample(29,size = sample.per.class)
rnd.ind.sample.T <- sample(29,size = sample.per.class) + 29
rnd.ind.all <- c(rnd.ind.sample.N, rnd.ind.sample.T)

dataset <- dataMatr[rnd.ind.all,]
class.df <- classes[rnd.ind.all,,drop=FALSE]
testSet <- dataMatr[-rnd.ind.all,]
class.TS <- classes[-rnd.ind.all,,drop=FALSE]

# save data
save.image("./NormalizedDatasets.RData")
