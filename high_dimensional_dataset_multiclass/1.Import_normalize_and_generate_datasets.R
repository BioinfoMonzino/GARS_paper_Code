library(DaMiRseq)
library(RColorBrewer)
library(limma)
library(ggplot2)
library(GARS)
library(pheatmap)
library(caret)
library(randomForest)
library(e1071)


# import file from GTEx
# download raw data from: https://storage.googleapis.com/gtex_analysis_v6p/rna_seq_data/GTEx_Analysis_v6p_RNA-seq_RNA-SeQCv1.1.8_gene_reads.gct.gz
# and save in "data"
# selecte only sample from the following brain regions:
#
# "Spinal Cord"
# "Nucleus accumbens"
# "Putamen"
# "Caudate"
# "Substantia Nigra"
# "Hypothalamus"
# "Hippocampus"
# "Amygdala"
# "Anterior Cingulate Cortex"
# "Frontal Cortex"
# "Cortex" 

dati <- read.delim("All.15.brain.tissues.txt")
#dati <- dati[,-1]
class_info <- read.delim("Covariates.txt")
dam.obj <- DaMiR.makeSE(dati,class_info)

# normalization and filtering
dam.norm <- DaMiR.normalization(dam.obj,
                                minCounts = 10,
                                fSample = 0.7)

class_info_all <- as.data.frame(colData(dam.norm))
dataMatr.tissues_all <- as.matrix(t(assay(dam.norm)))


## Leraning Set and tIndependent Test Set (11 classes)
sample.per.class <- 50
set.seed(123)
rnd.ind.sample.AMY <- sample(72, size = sample.per.class)
rnd.ind.sample.ACC <- sample(84, size = sample.per.class) + 72
rnd.ind.sample.CAU <- sample(117, size = sample.per.class) + 84 + 72
rnd.ind.sample.COR <- sample(114, size = sample.per.class) + 117 + 84 + 72
rnd.ind.sample.FRC <- sample(108, size = sample.per.class) + 114 + 117 + 84 + 72
rnd.ind.sample.HIP <- sample(94, size = sample.per.class) + 108 + 114 + 117 + 84 + 72
rnd.ind.sample.HYP <- sample(96, size = sample.per.class) + 94 + 108 + 114 + 117 + 84 + 72
rnd.ind.sample.NAC <- sample(113, size = sample.per.class) + 96 + 94 + 108 + 114 + 117 + 84 + 72
rnd.ind.sample.PUT <- sample(97, size = sample.per.class) + 113 + 96 + 94 + 108 + 114 + 117 + 84 + 72
rnd.ind.sample.SPC <- sample(71, size = sample.per.class) + 97 + 113 + 96 + 94 + 108 + 114 + 117 + 84 + 72
rnd.ind.sample.SNI <- sample(63, size = sample.per.class) + 71 + 97 + 113 + 96 + 94 + 108 + 114 + 117 + 84 + 72


rnd.ind.all <- c(rnd.ind.sample.AMY,
                 rnd.ind.sample.ACC,
                 rnd.ind.sample.CAU,
                 rnd.ind.sample.COR,
                 rnd.ind.sample.FRC,
                 rnd.ind.sample.HIP,
                 rnd.ind.sample.HYP,
                 rnd.ind.sample.NAC,
                 rnd.ind.sample.PUT,
                 rnd.ind.sample.SPC,
                 rnd.ind.sample.SNI)

dataset <- dataMatr.tissues[rnd.ind.all,]
class.df <- class_info[rnd.ind.all,]
testSet <- dataMatr.tissues[-rnd.ind.all,]
class.TS <- class_info[-rnd.ind.all,]

###########################
# Generate Sub datasets
tissue2filt <- c("Spinal Cord", "Nucleus accumbens") # 9 tissues
tissue2filt <- c("Spinal Cord", "Nucleus accumbens","Putamen","Caudate") # 7 tissues
tissue2filt <- c("Spinal Cord", "Nucleus accumbens","Putamen","Caudate","Substantia Nigra","Hypothalamus") # 5 tissues
tissue2filt <- c("Spinal Cord", "Nucleus accumbens","Putamen","Caudate","Substantia Nigra","Hypothalamus","Hippocampus","Amygdala") # 3 tissues

index2filt_tr <- which(class.df$class %in% tissue2filt)
index2filt_ts <- which(class.TS$class %in% tissue2filt)

dataset <- dataset[-index2filt_tr,]
class.df <- droplevels(class.df[-index2filt_tr,])
testSet <- testSet[-index2filt_ts,]
class.TS <- droplevels(class.TS[-index2filt_ts,])

