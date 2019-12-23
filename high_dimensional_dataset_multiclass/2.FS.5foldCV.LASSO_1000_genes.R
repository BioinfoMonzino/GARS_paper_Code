library(DaMiRseq)
library(RColorBrewer)
library(limma)
library(ggplot2)
library(GARS)
library(pheatmap)
library(caret)
library(randomForest)
library(e1071)
library(microbenchmark)
library(ROCR)
library(ineq)

## Load data
# 11 tissues
load("./NormalizedData_Brain_11_tissues.RData")
cv.sample.array <- c(1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,3,3,4,4,4,4,4,4,4,4,4,4,5,5,5,5,5,5,5,5,5,5,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,3,3,4,4,4,4,4,4,4,4,4,4,5,5,5,5,5,5,5,5,5,5,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,3,3,4,4,4,4,4,4,4,4,4,4,5,5,5,5,5,5,5,5,5,5,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,3,3,4,4,4,4,4,4,4,4,4,4,5,5,5,5,5,5,5,5,5,5,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,3,3,4,4,4,4,4,4,4,4,4,4,5,5,5,5,5,5,5,5,5,5,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,3,3,4,4,4,4,4,4,4,4,4,4,5,5,5,5,5,5,5,5,5,5,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,3,3,4,4,4,4,4,4,4,4,4,4,5,5,5,5,5,5,5,5,5,5,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,3,3,4,4,4,4,4,4,4,4,4,4,5,5,5,5,5,5,5,5,5,5,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,3,3,4,4,4,4,4,4,4,4,4,4,5,5,5,5,5,5,5,5,5,5,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,3,3,4,4,4,4,4,4,4,4,4,4,5,5,5,5,5,5,5,5,5,5,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,3,3,4,4,4,4,4,4,4,4,4,4,5,5,5,5,5,5,5,5,5,5)

# 9 tissues
load("./NormalizedData_Brain_9_tissues.RData")
cv.sample.array <- c(1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,3,3,4,4,4,4,4,4,4,4,4,4,5,5,5,5,5,5,5,5,5,5,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,3,3,4,4,4,4,4,4,4,4,4,4,5,5,5,5,5,5,5,5,5,5,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,3,3,4,4,4,4,4,4,4,4,4,4,5,5,5,5,5,5,5,5,5,5,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,3,3,4,4,4,4,4,4,4,4,4,4,5,5,5,5,5,5,5,5,5,5,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,3,3,4,4,4,4,4,4,4,4,4,4,5,5,5,5,5,5,5,5,5,5,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,3,3,4,4,4,4,4,4,4,4,4,4,5,5,5,5,5,5,5,5,5,5,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,3,3,4,4,4,4,4,4,4,4,4,4,5,5,5,5,5,5,5,5,5,5,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,3,3,4,4,4,4,4,4,4,4,4,4,5,5,5,5,5,5,5,5,5,5,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,3,3,4,4,4,4,4,4,4,4,4,4,5,5,5,5,5,5,5,5,5,5)

# 7 tissues
load("./NormalizedData_Brain_7_tissues.RData")
cv.sample.array <- c(1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,3,3,4,4,4,4,4,4,4,4,4,4,5,5,5,5,5,5,5,5,5,5,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,3,3,4,4,4,4,4,4,4,4,4,4,5,5,5,5,5,5,5,5,5,5,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,3,3,4,4,4,4,4,4,4,4,4,4,5,5,5,5,5,5,5,5,5,5,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,3,3,4,4,4,4,4,4,4,4,4,4,5,5,5,5,5,5,5,5,5,5,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,3,3,4,4,4,4,4,4,4,4,4,4,5,5,5,5,5,5,5,5,5,5,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,3,3,4,4,4,4,4,4,4,4,4,4,5,5,5,5,5,5,5,5,5,5,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,3,3,4,4,4,4,4,4,4,4,4,4,5,5,5,5,5,5,5,5,5,5)

# 5 tissues
load("./NormalizedData_Brain_5_tissues.RData")
cv.sample.array <- c(1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,3,3,4,4,4,4,4,4,4,4,4,4,5,5,5,5,5,5,5,5,5,5,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,3,3,4,4,4,4,4,4,4,4,4,4,5,5,5,5,5,5,5,5,5,5,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,3,3,4,4,4,4,4,4,4,4,4,4,5,5,5,5,5,5,5,5,5,5,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,3,3,4,4,4,4,4,4,4,4,4,4,5,5,5,5,5,5,5,5,5,5,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,3,3,4,4,4,4,4,4,4,4,4,4,5,5,5,5,5,5,5,5,5,5)

# 3 tissues
load("./NormalizedData_Brain_3_tissues.RData")
cv.sample.array <- c(1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,3,3,4,4,4,4,4,4,4,4,4,4,5,5,5,5,5,5,5,5,5,5,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,3,3,4,4,4,4,4,4,4,4,4,4,5,5,5,5,5,5,5,5,5,5,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,3,3,4,4,4,4,4,4,4,4,4,4,5,5,5,5,5,5,5,5,5,5)

# 2 tissues
load("./NormalizedData_Brain_2_tissues.RData")
cv.sample.array <- c(1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,3,3,4,4,4,4,4,4,4,4,4,4,5,5,5,5,5,5,5,5,5,5,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,3,3,4,4,4,4,4,4,4,4,4,4,5,5,5,5,5,5,5,5,5,5)


## Filter out Top 1000 less Variant Genes
cv_value <- as.data.frame(apply(dataset, 2, ineq, type = "var"))
colnames(cv_value) <- "cv_lab"
dataset <- dataset[,order(cv_value$cv_lab,decreasing = T)[1:1000]]

## create variables
predictionlist_TR <- list()
predictionlist_Val <- list()
predictionlist_TS <- list()
cl_mod <- "rf"
predictionlist_TR[[cl_mod]] <- list(predictions=list(),labels=list())
predictionlist_Val[[cl_mod]] <- list(predictions=list(),labels=list())
predictionlist_TS[[cl_mod]] <- list(predictions=list(),labels=list())
feats_range <- 15:25
idx.cv <- 0
model_cv <-list()
trainingSet_DM_list <- c()
validSet_DM_list <- c()
testSet_DM_list <- c()
stop_time <- c()
start_time <- c()
diff_time <- c()
tr_ctrl <- trainControl(method = "cv", number = 3, returnResamp = "all")
grid=10^seq(10,-2, length =10)


### custom function to extract important predictors
# from LASSO multi-class problem
find_multiclass_lasso_predictors <- function(x){
  idx_uniq <- 0
  for (i in seq_len(length(x))){
    
    idx_uniq <- c(idx_uniq,x[[i]]@i[-1])
  }
  idx_uniq <- idx_uniq[-1]
  idx_uniq <- unique(idx_uniq)
  return(idx_uniq)
}


###########################
# training and validation sets (5-folds CV)
for (crval in 1:5){
  
  start_time[crval] <- Sys.time()
  idx.cv <- which(cv.sample.array != crval)
  
  TR <- dataset[idx.cv,]
  class.TR <- class.df[idx.cv,]
  
  ValTS <- dataset[-idx.cv,]
  class.ValTS <- class.df[-idx.cv,]
  
  #############################################
  ########## LASSO: start
  res <- list()
  set.seed(101010)
  classe <- class.TR[,2,drop=FALSE]
  data_fs <- as.data.frame(cbind(classe,TR))
  
  ######### multi-class classification
  lasso_model <- train(class~.,
                         data=data_fs,
                         method = "glmnet",
                         trControl = tr_ctrl,
                         tuneGrid = expand.grid(alpha = 1, lambda = grid),
                         family="multinomial"
  )
     
  idx_predictors_lasso_all <- coef(lasso_model$finalModel, lasso_model$bestTune$lambda)
  idx_predictors_lasso <- find_multiclass_lasso_predictors(idx_predictors_lasso_all)

  
  #lasso_predictions <- predict(lasso_model,ValTS)
    
  data_reduced<-TR[,idx_predictors_lasso]
  data_reduced <- as.data.frame(data_reduced)
  
  stop_time[crval] <- Sys.time()
  diff_time[crval] <- stop_time[crval] - start_time[crval]
  ########## LASSO: stop
  ########################################################

  ## create formulas and datasets for next classifcation models
  # Training Set
  trainingSet_DM <- cbind(data_reduced,class.TR[,2,drop=FALSE])
  varNames <- colnames(data_reduced)
  colnames(trainingSet_DM) <- c(varNames, "classes")
  varNames1 <- paste(varNames, collapse = "+")
  formula_DM <- as.formula(paste("classes", varNames1, sep = " ~ "))
  
  # Validation Set
  validSet_DM <- cbind(ValTS[,varNames],class.ValTS[,2,drop=FALSE])
  colnames(validSet_DM) <- c(varNames, "classes")
  
  # Independent Test set
  testSet_DM <- cbind(testSet[,varNames],class.TS[,2,drop=FALSE])
  colnames(testSet_DM) <- c(varNames, "classes")
  
  trainingSet_DM_list[[crval]] <- trainingSet_DM
  validSet_DM_list[[crval]] <- validSet_DM
  testSet_DM_list[[crval]] <- testSet_DM
  
  ##########################################################################################
  # Calssification analysis
  set.seed(100)
  model_cv[[crval]] <- randomForest(formula = formula_DM,data = trainingSet_DM, ntree=1000, importance =TRUE)
  
}

