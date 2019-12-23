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
library(microbenchmark)
library(ROCR)

## Load data
load("./NormalizedDatasets.RData")

## cross validation sampling order (5-fold CV)
cv.sample.array <- c(1,1,1,1,1,2,2,2,2,2,3,3,3,3,3,4,4,4,4,4,5,5,5,5,5,1,1,1,1,1,2,2,2,2,2,3,3,3,3,3,4,4,4,4,4,5,5,5,5,5)

## create variables
predictionlist_TR <- list()
predictionlist_Val <- list()
predictionlist_TS <- list()
cl_mod <- "rf"
predictionlist_TR[[cl_mod]] <- list(predictions=list(),labels=list())
predictionlist_Val[[cl_mod]] <- list(predictions=list(),labels=list())
predictionlist_TS[[cl_mod]] <- list(predictions=list(),labels=list())
feats_range <- 5:20
idx.cv <- 0
model_cv <-list()
trainingSet_DM_list <- c()
validSet_DM_list <- c()
testSet_DM_list <- c()
stop_time <- c()
start_time <- c()
diff_time <- c()
control <- rfeControl(functions=lrFuncs, method="cv",number=3)

res <- matrix(nrow = 5,ncol = 9)
colnames(res) <- c("Accuracy", "Sensitivity","Specificity","PPV","NPV","Precision","Recall","AUC", "n_predictors")
res_TR <- res
res_Val <- res

### Training and validation sets
for (crval in 1:5){
  
  start_time[crval] <- Sys.time()

  idx.cv <- which(cv.sample.array != crval)
  TR <- dataset[idx.cv,]
  class.TR <- class.df[idx.cv,,drop=F]
  
  ValTS <- dataset[-idx.cv,]
  class.ValTS <- class.df[-idx.cv,,drop=F]
  
  #############################################
  ########## RFE: start
  res <- list()
  set.seed(101010)
 
    classe <- class.TR[,1,drop=FALSE]
    data_fs <- as.data.frame(cbind(classe,TR))
    results <- rfe(data_fs[,2:dim(TR)[2]], data_fs[,1], feats_range, rfeControl=control)
    data_reduced<-TR[,results$optVariables]
    data_reduced <- as.data.frame(data_reduced)  
	
    stop_time[crval] <- Sys.time()
    diff_time[crval] <- as.numeric(stop_time[crval] - start_time[crval])
  ########## RFE: stop
  ########################################################

  ## create formulas and datasets for next classifcation models
  # Training Set
  trainingSet_DM <- cbind(data_reduced,class.TR[,1,drop=FALSE])
  varNames <- colnames(data_reduced)
  colnames(trainingSet_DM) <- c(varNames, "classes")
  varNames1 <- paste(varNames, collapse = "+")
  formula_DM <- as.formula(paste("classes", varNames1, sep = " ~ "))
  
  # Validation Set
  validSet_DM <- cbind(ValTS[,varNames],class.ValTS[,1,drop=FALSE])
  colnames(validSet_DM) <- c(varNames, "classes")
  
  # Independent Test set
  testSet_DM <- cbind(testSet[,varNames],class.TS[,1,drop=FALSE])
  colnames(testSet_DM) <- c(varNames, "classes")
  
  trainingSet_DM_list[[crval]] <- trainingSet_DM
  validSet_DM_list[[crval]] <- validSet_DM
  testSet_DM_list[[crval]] <- testSet_DM
  
  ## Classification analysis by random forest
  set.seed(100)
  model_cv[[crval]] <- randomForest(formula = formula_DM,data = trainingSet_DM, ntree=1000, importance =TRUE)

  ## Performance Evaluation 
  trainingSet_DM <- trainingSet_DM_list[[crval]]
  validSet_DM <- validSet_DM_list[[crval]] 
  n_predictors <- length(colnames(trainingSet_DM_list[[crval]]))-1
  
  acc_model_TR <- caret::confusionMatrix(table(predict(model_cv[[crval]], trainingSet_DM),trainingSet_DM$classes),reference = trainingSet_DM$classes)$overall['Accuracy']
  acc_model_Val <- caret::confusionMatrix(table(predict(model_cv[[crval]], validSet_DM),validSet_DM$classes),reference = validSet_DM$classes)$overall['Accuracy']
  
  sen_model_TR <- caret::confusionMatrix(table(predict(model_cv[[crval]], trainingSet_DM),trainingSet_DM$classes),reference = trainingSet_DM$classes)$byClass['Sensitivity']
  sen_model_Val <- caret::confusionMatrix(table(predict(model_cv[[crval]], validSet_DM),validSet_DM$classes),reference = validSet_DM$classes)$byClass['Sensitivity']
  
  spe_model_TR <- caret::confusionMatrix(table(predict(model_cv[[crval]], trainingSet_DM),trainingSet_DM$classes),reference = trainingSet_DM$classes)$byClass['Specificity']
  spe_model_Val <- caret::confusionMatrix(table(predict(model_cv[[crval]], validSet_DM),validSet_DM$classes),reference = validSet_DM$classes)$byClass['Specificity']
  
  ppv_model_TR <- caret::confusionMatrix(table(predict(model_cv[[crval]], trainingSet_DM),trainingSet_DM$classes),reference = trainingSet_DM$classes)$byClass['Pos Pred Value']
  ppv_model_Val <- caret::confusionMatrix(table(predict(model_cv[[crval]], validSet_DM),validSet_DM$classes),reference = validSet_DM$classes)$byClass['Pos Pred Value']
  
  npv_model_TR <- caret::confusionMatrix(table(predict(model_cv[[crval]], trainingSet_DM),trainingSet_DM$classes),reference = trainingSet_DM$classes)$byClass['Neg Pred Value']
  npv_model_Val <- caret::confusionMatrix(table(predict(model_cv[[crval]], validSet_DM),validSet_DM$classes),reference = validSet_DM$classes)$byClass['Neg Pred Value']
  
  pre_model_TR <- caret::confusionMatrix(table(predict(model_cv[[crval]], trainingSet_DM),trainingSet_DM$classes),reference = trainingSet_DM$classes)$byClass['Precision']
  pre_model_Val <- caret::confusionMatrix(table(predict(model_cv[[crval]], validSet_DM),validSet_DM$classes),reference = validSet_DM$classes)$byClass['Precision']
  
  rec_model_TR <- caret::confusionMatrix(table(predict(model_cv[[crval]], trainingSet_DM),trainingSet_DM$classes),reference = trainingSet_DM$classes)$byClass['Recall']
  rec_model_Val <- caret::confusionMatrix(table(predict(model_cv[[crval]], validSet_DM),validSet_DM$classes),reference = validSet_DM$classes)$byClass['Recall']
  
  # AUC
  predictions_rf=predict(model_cv[[crval]], trainingSet_DM, type="prob")[,1]
  pred=prediction(predictions_rf,1-(as.numeric(trainingSet_DM$classes) - 1))
  perf_AUC=performance(pred,"auc")
  AUC_TR=perf_AUC@y.values[[1]]
  
  predictionlist_TR[[cl_mod]]$predictions[[crval]] <- predictions_rf
  predictionlist_TR[[cl_mod]]$labels[[crval]] <- 1-(as.numeric(trainingSet_DM$classes) - 1)
  
  predictions_rf=predict(model_cv[[crval]], validSet_DM, type="prob")[,1] 
  pred=prediction(predictions_rf,1-(as.numeric(validSet_DM$classes) - 1)) 
  perf_AUC=performance(pred,"auc") 
  AUC_Val=perf_AUC@y.values[[1]]
  
  predictionlist_Val[[cl_mod]]$predictions[[crval]] <- predictions_rf
  predictionlist_Val[[cl_mod]]$labels[[crval]] <- 1-(as.numeric(validSet_DM$classes) - 1)
  
  # Peformance summary tables
  res_TR[crval,1] <- acc_model_TR
  res_TR[crval,2] <- sen_model_TR
  res_TR[crval,3] <- spe_model_TR
  res_TR[crval,4] <- ppv_model_TR
  res_TR[crval,5] <- npv_model_TR
  res_TR[crval,6] <- pre_model_TR
  res_TR[crval,7] <- rec_model_TR
  res_TR[crval,8] <- AUC_TR
  res_TR[crval,9] <- n_predictors
  
  res_Val[crval,1] <- acc_model_Val
  res_Val[crval,2] <- sen_model_Val
  res_Val[crval,3] <- spe_model_Val
  res_Val[crval,4] <- ppv_model_Val
  res_Val[crval,5] <- npv_model_Val
  res_Val[crval,6] <- pre_model_Val
  res_Val[crval,7] <- rec_model_Val
  res_Val[crval,8] <- AUC_Val
  res_Val[crval,9] <- n_predictors
   
}

####################### test the best model on the Idependent test set
best_model_index <- 1

testSet_DM <- testSet_DM_list[[best_model_index]] 
acc_model_TS <- caret::confusionMatrix(table(predict(model_cv[[best_model_index]], testSet_DM),testSet_DM$classes),reference = testSet_DM$classes)$overall['Accuracy']
sen_model_TS <- caret::confusionMatrix(table(predict(model_cv[[best_model_index]], testSet_DM),testSet_DM$classes),reference = testSet_DM$classes)$byClass['Sensitivity']
spe_model_TS <- caret::confusionMatrix(table(predict(model_cv[[best_model_index]], testSet_DM),testSet_DM$classes),reference = testSet_DM$classes)$byClass['Specificity']
ppv_model_TS <- caret::confusionMatrix(table(predict(model_cv[[best_model_index]], testSet_DM),testSet_DM$classes),reference = testSet_DM$classes)$byClass['Pos Pred Value']
npv_model_TS <- caret::confusionMatrix(table(predict(model_cv[[best_model_index]], testSet_DM),testSet_DM$classes),reference = testSet_DM$classes)$byClass['Neg Pred Value']
pre_model_TS <- caret::confusionMatrix(table(predict(model_cv[[best_model_index]], testSet_DM),testSet_DM$classes),reference = testSet_DM$classes)$byClass['Precision']
rec_model_TS <- caret::confusionMatrix(table(predict(model_cv[[best_model_index]], testSet_DM),testSet_DM$classes),reference = testSet_DM$classes)$byClass['Recall']

predictions_rf=predict(model_cv[[best_model_index]], testSet_DM, type="prob")[,1] 
pred=prediction(predictions_rf,1-(as.numeric(testSet_DM$classes) - 1)) 
perf_AUC=performance(pred,"auc") #Calculate the AUC value
AUC_TS=perf_AUC@y.values[[1]]

