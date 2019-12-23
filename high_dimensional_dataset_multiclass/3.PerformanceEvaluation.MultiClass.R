library(randomForest)

for (crval in 1:5){
    trainingSet_DM <- trainingSet_DM_list[[crval]]
    validSet_DM <- validSet_DM_list[[crval]] 
    
    n_predictors <- length(colnames(trainingSet_DM_list[[crval]]))-1
    
    # Predict on Training Set
    caret::confusionMatrix(table(predict(model_cv[[crval]], trainingSet_DM),trainingSet_DM$classes),reference = trainingSet_DM$classes)
    caret::confusionMatrix(table(predict(model_cv[[crval]], trainingSet_DM),trainingSet_DM$classes),reference = trainingSet_DM$classes)[["byClass"]]
    
    # Predict on Validation Set
    caret::confusionMatrix(table(predict(model_cv[[crval]], validSet_DM),validSet_DM$classes),reference = validSet_DM$classes)
    caret::confusionMatrix(table(predict(model_cv[[crval]], validSet_DM),validSet_DM$classes),reference = validSet_DM$classes)[["byClass"]]

    }

# Test best model on Independent Test Set
index_best_model <- 4

testSet_DM <- testSet_DM_list[[index_best_model]] 
caret::confusionMatrix(table(predict(model_cv[[index_best_model]], testSet_DM),testSet_DM$classes),reference = testSet_DM$classes)
caret::confusionMatrix(table(predict(model_cv[[index_best_model]], testSet_DM),testSet_DM$classes),reference = testSet_DM$classes)[["byClass"]]

