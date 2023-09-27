##########function for caret################
run.this.shit <- function(mod,dat,obs,cores){
  mod <- cbind(mod,dat)
  library(caret)
  library(doMC)
  registerDoMC(cores)
  outcomeName <- 'dat'
  predictorsNames <- names(mod)[names(mod) != outcomeName]
  splitIndex <- createDataPartition(mod[,outcomeName], p = .75, list = FALSE, times = 1)
  train <- mod[ splitIndex,]
  test  <<- mod[-splitIndex,]
  
  objControl <- trainControl(method='boot', number=20, returnResamp='final',
                             classProbs = TRUE)
  
  nnetModel <<- train(train[,predictorsNames], train[,outcomeName],
                      method="nnet", maxit=2000,
                      trControl=objControl,
                      metric = "Accuracy",
                      preProc = c("center", "scale"))
  
  predictions <<- predict(object = nnetModel, test[,predictorsNames], type = 'raw')
  accu <<- postResample(pred = predictions, obs = test[,outcomeName])
  pred <<- predict(object = nnetModel, obs, type = 'prob')
  conf <<- confusionMatrix(data = predictions, test$dat)
  return(c(pred, accu))
}


##########function for keras################









