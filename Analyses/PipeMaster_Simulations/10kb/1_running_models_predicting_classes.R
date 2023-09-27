#Calculating demographic parameters for regions of the genome assigned to different models in diploS/HIC

library(PipeMaster)
library(ggplot2)
library(caret)
library(doMC)
library(ggpubr)
library(vroom)
library(keras)
library(matrixStats)
species="lipau_classes"
direc=paste("~/Dropbox/1_Chapman_fellowship/data/11_simulations/", "10kb", sep = "")
setwd(direc)

classes <- c("hard", "linkedHard", "soft", "linkedSoft", "neutral")

###read sumstats (reading msABC output .out and calculating mean and variance)
###There are too many windows and it is taking forever to simulate. I'm selecting 1 10km window every 100mb
#ls loci_summary_stats/ | awk 'NR % 5 == 0' | grep -v chrZ | grep -v chrW | parallel cp loci_summary_stats/{} loci_summary_stats_10kb_10/

load("1_models.Rdata.RData")

classes <- c("hard", "linkedHard", "soft", "linkedSoft", "neutral", "WG")

#cla <- c("WG")
for (cla in classes) {
  
if (file.exists(paste(species, "/", cla, "/stats_obs.csv", sep = ""))==F) {
  obs_data <- do.call(rbind, lapply(list.files(path = paste(species, "/", cla, "/loci_summary_stats/", sep = ""), full.names = T), read.csv, sep = "\t"))
  obs_data <- obs_data[,colnames(obs_data) != "X"]
  stats <- colnames(obs_data)
  stat=stats[1]
  mean_var_obs <- lapply(stats, function(stat) {
    stat_mean <- mean(obs_data[[stat]], na.rm = T)
    stat_vari <- var(obs_data[[stat]], na.rm = T)
    both <- cbind(stat_mean, stat_vari)
    colnames(both) <- c(gsub("s_", "s_average_", stat), gsub("s_", "s_variance_", stat))
    return(both)
  }
  )
  mean_var_obs <- do.call(cbind.data.frame, mean_var_obs)
  write.csv(mean_var_obs, paste(species, "/", cla, "/stats_obs.csv", sep = ""))
} else {
  mean_var_obs <- read.csv(paste(species, "/", cla, "/stats_obs.csv", sep = ""))
}

stat <- mean_var_obs

stat <- stat[,grep("X|thomson|ZnS|FayWuH|dvk|dvh|tajimasD" ,colnames(stat), invert = T)]

colnames(stat) <- colnames(models_caret)

if (species == "phleg_classes") {
  best_model <- "m1"
}
if (species == "xipho_classes") {
  best_model <- "m2"
}

if (species == "lipau_classes") {
  best_model <- "m1"
}

mdl2 <- models_caret[grep(best_model, data_caret),]
models_params <- get(paste(best_model, ".sim", sep = ""))
models_params <- models_params[1:40000, 1:c(grep("mean.rate", colnames(models_params))-1)]
models_params <- models_params[,-c(7,10)]
dat2 <- models_params
dat2 <- dat2[,c(1,2,3,4,5,6,8)]


#mean is used to normalize the data
mean_mdl2 <- colMeans(mdl2)
mean_dat2 <- colMeans(dat2)

#It might be a good idea to normalize the data. Keras have a function for that.
normalize = T
#Normalizing
if (normalize == T) {
  mdl2 <- apply(mdl2,2,function(x){x/mean(x)})
  dat2 <- apply(dat2,2,function(x){x/mean(x)})
  stat3 <- stat/mean_mdl2
}


# first convert the data set to array
data <- as.matrix(dat2[,1:ncol(dat2)])
models <- as.matrix(mdl2[,1:ncol(mdl2)])


models <- cbind(models, data)
dimnames(models) <- NULL
set.seed(12345)

#####
##Number of replicates. It is a good idea to run Keras multiple times and 
#calculate the mean of the stats and parameter estimation

#Parameters
#number of replicates
reps <- c(1:2)
#setting the proportion for train and test data
train = 0.75
test = 0.25
epochs = 300


##Running keras
parameters_keras <- function(rep) {
Res <- NULL
##here goes the loop for each variable 1:9
#i=8
for (i in 1:ncol(dat2)) {
# Determine sample size
ind <- sample(2, nrow(models), replace=TRUE, prob=c(train, test))
models.training <- as.matrix(models[ind==1, 1:c(ncol(models) - ncol(data))])
dim(models.training)
models.test <- as.matrix(models[ind==2, 1:c(ncol(models) - ncol(data))])
dim(models.test)
models.trainingtarget <- as.matrix(models[ind==1, c(ncol(mdl2)+i)])
dim(models.trainingtarget)
models.testtarget <- as.matrix(models[ind==2, c(ncol(mdl2)+i)])
dim(models.testtarget)

model <- keras_model_sequential() %>%
  layer_dense(units = 32, activation = "relu",
              input_shape = dim(models.training)[2]) %>%
  layer_dense(units = 32, activation = "relu") %>%
  layer_dense(units = 32, activation = "relu") %>%
  layer_dense(units = 1, activation = "relu")

model %>% compile(
  loss = "mae",
  optimizer = optimizer_rmsprop(),
  metrics = list("mean_absolute_error"))


model %>% fit(
  models.training,
  models.trainingtarget,
  epochs = epochs,
  nrow(models.trainingtarget)/600,
  validation_split = 0.05,
  verbose = 0)

mae %<-% (model %>% evaluate(models.test, models.testtarget, verbose = 0))

mae <- mae*mean_dat2[i]

test_pred <- model %>% predict(models.test)
test_pred <- test_pred*mean_dat2[i]
models.testtarget2 <- models.testtarget*mean_dat2[i]

cor.coef <- cor.test(test_pred, models.testtarget)$estimate

if (dir.exists(paste(species, "/", cla, "/Correlations_obs_vs_sim",sep=""))==F) {
  dir.create(paste(species, "/", cla, "/Correlations_obs_vs_sim",sep=""))
}
pdf(paste(species, "/", cla, "/Correlations_obs_vs_sim/",colnames(dat2)[i],"_", rep,"_IBD.pdf",sep=""))
plot(test_pred ~ models.testtarget2, xlab = "true", ylab="estimated", main = colnames(dat2)[i])
lines(c(0:1000000),(0:1000000), col=2)
dev.off()

stat_nnet <- as.matrix(stat3)
dimnames(stat_nnet) <- NULL

predictions <- model %>% predict(stat_nnet)
predictions <- predictions*mean_dat2[i]

z <- c(predictions, cor.coef, mae[1])
if (is.na(z[2])) {
  z[1:3] <- NA
}
names(z)<- c(paste("estimate_", rep, sep = ""), paste("r2_", rep, sep = ""), paste("MAE_", rep, sep = ""))
Res <- rbind(Res, z)

}
rownames(Res) <- colnames(dat2)[1:ncol(dat2)]
print(paste("replicate ",rep, " is done!", sep = ""))
return(Res)
}

#####Running
system.time({
  parameters_keras <- lapply(reps, parameters_keras)
})
parameters_keras <- do.call(cbind.data.frame, parameters_keras)

esti <- parameters_keras[,grep("estimate", colnames(parameters_keras))]
esti_mean <- as.data.frame(rowMeans(esti, na.rm = T))
esti_SD <- as.data.frame(rowSds(as.matrix(esti), na.rm = T))

r2 <- parameters_keras[,grep("r2", colnames(parameters_keras))]
r2_mean <- as.data.frame(rowMeans(r2, na.rm = T))

MAE <- parameters_keras[,grep("MAE", colnames(parameters_keras))]
MAE_mean <- as.data.frame(rowMeans(MAE, na.rm = T))

mean_keras_replicates <- cbind(esti_mean, esti_SD, r2_mean, MAE_mean)

colnames(mean_keras_replicates) <- c(paste("mean of ", max(reps), " replicates", sep = ""), "SD", "R^2", "MAE")

write.table(round(mean_keras_replicates, 3), paste(species, "/", cla, "/Parameter_estimation_Keras_Results.csv" ,sep=""), row.names = T, sep = "\t")

}

