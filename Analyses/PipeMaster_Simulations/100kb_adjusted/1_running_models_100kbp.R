###Running models pipemaster 
library(PipeMaster)
library(ggplot2)
library(caret)
library(doMC)
library(ggpubr)
library(vroom)
library(keras)
library(matrixStats)


direc="~/Dropbox/1_Chapman_fellowship/data/11_simulations/100kb_adjusted/"
setwd(direc)

#Load the R object with the drawn models (see 1_running_models_10kbp.R)
load("1_models.Rdata.RData")
setwd(direc)

for (species in c("xipho", "phleg", "lipau")) {

#Load new functuions to include intra locus recombination
source("New_functions/sim.msABC.sumstat_new2.R")
source("New_functions/run.carret.R")


#Read population assigment file
  pop.assign <- read.delim("pop_list.txt", header = FALSE)

###Desenhando os modelos

#Size of the windows
loci_size <- 100000


#Adding intralocus recombination rate parameter
if (exists("m1")==F) {
  m1 <- main.menu(m1)
  m2 <- main.menu(m2)
  m3 <- main.menu(m3)
}

#save.image("1_models.Rdata.RData")

if (exists("m1")==F) {
  m1 <- get.data.structure(m1,path.to.fasta = "../template_100K/", pop.assign=pop.assign)
  m2 <- get.data.structure(m2,path.to.fasta = "../template_100K/", pop.assign=pop.assign)
  m3 <- get.data.structure(m3,path.to.fasta = "../template_100K/", pop.assign=pop.assign)
}

#setwd(direc)
  
  if (file.exists("SIMS_m1.txt")==F) {
  
    
#This is where I state the recombination rate as estimated by relernn, to change the recomb rate edit the New_functions/sim.msABC.sumstat_new2
#reco.pars <- t(as.matrix(c("recomb", "rec", "1", "0", "0.000000003021924", "runif")))

#Running simulations for each model    
###m1###  
sim.msABC.sumstat_new2(m1, path="./",
                  #mu.rates = list("rtnorm",1225,5e-10,5e-10,1e-11),
                  nsim.blocks = 32, use.alpha = F,
                  output.name = "m1",
                  append.sims = F, ncores=8, block.size = 400)
    

###m1###  
sim.msABC.sumstat_new2(m2, path="./",
                  #mu.rates = list("rtnorm",1225,5e-10,5e-10,1e-11),
                  nsim.blocks = 32, use.alpha = F,
                  output.name = "m2",
                  append.sims = F, ncores=8, block.size = 400)

###m1###  
sim.msABC.sumstat_new2(m3, path="./",
                  #mu.rates = list("rtnorm",1225,5e-10,5e-10,1e-11),
                  nsim.blocks = 32, use.alpha = F,
                  output.name = "m3",
                  append.sims = F, ncores=8, block.size = 400)

print("It's Done")
}


##############################################################################################################

##reading simulations
if (exists("m1.sim")==F) {
setwd(direc)

m1.sim <- read.delim("SIMS_m1.txt")
m1.sim <- m1.sim[,grep("variance", colnames(m1.sim), invert = T)]
m2.sim <- read.delim("SIMS_m2.txt")
m2.sim <- m2.sim[,grep("variance", colnames(m2.sim), invert = T)]
m3.sim <- read.delim("SIMS_m3.txt")
m3.sim <- m3.sim[,grep("variance", colnames(m3.sim), invert = T)]

#concatenating the data set

#assigning models to simulations
#Subset to run PCA
data_pca <- c(rep("m1",1500),
              rep("m2",1500),
              rep("m3",1500))

#subset for caret/Keras
data_caret <- c(rep("m1",100000),
                rep("m2",100000),
                rep("m3",100000))



models_pca <- rbind(m1.sim[1:1500,c(grep("recomb", colnames(m1.sim))+1):ncol(m1.sim)],
                    m2.sim[1:1500,c(grep("recomb", colnames(m2.sim))+1):ncol(m2.sim)],
                    m3.sim[1:1500,c(grep("recomb", colnames(m3.sim))+1):ncol(m3.sim)])

models_caret <- rbind(m1.sim[1:100000,c(grep("recomb", colnames(m1.sim))+1):ncol(m1.sim)],
                      m2.sim[1:100000,c(grep("recomb", colnames(m2.sim))+1):ncol(m2.sim)],
                      m3.sim[1:100000,c(grep("recomb", colnames(m3.sim))+1):ncol(m3.sim)])


#colnames(models_caret) <- gsub("s_average_", "", colnames(models_caret))

#making sure stats and simulations have the same columns

colnames(models_caret) <- gsub("s_average_", "", colnames(models_caret))


stat <- read.csv(paste("stats_", species, "/chr1_1-100000.out", sep=""), sep = "\t")
stat <- stat[,-grep("Fay",colnames(stat))]
#stat <- stat[,-grep("_Zn",colnames(stat))]
stat <- stat[,-grep("thomson",colnames(stat))]
stat <- stat[,grep("X", colnames(stat), invert = T)]
colnames(stat) <- gsub("s_", "", colnames(stat))


data_caret<-data_caret[complete.cases(models_caret)]
models_caret<-models_caret[complete.cases(models_caret),]

data_pca<-data_pca[complete.cases(models_pca)]
models_pca<-models_pca[complete.cases(models_pca),]
colnames(models_pca) <- gsub("s_average_", "", colnames(models_pca))


# for (i in setdiff(colnames(stat), colnames(models_caret))) {
#   stat <- stat[,-grep(i ,colnames(stat))]
#   }
colnames(stat) <- colnames(models_caret)

if (identical(colnames(stat), colnames(models_caret))==T) {
  print("all good - obs stats and simulations match")
}

#save.image("1_models.Rdata.RData")



######Running Keras####
###Testing models ####
##testing just model 1 and 2
##
rm(model)
if (exists("model")==F) {
dat <- data_caret
mdl <- models_caret
n_models <- 2
if (n_models ==2) {
  dat <- dat[grep("m3", dat, invert = T)]
  mdl <- mdl[grep("m3", dat, invert = T),]
  
}


train = 0.75
test = 0.25

#calculating the mean of each variable to use as scaling point to normalize the data
mean_mdl <- colMeans(mdl)
mean_mdl <- t(as.data.frame(mean_mdl))
#It might be a good idea to normalize the data. Keras have a function for that.
normalize = T
#Normalizing
if (normalize == T) {
  #mdl <- apply(mdl,2,function(x){x/mean(x)})
  mdl <- sweep(mdl, 2, mean_mdl, '/')
  stat2 <- stat/mean_mdl
 # windowsData3 <- sweep(windowsData2, 2, mean_mdl, '/')
  
}
# first convert the data set to array
data <- as.numeric(as.factor(dat))-1
models <- as.matrix(mdl[,1:ncol(mdl)])

models <- cbind(models, data)
dimnames(models) <- NULL
set.seed(12345)
# Determine sample size
ind <- sample(2, nrow(models), replace=TRUE, prob=c(train, test))
models.training <- as.matrix(models[ind==1, 1:c(ncol(models) -1)])
dim(models.training)
models.test <- as.matrix(models[ind==2, 1:c(ncol(models) -1)])
dim(models.test)
models.trainingtarget <- as.matrix(models[ind==1, ncol(models)])
dim(models.trainingtarget)
models.testtarget <- as.matrix(models[ind==2, ncol(models)])
dim(models.testtarget)
##make sure that categories for the models start at 0 and not 1.
models.trainLabels <- as.matrix(to_categorical(models.trainingtarget))
models.testLabels <- as.matrix(to_categorical(models.testtarget))
#model
model <- keras_model_sequential()
# Add layers to the model

if (n_models ==2) {
model %>% 
  layer_dense(units = 36, activation = 'relu', input_shape = c(ncol(models)-1)) %>% 
  layer_dense(units = 36, activation = 'relu') %>% 
  layer_dense(units = 36, activation = 'relu') %>%
  layer_dense(units = 36, activation = 'relu') %>%
  layer_dense(units = 2, activation = 'softmax')
} else {
  model %>% 
    layer_dense(units = 36, activation = 'relu', input_shape = c(ncol(models)-1)) %>% 
    layer_dense(units = 36, activation = 'relu') %>% 
    layer_dense(units = 36, activation = 'relu') %>%
    layer_dense(units = 36, activation = 'relu') %>%
    layer_dense(units = 3, activation = 'softmax')
}
#Define the optimizer
sgd <- optimizer_sgd(lr = 0.1)
# Compile the model
model %>% compile(
  #loss = 'binary_crossentropy',
  loss = 'categorical_crossentropy',
  optimizer = 'adam',
  #optimizer = 'sgd',
  metrics = 'accuracy'
)
# Fit the model 
model %>% fit(
  models.training, 
  models.trainLabels, 
  epochs = 1000, 
  batch_size = nrow(models.trainLabels)/600, 
  validation_split = 0.05
)
# Predict the classes for the test data
classes <- model %>% predict_classes(models.test)
# Confusion matrix
conf <- table(models.testtarget, classes)
write.csv(conf, "1_confusion_matrix_Keras.csv")
# Evaluate on test data and labels
score <- model %>% evaluate(models.test, models.testLabels)
write.csv(score, "1_accuracy_Keras.csv")
# 
# save_model_hdf5(model, "my_model.h5")
# save_model_weights_hdf5(model, "my_model_weights.h5")
# model <- load_model_hdf5("my_model.h5")
# model %>% load_model_weights_hdf5("my_model_weights.h5")
save.image("1_models.Rdata.RData")
# Predict the classes for the test data
}


stat_nnet <- as.matrix(stat2)
dimnames(stat_nnet) <- NULL

# classes_obs <- model %>% predict_classes(stat_nnet)
# classes_obs <- classes_obs +1
# probs <- model %>% predict_proba(stat_nnet)
# colnames(probs) <- c("m1", "m2", "m3")

#write.csv(probs, "1_model_selection_table.csv")
if (normalize == T) {
  #mdl <- apply(mdl,2,function(x){x/mean(x)})
  #mdl <- sweep(mdl, 2, mean_mdl, '/')
  #stat2 <- stat/mean_mdl
  windowsData3 <- sweep(windowsData2, 2, mean_mdl, '/')
  
}

windowsData_nnet <- as.matrix(windowsData3)
dimnames(windowsData_nnet) <- NULL

classes_obs <- model %>% predict_classes(windowsData_nnet)
classes_obs <- classes_obs +1
probs <- model %>% predict_proba(windowsData_nnet)
if (n_models ==2) {
  colnames(probs) <- c("m1", "m2")
} else {
  colnames(probs) <- c("m1", "m2", "m3")
}
model_selection_table <- as.data.frame(cbind(scaffold, start, end, classes_obs, probs))
write.csv(model_selection_table, paste("1_model_selection_table_", species, ".csv", sep = ""))

colMeans(probs, na.rm = T)



#save.image("1_models.Rdata.RData")
######Estimating parameters########
#reading the classification probability of each model

#Selecting model with the highest probability
best_model <- paste("m", grep(max(colMeans(probs, na.rm = T)), colMeans(probs, na.rm = T)), sep = "")

if (best_model == "m1") {
  models_caret <- na.omit(m1.sim[1:100000,])
} 
if (best_model == "m2") {
  models_caret <- na.omit(m2.sim[1:100000,])
}
if (best_model == "m3") {
  models_caret <- na.omit(m3.sim[1:100000,])
}

mdl2 <- models_caret[grep(best_model, data_caret),]
mdl2 <- models_caret[,c(grep("recomb", colnames(models_caret))+1):ncol(models_caret)]
models_params <- models_caret[, 1:c(grep("mean.rate", colnames(models_caret))-1)]
models_params <- models_params[,-c(7,10)]
dat2 <- models_params



#mean is used to normalize the data
mean_mdl2 <- colMeans(mdl2)
mean_dat2 <- colMeans(dat2)

#It might be a good idea to normalize the data. Keras have a function for that.
normalize = T
#Normalizing
if (normalize == T) {
  mdl2 <- sweep(mdl2, 2, mean_mdl2, '/')
  dat2 <- sweep(dat2, 2, mean_dat2, '/')
  stat3 <- stat/mean_mdl2
  windowsData4 <- sweep(windowsData2, 2, mean_mdl2, '/')
}


# first convert the data set to array
data <- as.matrix(dat2[,1:ncol(dat2)])
models <- as.matrix(mdl2[,1:ncol(mdl2)])

stat_nnet2 <- as.matrix(stat3)
dimnames(stat_nnet2) <- NULL

windowsData_nnet2 <- as.matrix(windowsData4)
dimnames(windowsData_nnet2) <- NULL


models <- cbind(models, data)
dimnames(models) <- NULL
set.seed(12345)

#####
##Number of replicates. It is a good idea to run Keras multiple times and 
#calculate the mean of the stats and parameter estimation

#Parameters
#number of replicates
reps <- c(1:10)
#setting the proportion for train and test data
train = 0.75
test = 0.25
epochs = 300

# rep=1
# i=1
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

model2 <- keras_model_sequential() %>%
  layer_dense(units = 32, activation = "relu",
              input_shape = dim(models.training)[2]) %>%
  layer_dense(units = 32, activation = "relu") %>%
  layer_dense(units = 32, activation = "relu") %>%
  layer_dense(units = 1, activation = "relu")

model2 %>% compile(
  loss = "mae",
  optimizer = optimizer_rmsprop(),
  metrics = list("mean_absolute_error"))


model2 %>% fit(
  models.training,
  models.trainingtarget,
  epochs = epochs,
  nrow(models.trainLabels)/600,
  validation_split = 0.05,
  verbose = 0)

mae %<-% (model2 %>% evaluate(models.test, models.testtarget, verbose = 0))

mae <- mae*mean_dat2[i]

test_pred <- model2 %>% predict(models.test)
test_pred <- test_pred*mean_dat2[i]
models.testtarget2 <- models.testtarget*mean_dat2[i]

cor.coef <- cor.test(test_pred, models.testtarget)$estimate

if (dir.exists("Correlations_obs_vs_sim")==F) {
  dir.create("Correlations_obs_vs_sim")
}
pdf(paste("Correlations_obs_vs_sim/",colnames(dat2)[i],"_", rep,"_IBD.pdf",sep=""))
plot(test_pred ~ models.testtarget2, xlab = "true", ylab="estimated", main = colnames(dat2)[i])
lines(c(0:70000),(0:70000), col=2)
dev.off()


predictions <- model2 %>% predict(windowsData_nnet2)
predictions <- predictions*mean_dat2[i]

z <- cbind(predictions, cor.coef, mae[1])
#if (is.na(z[2])) {
#  z[1:3] <- NA
#}
colnames(z)<- c(paste(colnames(dat2)[1:ncol(dat2)][i], "_", rep, sep = ""), paste("r2_", rep, sep = ""), paste("MAE_", rep, sep = ""))
Res <- cbind(Res, z)

}
#rownames(Res) <- colnames(dat2)[1:ncol(dat2)]
print(paste("replicate ",rep, " is done!", sep = ""))
return(Res)
}

#####Running
system.time({
  parameters_keras <- lapply(reps, parameters_keras)
})
f <- do.call(cbind.data.frame, parameters_keras)
#e=colnames(dat2)[1]
parameters_keras_mean=NULL
for (e in colnames(dat2)) {
c <- grep(e , colnames(f))
r <- rowMeans(f[,c], na.rm = T)
c2 <- grep(e , colnames(f))+1
r2 <- rowMeans(f[,c2], na.rm = T)
c3 <- grep(e , colnames(f))+2
r3 <- rowMeans(f[,c3], na.rm = T)
r4 <- cbind(r,r2,r3)
colnames(r4) <- c(e, paste(e, "_R2", sep=""), paste(e, "_MAE", sep=""))
parameters_keras_mean <- cbind(parameters_keras_mean, r4)
}

parameters_keras_mean <- cbind(model_selection_table, parameters_keras_mean)

write.table(parameters_keras_mean, paste("Parameter_estimation_Keras_Results_", species, ".csv", sep=""), row.names = T, sep = "\t")


save.image("1_models.Rdata.RData")

#plotting PCA to test the fit of the simulations to the observed data

PCA <- prcomp(rbind(models_pca, stat), center = T, scale. = T, retx=T)
scores <- data.frame(PCA$x[,1:ncol(PCA$x)])
p1x2<-ggplot(scores, aes(x=PC1, y=PC2))+ theme(legend.position = "none") +  theme(legend.title=element_blank()) +
  geom_point(aes(colour=c(data_pca,"observed"), size=c(data_pca,"observed"), shape=c(data_pca,"observed")))+
  scale_size_manual(values=c(0.7,0.7,0.7,4))+
  scale_colour_manual(values = c("green", "blue", "red",  "black")) + scale_shape_manual(values = c(1,2,3,16))

p1x3<-ggplot(scores, aes(x=PC1, y=PC3))+ theme(legend.position = "none") + theme(legend.title=element_blank()) +
  geom_point(aes(colour=c(data_pca,"observed"), size=c(data_pca,"observed"), shape=c(data_pca,"observed")))+
  scale_size_manual(values=c(0.7,0.7,0.7,4))+
  scale_colour_manual(values = c("green", "blue", "red",  "black")) + scale_shape_manual(values = c(1,2,3,16))

p1x4<-ggplot(scores, aes(x=PC1, y=PC4))+ theme(legend.position = "none") + theme(legend.title=element_blank()) +
  geom_point(aes(colour=c(data_pca,"observed"), size=c(data_pca,"observed"), shape=c(data_pca,"observed")))+
  scale_size_manual(values=c(0.7,0.7,0.7,4))+
  scale_colour_manual(values = c("green", "blue", "red",  "black")) + scale_shape_manual(values = c(1,2,3,16))

p2x3<-ggplot(scores, aes(x=PC2, y=PC3))+ theme(legend.position = "none") + theme(legend.title=element_blank()) +
  geom_point(aes(colour=c(data_pca,"observed"), size=c(data_pca,"observed"), shape=c(data_pca,"observed")))+
  scale_size_manual(values=c(0.7,0.7,0.7,4))+
  scale_colour_manual(values = c("green", "blue", "red",  "black")) + scale_shape_manual(values = c(1,2,3,16))

p2x4<-ggplot(scores, aes(x=PC2, y=PC4))+ theme(legend.position = "none") + theme(legend.title=element_blank()) +
  geom_point(aes(colour=c(data_pca,"observed"), size=c(data_pca,"observed"), shape=c(data_pca,"observed")))+
  scale_size_manual(values=c(0.7,0.7,0.7,4))+
  scale_colour_manual(values = c("green", "blue", "red", "black")) + scale_shape_manual(values = c(1,2,3,16))

p3x4<-ggplot(scores, aes(x=PC3, y=PC4))+ theme(legend.position = "none") + theme(legend.title=element_blank()) +
  geom_point(aes(colour=c(data_pca,"observed"), size=c(data_pca,"observed"), shape=c(data_pca,"observed")))+
  scale_size_manual(values=c(0.7,0.7,0.7,4))+
  scale_colour_manual(values = c("green", "blue", "red", "black")) + scale_shape_manual(values = c(1,2,3,16))

PCA_Plot <- ggarrange(p1x2,p1x3, p1x4, p2x3, p2x4, p3x4, legend = "left", common.legend = T) + ggtitle(paste(species))
PCA_Plot
pdf(paste(species, "_PCA.pdf", sep = ""), width=12, height=6.69)
PCA_Plot
dev.off()


#Plotting variable contributions for summary stats
library("FactoMineR")
library("factoextra")

res.pca <- PCA(rbind(models_pca), graph = FALSE, scale.unit = TRUE, ncp=10)
eig.val <- get_eigenvalue(res.pca)
fviz_eig(res.pca, addlabels = TRUE, ylim = c(0, 50))
pdf("1_PC_contribution.pdf")
fviz_eig(res.pca, addlabels = TRUE, ylim = c(0, 50))
dev.off()

var <- get_pca_var(res.pca)

library("corrplot")
corrplot(var$cos2, is.corr=FALSE, tl.cex=0.5)
pdf("1_PCs_vs_variables_cor.pdf")
corrplot(var$cos2, is.corr=FALSE, tl.cex=0.5)
dev.off()





