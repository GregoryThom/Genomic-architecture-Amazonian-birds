library(ips)
library(parallel)
library(matrixStats)
library(ggplot2)
numCores <- detectCores()-2
species="lipau_classes"

setwd(paste("/home/gthom/Dropbox/1_Chapman_fellowship/data/11_simulations/10kb/", species, "/", sep = ""))

classes <- c("hard", "linkedHard", "soft", "linkedSoft", "neutral", "WG")

#classes <- c("WG")
plot_vari=NULL
for (i in c(1,4,7)) {
tables=NULL
for (cla in classes) {
t1 <- read.csv(paste(cla, "/Parameter_estimation_Keras_Results.csv", sep = ""), sep = "\t")[i,c(1,4)]

t1$class <- cla
tables <- rbind(tables, t1)
}
colnames(tables) <- c("mean", "SD", "class")
vari <- row.names(tables)[1]
tables$ord <- c("6", "4", "5", "3", "2", "1")
tables <- tables[order(tables$ord),]
rownames(tables) <- NULL
tables$class <- c("1_WholeGenome", "2_Neutral", "3_LinkedSoft", "4_Soft", "5_LinkedHard", "6_LinkedHard")
if (i == 1) {
  vari="Ne (Effective population size)"
}
if (i == 4) {
  vari="TMRCA"
}
if (i == 7) {
  vari="Gene flow (Mn)"
}
data_s <- tables
plot_vari[[i]] <- ggplot(data_s, aes(x=class, y=mean, group=1)) +
  geom_errorbar(width=.2, aes(ymin=mean-SD, ymax=mean+SD), colour="black") +
  geom_point(shape=21, size=3, fill="black") + ylab(vari) + xlab("") + theme_minimal(base_size = 30) +
  theme(axis.text.x = element_text(angle = 30, vjust = 0.5, hjust=0.5))
}


library(gridExtra)

pdf("1_dem_by_classes.pdf", height = 10, width = 30)
grid.arrange(plot_vari[[1]] ,plot_vari[[4]], plot_vari[[7]],nrow=1)
dev.off()

