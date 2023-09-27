# concatenates the prediction output files
species = "lipau"
setwd(paste("~/Dropbox/1_Chapman_fellowship/data/DiploShic/", species, "/observed/", sep=""))

x20 <- as.character(list.files(pattern=("predictions_20"), full.names=F))

if (species =="phleg") {
  a_rep <- read.table("chr1.predictions_20", sep="\t", header=T)
}

if (species =="xipho") {
  a_rep <- read.table("Pseudo1.predictions_20", sep="\t", header=T)
}

if (species =="lipau") {
  a_rep <- read.table("Pseudo1.predictions_20", sep="\t", header=T)
}

x20_output <- c()
#a=2
for(a in 1:length(x20)) {
	a_rep <- read.table(x20[a], sep="\t", header=T)
	x20_output <- rbind(x20_output, a_rep)
	}


write.table(x20_output, file="predictions_20kbp.txt", sep="\t", quote=F, row.names=F)
    
