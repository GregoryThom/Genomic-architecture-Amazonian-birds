library(ips)
library(parallel)
library(MASS)
library(plyr)
library(tidyverse)
library(fuzzyjoin)
numCores <- detectCores()-2
species="phleg"

for (species in c("phleg", "xipho", "lipau")) {
 
q <- c("100")


###########formating big table######
setwd(paste("~/Dropbox/1_Chapman_fellowship/data/2_windows_new/", species, sep = ""))
source("../../scripts/analyses/plot_twisst.R")

snps=q

big_table <- read.csv(paste(snps, "/1_all_stats_", species, "_", snps, ".csv", sep = ""), header = T)

#Reading and combining twist
#wheight <- read.csv(paste(snps, "/SNP_phased.weights.csv", sep = ""), sep="\t", header = T)
twist_windows <- read.csv(paste(snps,"/SNP-only_phased.phyml_bionj.w", snps, ".data.tsv", sep = ""), sep="\t", header = T)
trees <- read.csv(paste(snps, "/SNP-only_phased.phyml_bionj.w", snps, ".trees", sep = ""), header = F, sep = "\t")
twist <- cbind(twist_windows, trees)
twist$scaffold <- gsub("Pseudo", "chr", twist$scaffold )
twist$scaffold <- gsub("chr_", "chr", twist$scaffold)
twist$scaffold <- gsub("chrNC", "chrW", twist$scaffold )
twist <- na.omit(twist[grep("chr", x = twist$scaffold),])
twist <- twist[order(twist$scaffold, twist$start),]
levels(twist$scaffold) <- as.character(unique(twist$scaffold))

twist_match <- function(chr){
  chr_w = subset(dfg, scaffold == chr)
  chr_r <- subset(twist, scaffold == chr)
  twist_prop <- NULL
  if (nrow(chr_w)>0 && nrow(chr_r)>0) {
    for (z in 1:nrow(chr_w)) {
      for (i in 1:nrow(chr_r)) {
        if (chr_r$start[i] > chr_w$start[z] && chr_r$start[i] < chr_w$end[z]) {
          twist_prop <- rbind(twist_prop, rownames(chr_r)[i]) 
        } 
      }
    }
    return(twist_prop)
  }
}

#loop to extract the genetrees and coordinates of windows with specific thresholds for stats
#These are the stats
#d="recombRate"
for (d in c("recombRate", "dxy", "Fst", "pi", "fdM")) {
#these are the thresholds Lowers than 10 or upper than 90%  
#s="10"
#s=10
for (s in c("10", "90")) {
a <- as.numeric(grep(d, colnames(big_table)))[1]
big_table[a] <- rowMeans(big_table[grep(d, colnames(big_table))])
if (as.numeric(s) == 10) {
  dfg <- big_table[big_table[,a] < quantile(big_table[,a], 0.1, na.rm = T),]
} else { 
  dfg <- big_table[big_table[,a] > quantile(big_table[,a], 0.9, na.rm = T),]
}
dfg <- dfg[is.na(dfg$start)==F,]
twist <- twist[is.na(twist$start)==F,]
chrs <- levels(big_table$scaffold)
if (species == "lipau") {
  chrs <- chrs[-grep("chr_W", chrs)]
} else {
  chrs <- chrs[-grep("chrW", chrs)]
}

system.time({
  dfg_trees <- mclapply(chrs, twist_match, mc.cores = numCores)
})
dfg_trees <- do.call(rbind.data.frame, dfg_trees)
dfg_trees <- twist[as.vector(dfg_trees$V1),]
dfg_trees <- dfg_trees[order(dfg_trees$scaffold, dfg_trees$start),]
#Create a folder with the name
if (file.exists(paste(snps, "/", d, "_", s, "_trees", sep = "")) == F){
system(paste("mkdir ", snps, "/", d, "_", s, "_trees", sep = ""))
}
#exporting gene trees 
write.csv(dfg_trees$V1, paste(snps, "/", d, "_", s, "_trees/", d, "_", s, "_gentrees.txt", sep = ""), row.names = F, quote = F )
system(paste("sed -i 1d ", snps, "/", d, "_", s, "_trees/", d, "_", s, "_gentrees.txt", sep = ""))
#exporting bed file 
bed <- cbind(as.character(dfg_trees$scaffold), dfg_trees$start, dfg_trees$end)
colnames(bed) <- c("scaffold", "start", "end")
write.table(bed, paste(snps, "/", d, "_", s, "_trees/", d, "_", s, ".bed", sep = ""), col.names = T, row.names = F, quote = F, sep = "\t")
  }
}


#mtDNA and UCEs
#d="UCE"
for (d in c("UCE", "mtDNA")) {
    a <- as.numeric(grep(d, colnames(big_table)))[1]
    dfg <- big_table[grep("1", big_table[,a]),]
    dfg <- dfg[is.na(dfg$start)==F,]
      system.time({
      dfg_trees <- mclapply(chrs, twist_match, mc.cores = numCores)
    })
      dfg_trees <- do.call(rbind.data.frame, dfg_trees)
      dfg_trees <- twist[as.vector(dfg_trees$V1),]
      dfg_trees <- dfg_trees[order(dfg_trees$scaffold, dfg_trees$start),]
    #Create a folder with the name
    if (file.exists(paste(snps, "/", d, "_trees", sep = "")) == F){
      system(paste("mkdir ", snps, "/", d, "_trees", sep = ""))
    }
    #exporting gene trees 
    write.csv(dfg_trees$V1, paste(snps, "/", d, "_trees/", d, "_gentrees.txt", sep = ""), row.names = F, quote = F )
    system(paste("sed -i 1d ", snps, "/", d, "_trees/", d, "_gentrees.txt", sep = ""))
    #exporting bed file 
    bed <- cbind(as.character(dfg_trees$scaffold), dfg_trees$start, dfg_trees$end)
    colnames(bed) <- c("scaffold", "start", "end")
    write.table(bed, paste(snps, "/", d, "_trees/", d, ".bed", sep = ""), col.names = T, row.names = F, quote = F, sep = "\t")
    }


#chromossomes
#d="UCE"
chrs <- levels(big_table$scaffold)
if (species == "lipau") {
  chrs <- chrs[-grep("chr_W", chrs)]
} else {
  chrs <- chrs[-grep("chrW", chrs)]
}
d=chrs[1]
for (d in chrs) {
  dfg_trees <- subset(twist, scaffold == d)
  #Create a folder with the name
  if (file.exists(paste(snps, "/", d, "_trees", sep = "")) == F){
    system(paste("mkdir ", snps, "/", d, "_trees", sep = ""))
  }
  #exporting gene trees 
  write.csv(dfg_trees$V1, paste(snps, "/", d, "_trees/", d, "_gentrees.txt", sep = ""), row.names = F, quote = F )
  system(paste("sed -i 1d ", snps, "/", d, "_trees/", d, "_gentrees.txt", sep = ""))
  #exporting bed file 
  bed <- cbind(as.character(dfg_trees$scaffold), dfg_trees$start, dfg_trees$end)
  colnames(bed) <- c("scaffold", "start", "end")
  write.table(bed, paste(snps, "/", d, "_trees/", d, ".bed", sep = ""), col.names = T, row.names = F, quote = F, sep = "\t")
}

if (species == "phleg") {
  saa <- list.dirs(paste(snps, sep = ""), full.names = F)[grep("chr", list.dirs(paste(snps, sep = ""), full.names = F))]
  sa <- unique(gsub("\\.[1-2].*|_trees", "", saa))
  sa <- sa[-grep("chr5.ch1", sa)]
  #i=sa[1]
  for (i in sa[-1]) {
    system(paste("cat ", snps, "/", i, ".[1-2]*/*.txt >> ", snps, "/", i, "_trees/", i, "_gentrees.txt", sep=""))
    system(paste("rm -r ", snps, "/", i, ".*_trees", sep=""))
    }
}

#all trees
d="all"
if (file.exists(paste(snps, "/", d, "_trees", sep = "")) == F){
  system(paste("mkdir ", snps, "/", d, "_trees", sep = ""))
}
#exporting gene trees 
system(paste("cp ", snps, "/SNP-only_phased.phyml_bionj.w", snps, ".trees ", snps, "/", d, "_trees/", d, "_gentrees.txt" , sep = ""))
#exporting bed file 
bed <- cbind(as.character(twist$scaffold), twist$start, twist$end)
colnames(bed) <- c("scaffold", "start", "end")
write.table(bed, paste(snps, "/", d, "_trees/", d, ".bed", sep = ""), col.names = T, row.names = F, quote = F, sep = "\t")

setwd(paste("~/Dropbox/1_Chapman_fellowship/data/2_windows_new/", species, sep = ""))

if (file.exists("bed_for_tree") == F){
  system("mkdir bed_for_tree")
}
#folders="dxy_10_trees"
for (snps in q) {
  for (folders in grep("all", list.dirs("100/", full.names = F)[-1], invert = T, value = T)) {
    system(paste("cp ", snps, "/", folders, "/*.bed bed_for_tree/", snps, "_", folders, ".bed", sep=""))
  }
}

}
