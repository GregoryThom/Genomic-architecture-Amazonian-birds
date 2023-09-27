###5b####

species="lipau"

setwd(paste("~/Dropbox/1_Chapman_fellowship/data/DiploShic/", species, "/observed/", sep = ""))
options(scipen=999)


x <- read.table("total_mask.bed", sep="\t", stringsAsFactors=F)

if (species=="phleg") {
  
chrom <- c("chr1","chr1.1","chr1A","chr1A.1","chr2","chr2.1.1","chr2.1.2","chr3","chr3.1","chr4","chr4A",
"chr4A.1","chr5","chr5.1.1","chr5.1.2","chr5.ch1","chr6","chr7","chr8","chr8.1.1","chr8.1.2","chr9",
"chr10","chr10.1","chr11","chr12","chr12.1","chr14","chr15","chr15.2","chr17","chr18","chr18.1",
"chr20","chr20.1","chr21","chr22LGE","chr23","chr24","chr27","chr28","chrW","chrZ","chrZ.1.1","chrZ.1.2",
"chrZ.1.3")
}

if (species=="xipho") {
  chrom <- c("Pseudo1A", "Pseudo1B", "Pseudo2", "Pseudo3", "Pseudo4A", "Pseudo4", "Pseudo5", "Pseudo6", "Pseudo7", "Pseudo8", 
             "Pseudo9", "Pseudo10", "Pseudo11", "Pseudo12", "Pseudo13", "Pseudo14", "Pseudo15", "Pseudo17", "Pseudo18", "Pseudo19", "Pseudo20", 
             "Pseudo21", "Pseudo22", "Pseudo23", "Pseudo24", "Pseudo25", "Pseudo26", "Pseudo27", "Pseudo28", "PseudoLG2", "PseudoLGE22","PseudoW", "PseudoZ")
}

if (species=="lipau") {
  chrom <- c("Pseudo1", "Pseudo1A", "Pseudo1B", "Pseudo2", "Pseudo3", "Pseudo4A", "Pseudo4", "Pseudo5", "Pseudo6", "Pseudo7", "Pseudo8", 
             "Pseudo9", "Pseudo10", "Pseudo11", "Pseudo12", "Pseudo13", "Pseudo14", "Pseudo15", "Pseudo17", "Pseudo18", "Pseudo19", "Pseudo20", 
             "Pseudo21", "Pseudo22", "Pseudo23", "Pseudo24", "Pseudo25", "Pseudo26", "Pseudo27", "Pseudo28", "PseudoLG2", "PseudoLGE22","PseudoW", "PseudoZ")
}

output_name <- "total_mask_sorted.bed"

for(a in 1:length(chrom)) {
  a_rep <- x[x[,1] == chrom[a], ]
  if(a == 1) {
    write.table(a_rep, file=output_name, sep="\t", quote=F, col.names=F, row.names=F)
  } else {
    write.table(a_rep, file=output_name, sep="\t", quote=F, col.names=F, row.names=F, append=T)
  }
  
}

