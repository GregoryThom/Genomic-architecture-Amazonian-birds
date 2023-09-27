library(seqinr)
library(parallel)
setwd("~/Dropbox/1_Chapman_fellowship/data/2_windows_new/GC/")
a <- read.fasta("rhegma_ragout_complete.fasta")
chrs <- grep("chr", names(a))
#chrs <- seq(1:length(a))
#chr=chrs[1]
  
gc_window <- function(chr) {
name <- names(a)[chr]
b <- a[[chr]]
starts <- seq(1, length(b)-100000, by = 100000)
n <- length(starts)
tab <- NULL
#i=1
for (i in 1:n) {
  chunk <- b[starts[i]:(starts[i]+99999)]
  chunkGC <- GC(chunk)
  s <- as.data.frame(cbind(name, starts[i], c(starts[i]+99999), chunkGC))
  tab <- rbind(tab, s)
}
return(tab)
}

system.time({
  GCcont <- mclapply(chrs, gc_window, mc.cores = 10)
})
GCcont <- do.call(rbind.data.frame, GCcont)
colnames(GCcont) <- c("scaffold", "start", "end", "GC_content")
GCcont <- GCcont[order(GCcont$scaffold, GCcont$start, decreasing = F),]

write.csv(GCcont, "1_lipau_GC.csv")
