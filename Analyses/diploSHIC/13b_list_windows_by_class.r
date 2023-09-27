
sel <- read.csv("predictions_20kbp.txt", sep="\t")

big_table2 <- sel

chrs=as.character(unique(sel$chrom))

#d=6
for (d in grep("prob", colnames(big_table2))) {
  e <- colnames(big_table2[d])
  red <- big_table2[big_table2[,d] > 0.8,]
  if (nrow(red) < 1000) {
    red <- big_table2[big_table2[,d] > 0.5,]
  }
  wndw <- red[,1:3]
  colnames(wndw) <- c("scaffold", "start", "end")
  if (nrow(wndw) > 1000) {
  wndw <- wndw[sample(c(1:nrow(wndw)), 1000, replace = F),]
  }
  wndw$end <-  wndw$end - 10000
  print(paste(e, " ", nrow(wndw)))
  wndw <- wndw[order(wndw$scaffold, wndw$start),]
  write.table(wndw, paste("1_list_windows_", e, "csv", sep=""), sep="\t", row.names = F, col.names = F, quote = F)
}
