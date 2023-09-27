library(ggplot2)
source("~/Dropbox/1_Chapman_fellowship/data/scripts/analyses/plot_twisst.R")
###Running Astral
species="xipho"
snps="100"


setwd(paste("~/Dropbox/1_Chapman_fellowship/data/10_astral_new/", species, sep = ""))
vars <- list.dirs(paste("../../2_windows_new/", species, "/", snps, "/", sep = ""), full.names = F)[-1]
#var=vars[2]
##all was too big --->>> shuf -n 50000 all_trees/all_gentrees.txt > all_trees/all_gentrees_reduced.txt
astral <- function(var){
  nam <- gsub("_trees", "", var)
  if (file.exists(paste(nam, "species.tree", sep = ""))==F) {
    
  system(paste("java -Xmx10G -jar  ~/Astral/astral.5.7.3.jar -i ../../2_windows_new/",
                      species, "/", snps, "/", var, "/", nam,
                      "_gentrees.txt -q st1_to_score.tree -t 2 -a pop_list.txt -o ",
                      nam, "species.tree 2>", var, ".out", sep = ""), intern = T )
  }

out <- read.table(paste(nam, "species.tree", sep = ""), header = F, sep = "\t")
out <- gsub("\\(.*\\[(.*)\\].*\\[.*", "\\1", as.character(out[1,]))
#out <- gsub("\\].*", "", out)
out <- strsplit(out, ";")
names(out) <- nam
return(out)
}

system.time({
  all <- lapply(vars, astral)
})



all2 <- t(do.call(cbind.data.frame, all))
colnames(all2) <- gsub("=.*", "", all2[1,])
all2 <- as.data.frame(gsub(".*=", "", all2))
i <- c(1:ncol(all2)) 
all2[ , i] <- apply(all2[ , i], 2,          
                    function(x) as.numeric(as.character(x)))
all2$var <- rownames(all2)


data=NULL
#i="pp1"
for (i in c("pp1", "pp2", "pp3")) {
 d <- grep(i , colnames(all2))
   s <- as.data.frame(all2[,d])
   colnames(s) <- "pp"
   s$var <- i
   s$var2 <- row.names(all2)
   data <- rbind(data, s)
}

data$var <- gsub("pp2", "topo3", data$var)
data$var <- gsub("pp3", "topo2", data$var)
data$var <- gsub("pp1", "topo1", data$var)
#data$var2 <- rep(c("1_all", "5-dxy10", "5-dxy90", "8-ABBA10", "8-ABBA90", "3-Fst10", "3-Fst90",
#                   "6-mtDNA", "4-Pi10", "4-Pi90", "2-recom10", "2-recom90", "7-UCE"), 3)

data1 <- data[grep("chr|all", data$var2),]
data1$V3 <- as.numeric(data1$V3)
data1$var2 <- gsub("all", "Whole genome", data1$var2)
data1 <- data1[order(data1$V3),]
data1$var2 <- factor(data1$var2, levels = unique(data1$var2))


if (species == "phleg") {
  fg <- "Phlegopsis nigromaculata" 
} 
if (species == "xipho") {
  fg <- "Xiphorhynchus spixii"
} 
if (species == "lipau") {
  fg <- "Lipaugus vociferans"
}


p<-ggplot(data1, aes(x=var2, y=pp, group=var)) +
  geom_point(aes(color=var), size=6, position=position_jitter(w=0.1), alpha=0.9) +
  theme_minimal() + scale_color_manual(values = topo_cols) + 
  geom_hline(yintercept=0.95, linetype="dashed", color = "red") +
  theme_minimal(base_size = 21) + 
  theme(axis.text.x = element_text(angle = 45, size = 12, hjust = 1), 
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), legend.position = "none") +
  ggtitle(fg) + labs(y="Posterior probability", x = "")

p

pdf("1_ST_support_CHR.pdf")
p
dev.off()


#####
data2<- data
if (data2$pp[1] ==0 ) {
  data2$pp[1] <- 0.00001
}
#data2$pp <- as.character(data2$pp)
data2 <- data2[-grep("chr", data2$var2),]

data2 <- data2[order(data2$var, data2$var3),]
data2$var2 <- rep(c("Whole genome", "Low recombination", "High recombination", "Low Fst", "High Fst",
"Low pi", "High pi", "Low Dxy", "High Dxy", "mtDNA", "UCE", "Low fdM", "High fdM"), 3)
data2$var2 <- factor(data2$var2, levels = unique(data2$var2))


p2<-ggplot(data2, aes(x=var2, y=pp)) +
  geom_point(aes(color=var), size=6, position=position_jitter(w=0.1), alpha=0.9) +
  scale_color_manual(values = topo_cols) + 
  geom_hline(yintercept=0.95, linetype="dashed", color = "red") +
  theme_minimal(base_size = 21) + 
  theme(axis.text.x = element_text(angle = 45, size = 15, hjust = 1), 
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), legend.position = "none") +
  ggtitle(fg) + labs(y="Posterior probability", x = "")

p2

pdf("1_ST_support_stats.pdf")
p2
dev.off()

