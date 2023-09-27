species = "lipau"
setwd(paste("~/Dropbox/1_Chapman_fellowship/data/DiploShic/", species, "/observed/", sep=""))


options(scipen=999)
require(stats)


data <- read.table("predictions_20kbp.txt", sep="\t", header=T, stringsAsFactors=F)

if (species == "phleg") {
  data$chrom <- gsub("\\.[1-3]*", "", data$chrom)
  data <- data[grep("chrW|chr5ch1", data$chrom, invert = T),]
  data$chrom_num <- gsub("chr|LGE", "", data$chrom)
  data$chrom_num <- gsub("A", ".1", data$chrom_num)
  data$chrom_num <- gsub("Z", "29", data$chrom_num)
  data$chrom_num <- as.numeric(data$chrom_num)
  data <- data[order(data$chrom_num),]
} else  { 
  data$chrom <- gsub("Pseudo", "chr", data$chrom)
  data <- data[grep("chrLG2", data$chrom, invert = T),]
  data <- data[grep("chr1B", data$chrom, invert = T),]
  #data <- data[data$chrom != "chr1",]
  data$chrom_num <- gsub("chr|LGE", "", data$chrom)
  data$chrom_num <- gsub("A", ".1", data$chrom_num)
  #data$chrom_num <- gsub("4A", "4.1", data$chrom_num)
  #data$chrom_num <- gsub("1A", "1", data$chrom_num)
  data$chrom_num <- gsub("Z", "29", data$chrom_num)
  data$chrom_num <- as.numeric(data$chrom_num)
  data <- data[order(data$chrom_num),]
}

x1 <- data

for (i in 1:nrow(x1)) {
d <- max(x1[i,6:10])
if (d < 0.6) {
  x1$predClass[i] <- "low.prob"
}
}


# window size:
kbp_windows <- 20000

# chroms 
chroms <- unique(x1[,1])

# rename rownames
rownames(x1) <- seq(from=1, to=nrow(x1), by=1)


# add to counts of the window locations for plotting
count <- 0
for(a in 1:length(chroms)) {
	# length of this chromosome 
	a_rep_max <- max(x1[x1[,1] == chroms[a],3])
	
	# add the distance already covered by count
	x1[x1[,1] == chroms[a], 2] <- x1[x1[,1] == chroms[a], 2] + count
	x1[x1[,1] == chroms[a], 3] <- x1[x1[,1] == chroms[a], 3] + count
	
	
	# add length of this chromosome to the count
	count <- count + a_rep_max
}

#check that count addition worked
tail(x1)



# plot pie charts
par(mfrow=c(1,1))
par(mar=c(1,1,1,1))
# pie(table(x2$predClass), labels=NA)
pie(table(x1$predClass))


# how many windows are neutral in both lineages?
neutral1 <- x1[x1$predClass == "neutral", ]
# neutral2 <- x2[x2$predClass == "neutral", ]
locations_neutral1 <- paste(neutral1[,1], neutral1[,2], sep="_")
# locations_neutral2 <- paste(neutral2[,1], neutral2[,2], sep="_")
# table(locations_neutral1 %in% locations_neutral2)




# plot neutrality and classifications


five_class_colors <- c("purple", "lightblue", "lightcoral", "palegreen4", "red3", "darkgray")


window_size <- 10


# what are the unique chromosomes and their bounding areas for plotting?
chr <- unique(x1[,1])
chr_polygons_1 <- list()
chr_polygons_2 <- list()
# make the plotting polygons
for(a in 1:length(chr)) {
	a1 <- rownames(x1)[x1[,1] == chr[a]]
	a2 <- a1[length(a1)]
	a1 <- a1[1]
	chr_polygons_1[[a]] <- rbind(c(a1, 0), c(a2, 0), c(a2, 1), c(a1, 1), c(a1, 0))
	# a1 <- rownames(x2)[x2[,1] == chr[a]]
	# a2 <- a1[length(a1)]
	# a1 <- a1[1]
	# chr_polygons_2[[a]] <- rbind(c(a1, 0), c(a2, 0), c(a2, 1), c(a1, 1), c(a1, 0))
}




pdf("1_main_fig.pdf", width = 25, height = 10)
# set up plotting dimensions
par(mfrow=c(2,1))
par(mar=c(0.5,5,1,0))


# plot neutrality 
plot(c(-1,-1), ylim=c(0,1), xlim=c(1, nrow(x1)), xaxt="n", col="white", bty="n", cex.axis=1.1, cex.lab=1.3, ylab="Prob. Neutral")
odd <- 0
for(a in 1:length(chr_polygons_1)) {
	if(odd == 1) {
		polygon(chr_polygons_1[[a]], col="gray", border="white")
		odd <- 0	
	} else {
		odd <- 1
	}
}
# plot
#points(rownames(x1), x1$prob.neutral., pch=19, cex=0.2, col="steelblue3")	

# sliding windows
total_rep <- c()
place_rep <- c()
sliding_windows <- ceiling(as.numeric(rownames(x1)[nrow(x1)]) / window_size)
location <- window_size/2
for(b in 1:sliding_windows) {
	
	# get center location
	b_rep <- x1[location,2]
	# take surrounding lines in matrix
	lines_to_take <- (location-(window_size/2 - 1)):(location+(window_size/2))
	lines_to_take <- lines_to_take[lines_to_take >= 1 & lines_to_take <= nrow(x1)]
	b_rep2 <- x1[lines_to_take,]
	# see if they are within the needed window size
	b_rep2 <- b_rep2[b_rep2[,2] >= (b_rep - kbp_windows * (window_size/2)) & b_rep2[,2] <= (b_rep + kbp_windows * (window_size/2)),]
	
	# get values of window
	total_rep <- c(total_rep, mean(b_rep2$prob.neutral.))
	# get value for x axis
	place_rep <- c(place_rep, location)
	# add to row locator
	location <- location + window_size

}
lines(place_rep, total_rep, lwd=0.6, col="blue")

# plot classification of window
plot(c(-1,-1), ylim=c(0,1.2), xlim=c(1, max(x1[,3])), xaxt="n", yaxt="n", col="white", bty="n", cex.axis=1.1, cex.lab=1.3, ylab="Window Classification")


for(a in 1:nrow(x1)) {
	if(x1$predClass[a] == "neutral") {
		plot_color <- five_class_colors[1]
		line_points <- rbind(c(x1$classifiedWinStart[a], 0), c(x1$classifiedWinStart[a], 0.2))
	} else if(x1$predClass[a] == "linkedSoft") {
		plot_color <- five_class_colors[2]
		line_points <- rbind(c(x1$classifiedWinStart[a], 0.2), c(x1$classifiedWinStart[a], 0.4))
	} else if(x1$predClass[a] == "linkedHard") {
		plot_color <- five_class_colors[3]
		line_points <- rbind(c(x1$classifiedWinStart[a], 0.6), c(x1$classifiedWinStart[a], 0.8))
	} else if(x1$predClass[a] == "soft") {
		plot_color <- five_class_colors[4]
		line_points <- rbind(c(x1$classifiedWinStart[a], 0.4), c(x1$classifiedWinStart[a], 0.6))
	} else if(x1$predClass[a] == "hard"){
		plot_color <- five_class_colors[5]
		line_points <- rbind(c(x1$classifiedWinStart[a], 0.8), c(x1$classifiedWinStart[a], 1))
	} #else if(x1$predClass[a] == "low.prob"){
	  # plot_color <- five_class_colors[6]
	  # line_points <- rbind(c(x1$classifiedWinStart[a], 1.0), c(x1$classifiedWinStart[a], 1.2))
	# }
	
	lines(line_points, col=plot_color, lwd=0.1)
}
dev.off()


probs=NULL
neutral <- c(nrow(x1[x1$predClass=="neutral",])/nrow(x1))*100
linkedSoft <- c(nrow(x1[x1$predClass=="linkedSoft",])/nrow(x1))*100
soft <- c(nrow(x1[x1$predClass=="soft",])/nrow(x1))*100
linkedHard <- c(nrow(x1[x1$predClass=="linkedHard",])/nrow(x1))*100
hard <- c(nrow(x1[x1$predClass=="hard",])/nrow(x1))*100
lowPprob <- c(nrow(x1[x1$predClass=="low.prob",])/nrow(x1))*100


probs <- cbind(neutral, linkedSoft, soft, linkedHard, hard, lowPprob)
write.csv(probs, "1_proportion_of_classified_windows.csv")




