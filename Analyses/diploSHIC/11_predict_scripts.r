species = "lipau"
setwd(paste("~/Dropbox/1_Chapman_fellowship/data/DiploShic/", species, "/observed/", sep=""))


if (species=="phleg") {
  fasta_index <- read.table("../phleg/rhegma_ragout_complete.fasta.fai", sep="\t", stringsAsFactors=F)
  diplo_vcf_list1 <- fasta_index[grep("chr", fasta_index$V1), 1]
  vcf="tap.recode.vcf"
  }

if (species=="xipho") {
  fasta_index <- read.table("../pseudochromosomes.fasta.fai", sep="\t", stringsAsFactors=F)
  fasta_index$V1 <- gsub("chr", "Pseudo", fasta_index$V1)
  diplo_vcf_list1 <- fasta_index[grep("Pseudo", fasta_index$V1), 1]
  vcf="tap.recode.vcf"
}

if (species=="lipau") {
  fasta_index <- read.table("../pseudochromosomes.fasta.fai", sep="\t", stringsAsFactors=F)
  fasta_index$V1 <- gsub("chr", "Pseudo", fasta_index$V1)
  diplo_vcf_list1 <- fasta_index[grep("Pseudo", fasta_index$V1), 1]
  vcf="tap.recode.vcf"
}



a=1
a_script_name <- "1_running_predictions.sh"
write("", file= a_script_name, ncolumns=1, append=T)
write("conda activate selection3", file= a_script_name, ncolumns=1, append=T)
# loop for each chromosome in the analysis
for(a in 1:length(diplo_vcf_list1)) {
	
	# albescens
	

	a_vcf <- vcf
	a_diplo_name <- diplo_vcf_list1[a]
	a_chrom <- diplo_vcf_list1[a]
	a_length <- fasta_index[fasta_index[,1] == a_chrom,2]
	a_fvec <- paste(a_diplo_name, ".fvec", sep="")
	a_output_name <- paste(a_chrom, ".predictions_20", sep="")
	
	write("", file= a_script_name, ncolumns=1, append=T)
	write("python ../../diploSHIC/diploSHIC.py predict \\", file= a_script_name, ncolumns=1, append=T)
	
	
	###Change here for the best model of the 5 replicates
	if (species=="phleg") {
	write("../phlegModel3.json ../phlegModel3.weights.hdf5 \\", file= a_script_name, ncolumns=1, append=T)
	}
	if (species=="xipho") {
	  write("../xiphoModel3.json ../xiphoModel3.weights.hdf5 \\", file= a_script_name, ncolumns=1, append=T)
	}
	# if (species=="lipau") {
	#   write("../lipauModel1.json ../lipauModel1.weights.hdf5 \\", file= a_script_name, ncolumns=1, append=T)
	# }
	
	if (species=="lipau") {
	  write("../../phleg/phlegModel3.json ../../phleg/phlegModel3.weights.hdf5 \\", file= a_script_name, ncolumns=1, append=T)
	}
	write(paste(a_fvec, " ",  a_output_name, sep=""), file= a_script_name, ncolumns=1, append=T)
	write("", file= a_script_name, ncolumns=1, append=T)

}

###execute the .sh script manually



