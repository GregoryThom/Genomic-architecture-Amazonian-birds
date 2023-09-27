species = "xipho"

setwd(paste("~/Dropbox/1_Chapman_fellowship/data/DiploShic/", species, "/observed/", sep =""))

if (species == "phleg") {
fasta_index <- read.table("../rhegma_ragout_complete.fasta.fai", sep="\t", stringsAsFactors=F)
diplo_vcf_list1 <- fasta_index[grep("chr", fasta_index$V1), 1]
}

if (species == "xipho") {
  fasta_index <- read.table("Xiphorhynchus_elegans_OUT-0059_pseudo_chrom_only.fasta.fai", sep="\t", stringsAsFactors=F)
diplo_vcf_list1 <- fasta_index[grep("Pseudo", fasta_index$V1), 1]
  }



vcf="tap.recode.vcf"

a=1
# loop for each chromosome in the analysis

for(a in 1:length(diplo_vcf_list1)) {
	
	# albescens
	a_vcf <- vcf
	a_diplo_name <- diplo_vcf_list1[a]
	a_output_name <- paste("a_", a_diplo_name, ".job", sep="")
	a_chrom <- diplo_vcf_list1[a]
	a_length <- fasta_index[fasta_index[,1] == a_chrom,2]
	a_fvec <- paste(a_diplo_name, ".fvec", sep="")
	
	# write
	write("#!/bin/sh", file=a_output_name, ncolumns=1)
	write("#PBS -l ncpus=1:mem=20gb", file=a_output_name, ncolumns=1, append=T)
	write("#PBS -l walltime=5000:00:00", file=a_output_name, ncolumns=1, append=T)
	write("#PBS -N discoal_sim", file=a_output_name, ncolumns=1, append=T)
	write(paste("#PBS -j oe", a_diplo_name, sep=""), file=a_output_name, ncolumns=1, append=T)
	write("#PBS -m ae", file=a_output_name, ncolumns=1, append=T)
	write("#$ -pe sm 4", file=a_output_name, ncolumns=1, append=T)
	write("#$ -P quanah", file=a_output_name, ncolumns=1, append=T)
	write("#PBS -k oe", file=a_output_name, ncolumns=1, append=T)
	write("", file=a_output_name, ncolumns=1, append=T)
	write("source activate selection2", file=a_output_name, ncolumns=1, append=T)
	write("", file=a_output_name, ncolumns=1, append=T)
	write("cd /home/gthomesilva/mendel-nas1/discoal/observed", file=a_output_name, ncolumns=1, append=T)
	write("", file=a_output_name, ncolumns=1, append=T)
	write("python ../diploSHIC/diploSHIC.py fvecVcf diploid \\", file=a_output_name, ncolumns=1, append=T)
	write(paste(a_vcf, " ", a_chrom, " ", a_length, " \\", sep=""), file=a_output_name, ncolumns=1, append=T)
	write(paste(a_fvec, " --targetPop tap --sampleToPopFileName pop_list.txt --winSize 220000 \\", sep=""), 
		file=a_output_name, ncolumns=1, append=T)
	write("--maskFileName diplo_mask.fasta.gz", file=a_output_name, ncolumns=1, append=T)
	write("", file=a_output_name, ncolumns=1, append=T)
}
