#!/bin/bash
#this is to subset the 10 kb loci every 5 windows to reduce the size of the dataset and LD between loci
module load parallel-20171122 

mkdir loci_pipemaster

##sample the list of 10kb windows every 5 so every window will be 40kb apart
awk 'NR % 5 == 0' 1_10kb_windows.txt | perl -p -w -e 's/ (\w+\n)/-\1/g; s/ (\w+-)/_\1/g' > loci_pipemaster/1_loci_list.txt

cat loci_pipemaster/1_loci_list.txt | parallel -j 10 cp loci/{}.fasta loci_pipemaster/

screen 
qsub -I ../../../1_PBS_Processing_genomes.sh
cd ~/nas4/phleg/11_simulations/loci_pipemaster/
conda activate phyluce
mkdir loci_w_outgroup
mv *.fasta loci_w_outgroup/
phyluce_align_extract_taxa_from_alignments --alignments loci_w_outgroup/ --output loci_no_outgroup --input-format fasta --output-format fasta --cores 20 --exclude rgym-A11932_2 rgym-A11932_1


module load R-3.6.3 
R

#This goes in R
library(PipeMaster)
library(parallel)
load("../1_models.Rdata.RData")
pop.assign <- read.delim("../pop_list.txt", header = FALSE)
m1 <- get.data.structure(m1,path.to.fasta = "loci_no_outgroup", pop.assign=pop.assign)
stat <- obs.sumstat.ngs(model = m1, path.to.fasta = "loci_no_outgroup/" , pop.assign = pop.assign, moments = T)
write.table(stat, file = "stat.txt")