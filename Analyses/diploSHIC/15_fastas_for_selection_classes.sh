
#For each folder


#cd loci_pipemaster/

rename 's/\./_/g' *fasta
rename 's/_fasta/.fasta/g' *fasta



for i in `seq 1 40`; do mkdir -p "folder$i"; find . -type f -maxdepth 1 | head -n 25 | xargs -i mv "{}" "folder$i"; done

#screen 
#cd nas4/phleg/11_simulations/
#qsub -I ../../1_PBS_Processing_genomes.sh
#cd nas4/phleg/11_simulations/
module load R-3.6.3 
R
library(PipeMaster)
library(parallel)
load("../../1_models.Rdata.RData")
pop.assign <- read.delim("../../pop_list.txt", header = FALSE)


numCores=40
folders=c(1:40)

mc_sum_stat <- function(fold){
stat <- obs.sumstat.ngs(model = m1, path.to.fasta = paste("folder", fold, sep ="") , pop.assign = pop.assign, moments = F )
}

 system.time({
  mclapply(folders, mc_sum_stat, mc.cores = numCores)
})
q(save="no")


#on system
mkdir loci_summary_stats/
seq 40 | parallel cp folder{}/*.out loci_summary_stats/