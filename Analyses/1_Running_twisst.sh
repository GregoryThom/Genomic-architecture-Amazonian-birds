#Running twist
Please check the original repository for twisst: https://github.com/simonhmartin/twisst 

##compress and index filtered vcf (wo/ monomorphic sites)
bgzip input.vcf
tabix input.vcf.gz
conda activate twist

#phase the vcd
java -Xmx12g -jar ../../programs/beagle/beagle.27Jan18.7e1.jar gt=input.vcf.gz out=output.vcf.gz impute=true nthreads=20 window=10000 overlap=1000 gprobs=false

##Converting VCF into Geno

python /home/gthomesilva/nas4/programs/genomics_general-master/VCF_processing/parseVCF.py -i SNP-only_phased.vcf.gz --skipIndels | gzip > SNP_phased.geno.gz

#Neighbour joining trees for snp windows.

conda activate twist

export PYTHONPATH=$PYTHONPATH:/home/gthomesilva/nas4/programs/genomics_general-master/

python ../../programs/genomics_general-master/phylo/phyml_sliding_windows.py -T 16 -g SNP_phased.geno.gz --prefix 100SNPs/SNP-only_phased.phyml_bionj.w100 -w 100 --windType sites --model GTR --phyml /home/gthomesilva/nas4/programs/phyml-3.3.20190909/src/phyml --log phyml.log

#duplicate the number of individuals in pop_list (phased haplotypes)

#Running twists

python ../../programs/twisst-master/twisst.py -t 100SNPs/SNP-only_phased.phyml_bionj.w100.trees.gz -w 100SNPs/SNP_phased.weights.csv.gz --outputTopos topologies.trees -g bel -g tap -g xin -g out --groupsFile pop_list.txt 
