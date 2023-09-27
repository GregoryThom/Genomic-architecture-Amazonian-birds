#!/bin/bash

#this is to create fasta alignments for 10kbp or 100kbp windows from whole genome vcf files

###this is for sampling 100kbp. I'm using the same windows as the ones used for summary statistics using https://github.com/simonhmartin/genomics_general 
zcat ../windows/window_stats_100k.csv.gz | awk -F',' '{print $2,$3,$4}' > 1_100kb_windows.txt

##This for sampling 10kb every 100kb
zcat ../windows/window_stats_10k.csv.gz | awk -F',' '{print $2,$3,$4}' | awk 'NR % 10 == 0' | grep -v Z | grep -v W > 1_10kb_windows.txt


screen 
qsub -I 1_PBS_sampling_windows.sh
cd ~/nas4/xipho/11_simulations/
module load parallel-20171122 

#to create 100k vcf files
mkdir 1_windows_100kb.vcf
cat 1_100kb_windows.txt | perl -p -w -e 's/^/vcftools --gzvcf SNP-only_phased.vcf.gz --chr /g; s/ (\w+ \w+\n)/ --from-bp \1/g; s/ (\w+\n)/ --to-bp \1/g; s/(.*--chr )(.*)( --from-bp )(\w+)( --to-bp )(\w+)/\1\2\3\4\5\6 --recode --out 1_windows_100kb.vcf\/\2_\4-\6/g' | parallel -j 24 {}


###now create the reference for each loci using the 1_100kb_windows.txt as a bed file to clip those 
#regions from the reference

#reference file
ref=Cephalopterus_ornatus_B10K001_pseudochromossome_only.fasta

mkdir refs_100K
mkdir loci_100k
module load bedtools-2.26.0
perl -pi -w -e 's/ /\t/g' 1_100kb_windows.txt
sed -i '1d' 1_100kb_windows.txt

bedtools getfasta -fi $ref -bed 1_100kb_windows.txt -fo refs_100K/all_loci.fasta
grep \> refs_100K/all_loci.fasta | parallel "grep -A1 {} refs_100K/all_loci.fasta > refs_100K/{}.fasta"
rm refs_100K/all_loci.fasta
rename 's/>//g' refs_100K/*
rename 's/:/_/g' refs_100K/*
module load bcftools-1.10.2
ls 1_windows_100kb.vcf/*.vcf | parallel bgzip {}
ls 1_windows_100kb.vcf/ | parallel tabix 1_windows_100kb.vcf/{}
perl -p -w -e 's/\t(\w+\n)/-\1/g; s/\t(\w+-)/_\1/g' 1_100kb_windows.txt > windows_100kb.txt
# for some reason the getfasta is removing the first base pair of the reference, I'm adding 
# a N in the begining so the reference matches the vcf
perl -pi -w -e 's/^(\w)/N\1/g' refs_10K/*
ls 1_windows_100kb.vcf/*.gz | head -1 | parallel zcat {} | grep "#CHROM" | perl -p -w -e 's/.*FORMAT\t//g; s/\t/\n/g' > ind.txt


for i in $(cat windows_100kb.txt); do cat ind.txt | parallel \
"bcftools consensus -H 1 -s {} -f refs_100K/$i.fasta 1_windows_100kb.vcf/$i.recode.vcf.gz | \
perl -p -w -e 's/>.*/>{}_1/g'" > loci_100k/$i.fasta; cat ind.txt | parallel \
"bcftools consensus -H 2 -s {} -f refs_100K/$i.fasta 1_windows_100kb.vcf/$i.recode.vcf.gz | \
perl -p -w -e 's/>.*/>{}_2/g'" >> loci_100k/$i.fasta; done
