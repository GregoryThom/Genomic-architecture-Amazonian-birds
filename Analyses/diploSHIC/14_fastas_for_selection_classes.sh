#!/bin/bash

#This should be run on the server
#first make a folder and drop the lists for each class
module load parallel-20171122

species=phleg

ls *.csv | sed 's/.csv//g' | parallel "mkdir {} | mv {}.csv {}/"


if [ $species == phleg ]; then
ref=../rhegma_ragout_complete.fasta
elif [ $species == xipho ]; then
ref=../Xiphorhynchus_elegans_OUT-0059_pseudo_chrom_only_renamed.fasta
else
ref=../cepha
fi

#qsub -I 1_PBS_sampling_windows.sh
#cd ~/nas4/xipho/11_simulations/
#module load parallel-20171122 

#for i in $(ls); do
#mkdir $i/1_windows_10kb.vcf
#done

for i in $(ls); do
mkdir $i/1_vcf
done

for i in $(ls); do
#to create 100k vcf files
cat $i/$i.csv | perl -p -w -e "s/^/vcftools --gzvcf \.\.\/SNP-only_phased.vcf.gz --chr /g; s/\t(\w+\t\w+\n)/ --from-bp \1/g; s/\t(\w+\n)/ --to-bp \1/g; s/(.*--chr )(.*)( --from-bp )(\w+)( --to-bp )(\w+)/\1\2\3\4\5\6 --recode --out $i\/1_vcf\/\2_\4-\6/g" | parallel -j 16 {}
done


###now create the reference for each loci using window coordinates as a bed file to clip those 
#regions from the reference

#i=1_list_windows_prob.hard

for i in $(ls --hide=ind.txt); do
mkdir $i/refs_10K
mkdir $i/loci_10k
module load bedtools-2.26.0

bedtools getfasta -fi $ref -bed $i/$i.csv -fo $i/refs_10K/all_loci.fasta
grep \> $i/refs_10K/all_loci.fasta | parallel "grep -A1 {} $i/refs_10K/all_loci.fasta > $i/refs_10K/{}.fasta"
rm $i/refs_10K/all_loci.fasta
rename 's/>//g' $i/refs_10K/*
rename 's/:/_/g' $i/refs_10K/*
module load bcftools-1.10.2
ls $i/1_vcf/*.vcf | parallel bgzip {}
ls $i/1_vcf/ | parallel tabix $i/1_vcf/{}
perl -p -w -e 's/\t(\w+\n)/-\1/g; s/\t(\w+-)/_\1/g' $i/$i.csv > $i/windows_10kb.txt
# for some reason the getfasta is removing the first base pair of the reference, I'm adding 
# a N in the begining so the reference matches the vcf
perl -pi -w -e 's/^(\w)/N\1/g' $i/refs_10K/*
ls $i/1_vcf/*.gz | head -1 | parallel zcat {} | grep "#CHROM" | perl -p -w -e 's/.*FORMAT\t//g; s/\t/\n/g' > ind.txt


for f in $(cat $i/windows_10kb.txt); do cat ind.txt | parallel \
"bcftools consensus -H 1 -s {} -f $i/refs_10K/$f.fasta $i/1_vcf/$f.recode.vcf.gz | \
perl -p -w -e 's/>.*/>{}_1/g'" > $i/loci_10k/$f.fasta; cat ind.txt | parallel \
"bcftools consensus -H 2 -s {} -f $i/refs_10K/$f.fasta $i/1_vcf/$f.recode.vcf.gz | \
perl -p -w -e 's/>.*/>{}_2/g'" >> $i/loci_10k/$f.fasta; done
done




###This is to remove individuals from alignemnts and only use the loci that are fully aligned
conda activate twist


if [ $species == phleg ]; then
out=$(echo rgym-A11932_2 rgym-A11932_1)
elif [ $species == xipho ]; then
out=$(echo rgym-A11932_2 rgym-A11932_1)
else
out=$(echo rgym-A11932_2 rgym-A11932_1)
fi


for i in $(ls --hide=ind.txt); do
mkdir $i/loci_pipemaster
ls $i/loci_10k | parallel -j 30 "python ~/nas4/Avian_Tree_of_Life/tools/AMAS/amas/AMAS.py remove -x $out -u fasta -g $i/loci_pipemaster/ -c 30 -i $i/loci_10k/{} -f fasta -d dna -e"
rename 's/-out.fas//g' $i/loci_pipemaster/*
done


