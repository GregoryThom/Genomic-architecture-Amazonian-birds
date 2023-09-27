# Genomic architecture drives population structuring in Amazonian birds
This repository contains:
Supplementary datasets: Tables S1, S3-S9, and S21.

Analyses: The scripts used to run phylogenetic and population genetic analyses

Data: (link to Dryad)

Description of the data available in Dryad:

Concatenated_matrix/ 
Whole-genome concatenated SNPs alignment in nexus format for Phlegopsis nigromaculata (P_nigromaculata.nexus.gz), Xiphorhynchus spixii (X_spixii.nexus.gz), and Lipaugus vociferans (L_vociferans.nexus.gz). “To estimate phylogenetic relationships between individuals, we estimated supermatrix trees concatenating all SNPs using IQTree2 (Minh et al. 2020). We converted vcf files to phylip format using vcf2phylip.py (Ortiz 2019), randomly resolving heterozygous genotypes, and keeping SNPs present in at least 80% of the individuals.”


diploSHIC/
Observed data, simulations, population assignment, and genomic mask used in the diploS/HIC  (Kern and Schrider 2018) analyses for Phlegopsis nigromaculata, Xiphorhynchus spixii, and Lipaugus vociferans. 


Simulations/
    Simulations_10kb/
    Genetic simulations for the three alternative topologies. We simulated data under three alternative topologies, matching the unrooted trees tested in our phylogenetic approach: topology 1) (out,(Belem,(Xingu,Tapajos))); topology 2) (out,(Tapajos,(Xingu,Belem))); and topology 3) (out,(Xingu,(Tapajos,Belem))). We simulated 5,000 loci of 10kb, using uniform and wide priors for all parameters (Table S19), and performed one million simulations per model.
    
Simulations_100kb/
    Genetic simulations for the two alternative topologies. We simulated data under two alternative topologies, matching the unrooted trees tested in our phylogenetic approach: 
    topology 1) (out,(Belem,(Xingu,Tapajos))); topology 2) (out,(Tapajos,(Xingu,Belem))).We performed 100,000 simulations per model, and used the same uniform priors for all parameters as implemented above.
    
P_nigromaculata/ 
    Observed data (fasta files) for Phlegopsis nigromaculata for the 10kb and 100kb approaches. 
    
X_spixii/ 
    Observed data (fasta files) for Xiphorhynchus spixii  for the 10kb and 100kb approaches. 

L_vociferans/ 
    Observed data (fasta files) for Lipaugus vociferans for the 10kb and 100kb approaches. 
    

Species_tree_Astral_inputs/ 
Gene trees for chromosomes and summary statistics partitions for the three studies species used in our species tree approach. “To estimate the posterior probability of unrooted species trees, we used Astral-III v5.1.1 (Zhang et al. 2018; Rabiee et al. 2019), using the gene trees produced with phyml as inputs.”

Twisst/
Inputs and outputs from topology weighting analyses for the three studies species. To estimate unrooted topology weight for each window across the genome, we used Twisst (Martin and Van Belleghem 2017). “*.vcf.gz” are phased VCF files, and “*.trees.gz” are input gene trees.

