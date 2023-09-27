#README

This folder contains all scripts used to run diploS/HIC (#1-13). They were adapted from Manthey et al. 2021 (https://doi.org/10.1093/gbe/evab120). I recomend checking and citing the original source (github.com/jdmanthey/certhia_genomes1).
In this folder you will also find an additional set of scripts (#13b-15) to select genomic windows classified under distinct models, and covert them into fasta alignemnts to run the model selection approach implemented with PipeMaster.

All priors and parameter convertions are described in discoal_notes_SPECIES.txt

01 - Runs the discoal simulations
02 - Concatenate all the output 
03 - Get 2000 of the simulations for each scenario
03b - Change the header so it is clear there were 2000 simulations 
4 - Zip the simulations for use in diploS/HIC.
5 - Applies a genomic mask to the analyzis
6 - Calculate summary statistics and convert them into feature vectors
7 - Organize the data for trainning the neural network
8 - Diplos/HIC trainning
10 - Create cluster submission jobs for turning the vcf files into feature vectors (20 kbp windows)
11 - Predict models for observed data.
12 - Concatenate predictions
13 - Plotting the results
13b - listing windows by model
14 - building fasta alignemts for Pipemaster
15 - loading Keras model (Pipemaster approach) in R and predicting models and parameter for windows.
16 - Combining results