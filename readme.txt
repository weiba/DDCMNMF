This file contain the original datasets(in data file) used in this research ,source codes for running the SNMNMF(in SNMNMF_source file), the preprocess steps and driver gene smoothing steps(in R script of preprocess_SNMNMF), final clustering consensus steps(in R script of clustering_consensus)

Running the DDCMNMF need to install the matlab and R environment in advance.

The whole process to employ the DDCMNMF can be divided into following three steps:
In the R environment
1. preprocess the data for running the SNMNMF:
a. set the right location of running dataset in preprocess_SNMF.R file, such as patMutMatrix=read.table('Lung(102)106/LUSC_Mutation.txt')
b. source('preprocess_SNMF.R')
c. run the main_function()

2.Run the SNMNMF model
 Use the results which generated from the first step as input to run the SNMNMF model under Matlab environment.
a. set the running location in the('SNMNMF_source/SNMNMF')  
b. load the results of last step. such as load('SNMNMF_source/InputData/Lung_SNMNMF.mat')
c.  run the Run_SNMNMF(Input)

3. final consensus clustering
In the R environment
a. set the right location of the W which generated from step 2 such as path='SNMF_source/SNMNMF/W/'. and the original mutation matrix such as patMutMatrix=read.table('data/Breast(104)105/BREAST_Mutation.txt',header = T)

b.source the clustering_consensus.R script

the output of final results is in the list of final_ours.
