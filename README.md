# BayesC

## Description
Implementation of GIBBS sampler for Bayesian penalized regression with a spike and slab prior.

The model is 
``y=Xb+e``

where ``y`` is the vector of phenotypes, ``X`` is the genotype matrix, ``b`` is the vector with the SNP effect sizes and ``e`` is the error.


## Installation
To install simply type in the command line:    
make

and run with 

./BayesC

adding any of the input options that are given below.


## Requirements
To correctly compile the program you need to have installed the boost and Eigen libraries.   

In ubuntu linux you type the following:

sudo apt-get install libboost-all-dev     
sudo apt install libeigen3-dev

Alternatively you can download the headers and transfer them to your /usr/include/ folder.

## Input and options 

If an input is not provided then it will simulate the genotype and phenotype data and perform estimation using the simulated data.

Possible options:

--M : No. of markers    
--N : No. of individuals    
--iter : No. of Gibbs sampling iterations    
--pNZ : Proportion of nonzero markers assumed in simulation mode (when input file is not provided).   
--h2 : Assumed heritability (when input file is not provided).    
--input : input tag for genotype matrix file and for phenotype vector file. The genotype matrix is assumed to have the suffix ".X" after the tag provided with input. Similarly for the file containing the phenotypes of the individuals the suffix ".Y" is assumed. If an input is not provided then the program will simulate the data to do the inference.    
--output : Results of the gibbs sampler are outputed using the tag provided here with the ending "_estimates.txt". If the program has simulated the data it will also output the simulated genotype matrix with suffix "_simulated_X.txt", the phenotype data with suffix "_simulated_Y.txt" and the true effect sizes (betas) with suffix "_simulated_betatrue.txt"    
