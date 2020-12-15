###################################################################################
## This code perform SASOM and competing method to SKCM data.
## The preparation of the dataset from raw data is given in data_preparation.R
###################################################################################

## load required functions and packages
source('./auxfuns/pkg.R')
source('./auxfuns/competing_methods.R')
## install SASOM from local
install.packages("./auxfuns/SASOM_0.1.0.tar.gz", repos = NULL, type ="source")
library(SASOM)

## load data
dat = readRDS("../data/SKCM_data_for_H_pathway.rds")
n = length(dat)

########### Burden ############
res_burden = NULL
for(i in 1:n)
    res_burden = rbind(res_burden, burden(dat[[i]]$y,dat[[i]]$X,dat[[i]]$G))

########### uSKAT ############
res_uSKAT = NULL
for(i in 1:n)
    res_uSKAT = rbind(res_uSKAT, uSKAT(dat[[i]]$y,dat[[i]]$X,dat[[i]]$G))

########### DKAT ############
res_dkat = NULL
for(i in 1:n)
    res_dkat = rbind(res_dkat, DKAT(dat[[i]]$y,dat[[i]]$X,dat[[i]]$G,3))

########### mRand ############
res_mRand = NULL
for(i in 1:n)
    res_mRand = rbind(res_mRand, mRand(dat[[i]]$y,dat[[i]]$X,dat[[i]]$G))

########### uMiST ############
res_uMiST = NULL
for(i in 1:n)
    res_uMiST = rbind(res_uMiST,uMiST(dat[[i]]$y,dat[[i]]$X,dat[[i]]$G,dat[[i]]$W))

########### SASOM ############
res_sasom = NULL
for(i in 1:n)
    res_sasom = rbind(res_sasom, SASOM(dat[[i]]$y,dat[[i]]$X,dat[[i]]$G,dat[[i]]$W,"all"))

########## summarize results ###########
res = cbind(res_burden,res_uSKAT,res_dkat,res_mRand,res_uMiST,res_sasom[,3:5])
colnames(res) = c('Burden','uSKAT','DKAT-I','DKAT-T','mRand','uMiST','SASOM-F','SASOM-T','SASOM-D')
rownames(res) = names(dat)





