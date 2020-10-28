require(parallel)
require(pbapply)
source('auxfuns/simudat.R')
## install SASOM from local
install.packages("auxfuns/SASOM_0.1.0.tar.gz", repos = NULL, type ="source")
library(SASOM)

## load seeds
seeds = readRDS('../data/seed1M.rds')

## a function for parallel computing
pcf <- function(i,seed,n,paras){
    dat = simudat(seed[i],n,paras$W,paras$alpha,paras$theta,paras$xi,paras$rlevel)
    out = SASOM(dat$y,dat$X,dat$G,dat$W,"all")
    return(out[3:5])
}


##############################
######## Type 1 Error ########
##############################
## processor: 3.1 GHz Dual-Core Intel Core i5
## number of cores: 4
## approximate time required for each secnario 
## p = 10, n = 200 -- 6 hrs
## p = 10, n = 300 -- 14 hrs
## p = 10, n = 500 -- 33 hrs
## p = 10, n = 1000 -- 108 hrs
## p = 20, n = 200 -- 8 hrs
## p = 20, n = 300 -- 15 hrs
## p = 20, n = 500 -- 43 hrs
## p = 20, n = 1000 -- 132 hrs

B = 10^6
t1e = NULL
a = 0.001
for(p in c(10,20)){
    paras = readRDS(paste0('../data/S0_p',p,'.rds'))
    for(n in c(200,300,500,1000)){
        out = do.call(rbind,pblapply(1:B,pcf,seeds,n,paras))
        t1e = cbind(t1e,colMeans(out< a))
    }
}

colnames(t1e) = c('p = 10, n = 200','p = 10, n = 300','p = 10, n = 500','p = 10, n = 1000',
                  'p = 20, n = 200','p = 20, n = 300','p = 20, n = 500','p = 20, n = 1000')
rownames(t1e) = c('SASOM-F','SASOM-T','SASOM-D')
t1e


##############################
########### Power ############
##############################
## s = 1 -- Scenario I
## s = 2 -- Scenario II
## s = 3 -- Scenario III
## s = 4 -- Supplementary SI
## s = 5 -- Supplementary SII
## processor: 3.1 GHz Dual-Core Intel Core i5
## number of cores: 4
## approximate time required for each secnario 
## n = 250 -- 7 mins
## n = 300 -- 9 mins

s = 1 # change the value of s to obtain power table for each scenario
B = 10^2
a = 0.05
power = NULL
for(p in c(10,20)){
    paras = readRDS(paste0('../data/S',s,'_p',p,'.rds'))
    for(n in c(250,300)){
        out = do.call(rbind,pblapply(1:B,pcf,seeds,n,paras))
        power = cbind(power,colMeans(out< a))
    }
}

colnames(power) = c('p = 10, n = 250','p = 10, n = 300','p = 20, n = 250','p = 20, n = 300')
rownames(power) = c('SASOM-F','SASOM-T','SASOM-D')
power

