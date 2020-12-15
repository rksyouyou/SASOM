########################################################################
## The raw RNAseq, phenotype, and somatic mutation datasets were download from https://xenabrowser.net/datapages/?cohort=GDC%20TCGA%20Liver%20Cancer%20(LIHC)&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443
########################################################################

### load required packages 
require(qusage)
require(reshape2)

########################################################################
########################### Covariate data #############################
########################################################################
## the raw data contains 477 subjects and 2 columns (subject name and age)
pdat = read.csv('../data/raw/TCGA-LIHC.GDC_phenotype.tsv',sep = '\t',header=TRUE,fill=TRUE)[,c('submitter_id.samples','age_at_initial_pathologic_diagnosis')] 
colnames(pdat) = c('sample','age')
pdat = pdat[-which(is.na(pdat[,2])),]

########################################################################
############################ Somatic data ##############################
########################################################################
## somatic data: each row corresponds to a somatic mutation in a sample
sdat = read.csv('../data/raw/TCGA-LIHC.muse_snv.tsv',fill=TRUE,header=TRUE,sep='\t',check.names = FALSE)
## mutation effect serverity map
effect_tbl = readRDS('../data/raw/effect_table.rds')
## in somatic data, mutations with effect belong to grp 0 (according to effect_tbl) were removed
sdat = sdat[-which(sdat$effect%in%rownames(effect_tbl)[which(effect_tbl==0)]),]
## re-organize somatic data from vector format to matrix format, each row represent a subject, each columnrepresents a mutation (defined by location). Given a subject and mutation, if the mutation occures on the subject, the corresponding element value is 1 and 0 otherwise. 
tmp = cbind(sdat[,c(1,2,4,5)],val=1)
smat = dcast(tmp,Sample_ID~gene+start+end,value.var = "val",fill=0) 
rownames(smat) = smat[,1]
smat = as.matrix(smat[,-1])

########################################################################
################## Identify cancer subtypes (response) #################
########################################################################
## conversion map between ensembl gene ID and gene symbol
gmap = read.csv('../data/raw/gencode.v22.annotation.gene.probeMap',fill=TRUE,header=TRUE,sep='\t',check.names = FALSE,row.names=1)
## RNAseq data: note that the RNAseq data is too big to upload to GitHub, this data will be download from the website. This will take a few minites. 
con <- gzcon(url("https://gdc.xenahubs.net/download/TCGA-LIHC.htseq_counts.tsv.gz"))
txt <- readLines(con)
gdat <- read.csv(textConnection(txt),fill=TRUE,header=TRUE,sep='\t',check.names = FALSE,row.names=1)
## pre-specified genes that related to immune response
ref_gene = c("BLK", "CD19", "FCRL2", "MS4A1", "KIAA0125", "TNFRSF17", "TCL1A", "SPIB", "PNOC")
## subset from the raw data and formating 
gdat1 = t(gdat[which(rownames(gdat)%in%rownames(gmap)[(gmap[,1]%in%ref_gene)]),])
gmap1 = gmap[(gmap[,1]%in%ref_gene),]
gmap1 = gmap1[match(colnames(gdat1),rownames(gmap1)),] # match order 
colnames(gdat1) = gmap1$gene
## hierarchical clustering
set.seed(500)
hc = hclust(dist(gdat1),method='ward.D') 
ymat = as.matrix(cutree(hc,3))

########################################################################
##################### mutation effect map  #############################
########################################################################
mmap = data.frame(mut = paste0(sdat$gene, '_', sdat$start, '_', sdat$end),gene = sdat$gene,effect = sdat$effect)
mmap = mmap[!duplicated(mmap[, 1]),] # remove duplicated rows 
mmap = mmap[mmap[, 1] %in% colnames(smat),] # subsample 
mmap = mmap[match(colnames(smat), mmap[, 1]),] # match the order
mmap = cbind(mmap,grp = effect_tbl[match(mmap$effect,rownames(effect_tbl)),])


########################################################################
########################## data formatting  ############################
########################################################################
## load pathway: each list contains a set of gene names that belong to the respective pathway 
mdat = read.gmt('../data/raw/h.all.v7.0.symbols.gmt')
npath = length(mdat)
## subject with complete information (covariate, somatic mutation and response) were used 
sid = intersect(intersect(rownames(ymat), pdat[, 1]), rownames(smat))

## covariante matrix (X): same for all pathways 
X = pdat[pdat$sample %in% sid,]
X = X[match(sid, X$sample),] # match order 
Xmat = model.matrix(~ X[, 2])
colnames(Xmat) = c('Intercept', 'age')
rownames(Xmat) = X[, 1]

## construct genotype/somatic mutation  matrix (G) and variate character matrix (W)
gwlist = list() # a list of length npath, each list contains G and M matrices
thres = 2 # mutations with at least 3 mutants were preserved
for (i in 1:npath) {
    ## genotype matrix
    G = smat[, which(mmap[, 2] %in% mdat[[i]])]
    rownames(G) = rownames(smat)
    G = G[match(sid, rownames(G)), ]
    if (sum(colSums(G) >= thres) > 0) {
        G = G[, colSums(G) >= thres, drop = FALSE]
    } else G = NULL
    if(!is.null(G))
        if(length(unique(do.call(rbind,strsplit(colnames(G),'_'))[,1]))<5)
            G = NULL
    ## pathways that contain at least ten genes with somatic mutations were analyzed
    if(!is.null(G)){
        ## character variant matrix (W)
        tmp = mmap[mmap[, 1] %in% colnames(G), ]
        tmp = tmp[match(colnames(G), tmp$mut), ]
        suppressMessages({W = dcast(tmp, mut ~ grp, length)})
        W = W[match(colnames(G), W[,1]), , drop = FALSE] # match order 
        rownames(W) = W[, 1]
        W = W[,-1,drop=FALSE]
        W[,1] = 1
        W = t(W)
    } else W = NULL
    gwlist[[i]] = list(G = G,W = W)
}

## response value 
ymat = ymat[rownames(ymat)%in%sid,,drop=FALSE]
ymat = ymat[match(sid,rownames(ymat)),]

## data wrap up
dlist = list()
dname = NULL
i = 0
j = 1
while(i < npath){
    i = i + 1
    if(!is.null(gwlist[[i]]$G)){
        dlist[[j]] = list(y = ymat,X = Xmat,G = gwlist[[i]]$G,W = gwlist[[i]]$W)
        dname = c(dname,names(mdat)[i])
        j = j+1
    }
}
names(dlist) = dname
## saveRDS(dlist,'../data/LIHC_data_for_H_pathway.rds')




