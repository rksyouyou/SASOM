
simudat <- function(seed,n,W,alpha,theta,xi,rlevel){
    set.seed(seed)

    p = ncol(W)
    pvec_G = rep(0.05,p)
    X = cbind(1,replicate(1,rbinom(n,1,0.5)),replicate(1,rnorm(n)))
    colnames(X) = paste0('X',0:2)
    G <- sapply(pvec_G,function(x) rbinom(n,1,x)) # variant matrix n*p
    Gw <- G%*%t(W)
    eta <- X%*%alpha + Gw%*%theta + G%*%xi
    seta <- rowSums(exp(eta))+1
    if(rlevel==1) pmat <- cbind(1/seta,exp(eta[,1])/seta, exp(eta[,2])/seta)
    if(rlevel==2) pmat <- cbind(exp(eta[,1])/seta, 1/seta,exp(eta[,2])/seta) 
    if(rlevel==3) pmat <- cbind(exp(eta[,1])/seta, exp(eta[,2])/seta,1/seta) 
    ymat <- t(apply(pmat,1,function(x) rmultinom(1,1,x))) # response matrix

    return(list(y=ymat,X=X,G=G,W=W,rlevel=rlevel))
}



