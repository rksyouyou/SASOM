### main code for SASOM method
## load required packages for SASOM
require(nnet)
require(CompQuadForm) 


## input
## y: matrix or a factor vector. When a matrix is given, each row represents a subject, each column corresponds to a subtype. A subject must belong to one of the subtypes, i.e. each row sums up to 1. When a vector is given, each element corresponds to the assignment of subtype for each subject.  
## rlevel - reference level, default is the first level. rlevel  must be a number if y is a matrix or a factor when y is a vector.
## X - covariant matrix, each row represents a subject.
## G - variant matrix, each row represents a subject.
## W - variant character matrix, each row represents a genome character, each column represents a variant. The column number of W is equal to the row number of G. 

varfun <- function(mu){
    ## construct variance matrix for response
    mu.s <- as.vector(mu)
    n <- nrow(mu)
    s.J <- ncol(mu)
    D <- matrix(0,n*s.J,n*s.J)
    for(i in 1:(s.J-1))
        for(j in (i+1):s.J)
            diag(D[(n*(i-1)+1):(n*i),(n*(j-1)+1):(n*j)]) <- -mu[,i]*mu[,j]
    D <- D+t(D)
    diag(D)  <- mu.s*(1-mu.s)
    return(D)
}


SASOM <- function(y,X,G,W,rlevel=1){

    ## check the format of Y, if its a vector, convert it to matrix
    if(is.vector(y)|is.factor(y)) {
        y <- as.factor(y)
        J <- length(levels(y))
        n <- length(y)
        ymat <- matrix(0,n,J)
        ymat[cbind(1:n,as.numeric(y))] <- 1
        colnames(ymat) <- levels(y)
    } else if(is.matrix(y)) {
        if(any(rowSums(y)!=1)) stop("Error: summation of each row must equal to 1!")
        J <- ncol(y)
        n <- nrow(y)
        ymat <- y
        colnames(ymat) <- 1:J
    } else {
        stop("Error: y must be a vector or matrix!")
    }

    ## adjust the order of columns based on reference level
    ## put the column correspond to reference level as the first column
    rc <- which(colnames(ymat) == rlevel)
    ymat <- ymat[,c(rc,setdiff(1:J,rc))]


    ## check if there is intercept in X, if not add intercept column
    if(sum(X[,1] == 1)<n){
        print('An intercept column is added to covariant matrix.')
        X <- cbind(1,X)
    }

    ## check dimension of input matrices
    if(nrow(X)!=n) stop('Error: dimensions of covariant matrix and response are inconsistent!')
    if(nrow(G)!=n) stop('Error: dimensions of variant matrix and response are inconsistent!')
    if(ncol(G)!=ncol(W)) stop('Error: dimensions of variant and variant character matrices are inconsistent!')

    ## define new notations for futurn use
    q <- nrow(W)
    Gw <- G%*%t(W)
    Y.s <- matrix(as.vector(ymat[,-1]),ncol=1)
    X.s <- kronecker(diag(J-1),X)
    G.s <- kronecker(diag(J-1),G)
    S <- kronecker(diag(J-1),t(Gw))
    M <- cbind(X,Gw)
    M.s <- kronecker(diag(J-1),M)

    ## fit two models under H01 and H02
    fit1 <- multinom(ymat~0+X,trace=FALSE)
    fit2 <- multinom(ymat~0+X+Gw,trace=FALSE)
    mu1 <- fit1$fitted.values[,-1]
    mu2 <- fit2$fitted.values[,-1]
    mu1.s <- as.vector(mu1)
    mu2.s <- as.vector(mu2)

    ## construct variance matrix under H01 and H02
    D1 <- varfun(mu1)
    D2 <- varfun(mu2)

    ## score test for theta
    U.theta <- S%*%(Y.s-mu1.s) # (J-1)q*1
    DX <- D1%*%X.s
    V.theta <- S%*%(D1-DX%*%solve(t(X.s)%*%DX)%*%t(DX))%*%t(S)
    stat.theta <- t(U.theta)%*%solve(V.theta)%*%(U.theta)
    pval.theta <- pchisq(stat.theta,(J-1)*q,lower.tail=FALSE)

    ## score test for tau 
    U.tau <- t(G.s)%*%(Y.s-mu2.s)
    DM <- D2%*%M.s
    V.tau = NULL
    V.tau <- t(G.s)%*%(D2-DM%*%solve(t(M.s)%*%DM)%*%t(DM))%*%G.s
    stat.tau <- sum(U.tau^2)
    lam <- eigen(V.tau,symmetric=TRUE)$values
    pval.tau <- davies(stat.tau, lam,lim=5000,acc=1e-04)$Qq
    if(pval.tau<=0) pval.tau  = liu(stat.tau,lam)
        
    ## overall p-value
    pval.fisher <-  fisher(c(pval.theta,pval.tau))
    pval.tippet <-  tippet(c(pval.theta,pval.tau))
    pval.dapc <-  DAPC(c(pval.theta,pval.tau))

    out <- c(p.fisher=pval.fisher,p.tippet=pval.tippet,p.dapc=pval.dapc)
    return(out)
}


     
