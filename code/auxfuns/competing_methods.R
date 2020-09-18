####################################################################################
#### DKAT code is downloaded from https://github.com/xyz5074/DKAT/tree/master/R ###
####################################################################################


DKAT <- function(y,X,G,rlevel){
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

    fit = multinom(ymat~0+X,trace=FALSE)
    mu = fit$fitted.values
    rc <- which(colnames(ymat) == rlevel)
    re = (ymat-mu)[,-rc]
    KX = wlin.kernel(G,c(1,1)) # flat weights 
    KY1 = pheno.kernel(re,0.1,TRUE)
    KY2 = pheno.kernel(re,0.1,FALSE)
    pv = c(d_DKAT(KX,KY2),d_DKAT(KX,KY1))
    names(pv) = c('DKAT-I','DKAK-T')
    return(pv)
}

wlin.kernel <- function(X, W.beta) {
  MAF=colMeans(X)/2  ## colMeans is in package MSKAT
  a1=W.beta[1]
  a2=W.beta[2]
  w=dbeta(MAF,a1,a2)
  W=diag(w)
  return( X%*%W%*%W%*%t(X) )  
}


pheno.kernel <- function(Y, rho=0.1,glasso = TRUE){
  if(glasso){
    s=var(Y)
    N=glasso(s,rho=0.1)$wi
    Ky=Y%*%N%*%t(Y)
  } else {
    Ky=Y%*%t(Y)
  }
  return(Ky)  
}

d_DKAT <- function(K,L){
  tr=function(x){return(sum(diag(x)))}
  n=nrow(K)
  I.n=diag(1,n)
  I.1=rep(1,n)
  H=I.n-I.1%*%t(I.1)/n
  K=H%*%K%*%H
  L=H%*%L%*%H
  A=K/tr(K%*%K)  
  W=L/tr(L%*%L)
  Fstar=tr(A%*%W) ## 
  
  mean.krv=tr(A)*tr(W)/(n-1)## mean of DKAT 
  T=tr(A);T2=tr(A%*%A);S2=sum(diag(A)^2)
  Ts=tr(W);T2s=tr(W%*%W);S2s=sum(diag(W)^2)
  temp1=2*((n-1)*T2-T^2)*((n-1)*T2s-Ts^2)/(n-1)^2/(n+1)/(n-2)
  temp21=n*(n+1)*S2- (n-1)*(T^2+2*T2)
  temp22=n*(n+1)*S2s- (n-1)*(Ts^2+2*T2s)
  temp23=(n+1)*n*(n-1)*(n-2)*(n-3)
  temp2=temp21*temp22/temp23
  variance.krv=temp1+temp2## variance of DKAT
  
  T3=tr(A%*%A%*%A);S3=sum(diag(A)^3);U=sum(A^3);R=t(diag(A))%*%diag(A%*%A);B=t(diag(A))%*%A%*%diag(A)
  T3s=tr(W%*%W%*%W);S3s=sum(diag(W)^3);Us=sum(W^3);Rs=t(diag(W))%*%diag(W%*%W);Bs=t(diag(W))%*%W%*%diag(W)
  t1=n^2*(n+1)*(n^2+15*n-4)*S3*S3s
  t2=4*(n^4-8*n^3+19*n^2-4*n-16)*U*Us
  t3=24*(n^2-n-4)*(U*Bs+B*Us)
  t4=6*(n^4-8*n^3+21*n^2-6*n-24)*B*Bs
  t5=12*(n^4-n^3-8*n^2+36*n-48)*R*Rs
  t6=12*(n^3-2*n^2+9*n-12)*(T*S2*Rs+R*Ts*S2s)
  t7=3*(n^4-4*n^3-2*n^2+9*n-12)*T*Ts*S2*S2s
  t81=(n^3-3*n^2-2*n+8)*(R*Us+U*Rs);t82=(n^3-2*n^2-3*n+12)*(R*Bs+B*Rs)
  t8=24*(t81+t82)
  t9=12*(n^2-n+4)*(T*S2*Us+U*Ts*S2s)
  t10=6*(2*n^3-7*n^2-3*n+12)*(T*S2*Bs+B*Ts*S2s)
  t11=-2*n*(n-1)*(n^2-n+4)*((2*U+3*B)*S3s+(2*Us+3*Bs)*S3)
  t12=-3*n*(n-1)^2*(n+4)*((T*S2+4*R)*S3s+(Ts*S2s+4*Rs)*S3)
  t13=2*n*(n-1)*(n-2)*((T^3+6*T*T2+8*T3)*S3s+(Ts^3+6*Ts*T2s+8*T3s)*S3)
  t14=T^3*((n^3-9*n^2+23*n-14)*Ts^3+6*(n-4)*Ts*T2s+8*T3s)
  t15=6*T*T2*((n-4)*Ts^3+(n^3-9*n^2+24*n-14)*Ts*T2s+4*(n-3)*T3s)
  t16=8*T3*(Ts^3+3*(n-3)*Ts*T2s+(n^3-9*n^2+26*n-22)*T3s)
  t17=-16*(T^3*Us+U*Ts^3)-6*(T*T2*Us+U*Ts*T2s)*(2*n^2-10*n+16)
  t18=-8*(T3*Us+U*T3s)*(3*n^2-15*n+16)-(T^3*Bs+B*Ts^3)*(6*n^2-30*n+24)
  t19=-6*(T*T2*Bs+B*Ts*T2s)*(4*n^2-20*n+24)-8*(T3*Bs+B*T3s)*(3*n^2-15*n+24)
  t201=24*(T^3*Rs+R*Ts^3)+6*(T*T2*Rs+R*Ts*T2s)*(2*n^2-10*n+24)
  t202=8*(T3*Rs+R*T3s)*(3*n^2-15*n+24)+(3*n^2-15*n+6)*(T^3*Ts*S2s+T*S2*Ts^3)
  t203=6*(T*T2*Ts*S2s+Ts*T2s*T*S2)*(n^2-5*n+6)+48*(T3*Ts*S2s+T3s*T*S2)
  t20=-(n-2)*(t201+t202+t203)
  temp31=t1+t2+t3+t4+t5+t6+t7+t8+t9+t10+t11+t12+t13+t14+t15+t16+t17+t18+t19+t20
  temp32=n*(n-1)*(n-2)*(n-3)*(n-4)*(n-5)
  mom3=temp31/temp32
  skewness.krv= (mom3-3*mean.krv*variance.krv-mean.krv^3)/variance.krv^1.5 ## skewness of DKAT
  
  m1=mean.krv
  m2=variance.krv
  m3=skewness.krv
  shape=4/m3^2
  scale=sqrt(m2)*m3/2
  location=m1-2*sqrt(m2)/m3
  PIIIpars=list(shape,location,scale)
  pv=1-ppearsonIII(Fstar, params=PIIIpars) 
  if (is.na(pv)) {pv=1}  ## usually happens if var=0 in the denominator
  return(pv)
}

####################################################################################
################################ burden method #####################################
####################################################################################
burden <- function(y,X,G){
    fit1 <- multinom(y~0+X,trace=FALSE)
    fit2 <- multinom(y~0+X+rowSums(G),trace=FALSE)
    ares = anova(fit1,fit2)
    return(ares$`Pr(Chi)`[2])
}

####################################################################################
################################ uSKAT method ######################################
####################################################################################
uSKAT <- function(y,X,G,rlevel){
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
    id1 = which(ymat[,1]==1)
    id2 = which(ymat[,2]==1)
    id3 = which(ymat[,3]==1)
    X1 = X[id1,,drop=FALSE]
    X2 = X[id2,,drop=FALSE]
    X3 = X[id3,,drop=FALSE]
    ymat1 = ymat[id1,]
    ymat2 = ymat[id2,]
    ymat3 = ymat[id3,]
    G1 = G[id1,,drop=FALSE]
    G2 = G[id2,,drop=FALSE]
    G3 = G[id3,,drop=FALSE]
    y12 = rbind(ymat1,ymat2)
    y12 = y12[,apply(y12,2,sum)>0]
    y13 = rbind(ymat1,ymat3)
    y13 = y13[,apply(y13,2,sum)>0]
    suppressWarnings({
        obj12 <- SKAT_Null_Model(y12[,1]~0+rbind(X1,X2),out_type='D',Adjustment = FALSE)
        obj13 <- SKAT_Null_Model(y13[,1]~0+rbind(X1,X3),out_type='D',Adjustment = FALSE)
        p12 <- SKAT(rbind(G1,G2), obj12, kernel="linear")$p.value 
        p13 <- SKAT(rbind(G1,G3), obj13, kernel="linear")$p.value
    })
    return(min(min(p12,p13)*2,1))
}

####################################################################################
################################ uMiST method ######################################
####################################################################################
ref_mist = function(y1,y2,x1,x2,g1,g2,w){
  y = rbind(y1,y2)
  y = y[,apply(y,2,sum)>0]  
  x = rbind(x1,x2)
  g = rbind(g1,g2)
  suppressWarnings({out = try(unlist(logit.test(y[,1],x,g,w,"liu"))[5],silent = TRUE)})
  if(!is.numeric(out)) out = NA
  return(out)
}

uMiST <- function(y,X,G,W,rlevel){ # mist-univariate
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
    W = t(W)
    ## adjust the order of columns based on reference level
    ## put the column correspond to reference level as the first column
    rc <- which(colnames(ymat) == rlevel)
    ymat <- ymat[,c(rc,setdiff(1:J,rc))]
    id1 = which(ymat[,1]==1)
    id2 = which(ymat[,2]==1)
    id3 = which(ymat[,3]==1)
    X1 = X[id1,,drop=FALSE]
    X2 = X[id2,,drop=FALSE]
    X3 = X[id3,,drop=FALSE]
    ymat1 = ymat[id1,]
    ymat2 = ymat[id2,]
    ymat3 = ymat[id3,]
    G1 = G[id1,,drop=FALSE]
    G2 = G[id2,,drop=FALSE]
    G3 = G[id3,,drop=FALSE]
    tmp = c(ref_mist(ymat1,ymat2,X1,X2,G1,G2,W),ref_mist(ymat1,ymat3,X1,X3,G1,G3,W))
    return(min(min(tmp)*2,1))
}

####################################################################################
################################ mRand  method #####################################
####################################################################################
## input
## y - ymat of matrix n*J or yvec of n
## rlevel - reference level, default is the first level, must be a number if ymat is given or factor when yvec is given
## X - covariant matrix n*m
## G - variant matrix n*p
## W - variant character matrix q*p

varfun <- function(mu){
    ## construct variance matrix
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

mRand <- function(y,X,G,rlevel=1){

    ## check the format of Y, if its a vector change it to matrix
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

    ## define new notations for futurn use
    Y.s <- matrix(as.vector(ymat[,-1]),ncol=1)
    X.s <- kronecker(diag(J-1),X)
    G.s <- kronecker(diag(J-1),G)

    ## fit a model under H0
    fit1 <- multinom(ymat~0+X,trace=FALSE)
    mu1 <- fit1$fitted.values[,-1]
    mu1.s <- as.vector(mu1)
    ## construct variance matrix under H0
    D1 <- varfun(mu1)
    ## score test for theta 
    A = t(G.s)%*%(Y.s-mu1.s)
    L = t(A)%*%A
    DX = D1%*%X.s
    
    V.tau = t(G.s)%*%(D1 - DX%*%solve(t(X.s)%*%DX)%*%t(DX))%*%(G.s)
    lam <- eigen(V.tau,symmetric=TRUE)$values
    pval.tau <- davies(L, lam)$Qq
    if(pval.tau<0) pval.tau  = liu(L,lam)
        
    return(pval.tau)
}



