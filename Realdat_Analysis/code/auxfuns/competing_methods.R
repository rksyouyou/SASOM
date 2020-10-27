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
uSKAT <- function(y,X,G){
    
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

    X12 = X[ymat[,3]==0,]
    X13 = X[ymat[,2]==0,]
    X23 = X[ymat[,1]==0,]
    G12 = G[ymat[,3]==0,]
    G13 = G[ymat[,2]==0,]
    G23 = G[ymat[,1]==0,]
    y12 = ymat[ymat[,3]==0,1]
    y13 = ymat[ymat[,2]==0,1]
    y23 = ymat[ymat[,1]==0,2]

    out = rep(NA,3)
    suppressWarnings({
        obj12 <- SKAT_Null_Model(y12~0+X12,out_type='D',Adjustment = FALSE)
        obj13 <- SKAT_Null_Model(y13~0+X13,out_type='D',Adjustment = FALSE)
        obj23 <- SKAT_Null_Model(y23~0+X23,out_type='D',Adjustment = FALSE)
        p12 <- SKAT(G12, obj12, kernel="linear")$p.value 
        p13 <- SKAT(G13, obj13, kernel="linear")$p.value
        p23 <- SKAT(G23, obj23, kernel="linear")$p.value
        out = c(p12,p13,p23)
    })
    return(min(min(out,na.rm=TRUE)*3,1))
}


####################################################################################
################################ uMiST method ######################################
####################################################################################

uMiST <- function(y,X,G,W){ # mist-univariate
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

    X12 = X[ymat[,3]==0,]
    X13 = X[ymat[,2]==0,]
    X23 = X[ymat[,1]==0,]
    G12 = G[ymat[,3]==0,]
    G13 = G[ymat[,2]==0,]
    G23 = G[ymat[,1]==0,]
    y12 = ymat[ymat[,3]==0,1]
    y13 = ymat[ymat[,2]==0,1]
    y23 = ymat[ymat[,1]==0,2]

    out = rep(NA,3)
    suppressWarnings({
        p12 = as.numeric(try(unlist(logit.test(y12,X12,G12,W,"liu"))[5],silent = TRUE))
        p13 = as.numeric(try(unlist(logit.test(y13,X13,G13,W,"liu"))[5],silent = TRUE))
        p23 = as.numeric(try(unlist(logit.test(y23,X23,G23,W,"liu"))[5],silent = TRUE))
        out = c(p12,p13,p23)
    })
    return(min(min(out,na.rm=TRUE)*3,1))
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

mRand <- function(y,X,G){

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
    p = ncol(G)

    ## fit a model under H0
    fit1 <- multinom(ymat~0+X,trace=FALSE)
    mu1 <- fit1$fitted.values[,-1]
    mu1.s <- as.vector(mu1)
    ## construct variance matrix under H0
    D1 <- varfun(mu1)
    ## score test for theta
    A = rbind(kronecker(diag(J-1),diag(p)),kronecker(matrix(1,ncol=J-1,nrow=1),diag(p)))
    B = A%*%t(G.s)%*%(Y.s-mu1.s)
    L = t(B)%*%B
    DX = D1%*%X.s
    V.tau = A%*%t(G.s)%*%(D1 - DX%*%solve(t(X.s)%*%DX)%*%t(DX))%*%(G.s)%*%t(A)
    lam <- eigen(V.tau,symmetric=TRUE)$values
    pval.tau <- davies(L, lam)$Qq
    if(pval.tau<0) pval.tau  = liu(L,lam)
        
    return(pval.tau)
}



