## DKAT
DKATr <- function(y,X,G){
    ## check the format of Y, if its a vector, convert it to matrix
    if(is.vector(y)|is.factor(y)) {
        y <- as.factor(y)
        J <- length(levels(y))
    } else if(is.matrix(y)) {
        if(any(rowSums(y)!=1)) stop("Error: summation of each row must equal to 1!")
        J <- ncol(y)
        n <- nrow(y)
    } else {
        stop("Error: y must be a vector or matrix!")
    }
    tmp = NULL
    for(j in 1:J) tmp = rbind(tmp,DKAT(y,X,G,rlevel=levels(y)[j]))
    out1 = min(min(tmp[,1],na.rm = FALSE)*J,1)
    out2 = min(min(tmp[,2],na.rm = FALSE)*J,1)
    return(c(out1,out2))
}

## uSKAT
uSKATr <- function(y,X,G){
    ## check the format of Y, if its a vector, convert it to matrix
    if(is.vector(y)|is.factor(y)) {
        y <- as.factor(y)
        J <- length(levels(y))
    } else if(is.matrix(y)) {
        if(any(rowSums(y)!=1)) stop("Error: summation of each row must equal to 1!")
        J <- ncol(y)
        n <- nrow(y)
    } else {
        stop("Error: y must be a vector or matrix!")
    }
    tmp = NULL
    for(j in 1:J) tmp = c(tmp,uSKAT(y,X,G,rlevel=levels(y)[j]))
    out = min(min(tmp,na.rm = FALSE)*J,1)
    return(out)
}

## uMiST
uMiSTr <- function(y,X,G,W){
        ## check the format of Y, if its a vector, convert it to matrix
    if(is.vector(y)|is.factor(y)) {
        y <- as.factor(y)
        J <- length(levels(y))
    } else if(is.matrix(y)) {
        if(any(rowSums(y)!=1)) stop("Error: summation of each row must equal to 1!")
        J <- ncol(y)
        n <- nrow(y)
    } else {
        stop("Error: y must be a vector or matrix!")
    }
    tmp = NULL
    for(j in 1:J) tmp = c(tmp,uMiST(y,X,G,W,rlevel=levels(y)[j]))
    out = min(min(tmp,na.rm = FALSE)*J,1)
    return(out)
}


## mRand
mRandr <- function(y,X,G){
       ## check the format of Y, if its a vector, convert it to matrix
    if(is.vector(y)|is.factor(y)) {
        y <- as.factor(y)
        J <- length(levels(y))
    } else if(is.matrix(y)) {
        if(any(rowSums(y)!=1)) stop("Error: summation of each row must equal to 1!")
        J <- ncol(y)
        n <- nrow(y)
    } else {
        stop("Error: y must be a vector or matrix!")
    }
    tmp = NULL
    for(j in 1:J) tmp = c(tmp,mRand(y,X,G,rlevel=levels(y)[j]))
    out = min(min(tmp,na.rm = FALSE)*J,1)
    return(out)
}

## SASOM
SASOMr <- function(y,X,G,W){
    ## check the format of Y, if its a vector, convert it to matrix
    if(is.vector(y)|is.factor(y)) {
        y <- as.factor(y)
        J <- length(levels(y))
    } else if(is.matrix(y)) {
        if(any(rowSums(y)!=1)) stop("Error: summation of each row must equal to 1!")
        J <- ncol(y)
        n <- nrow(y)
    } else {
        stop("Error: y must be a vector or matrix!")
    }
    tmp = NULL
    for(j in 1:J) tmp = rbind(tmp,SASOM(y,X,G,W,rlevel=levels(y)[j]))
    out = pmin(apply(tmp,2,min,na.rm=FALSE)*J,1)
    return(out)
}


