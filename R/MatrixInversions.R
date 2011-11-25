SWsolve <- function(S,K,D,Dinv=NULL,b) {
    ## solve(a,b) where a has the form SKS' + D using the Sherman Morrison Woodbury identities

    if(is.matrix(K) & is.matrix(D) & !is.null(Dinv)) {
        ## Case 1 - all are matrices
        tSDi <- crossprod(S,Dinv)
        Kinv <- solve(K)
        ret <- solve(Kinv + tSDi %*% S, tSDi)
        ret <- Dinv - crossprod(tSDi,ret)
        if(!missing(b)) ret <- ret %*% b
        return(ret)
    }

    if(is.numeric(K) & !is.null(Dinv)) {
        tSDi <- crossprod(S,Dinv)
        ret <- solve(1/K * diag(ncol(S)) + tSDi %*% S, tSDi)
        ret <- Dinv - crossprod(tSDi,ret)
        if(!missing(b)) ret <- ret %*% b
        return(ret)
    }

    if(is.numeric(D) & is.matrix(K)) {
        ret <- 1/D * diag(nrow(S)) - 1/D^2 * S %*% solve(solve(K) + 1/D * crossprod(S),t(S))
        if(!missing(b)) ret <- ret %*% b
        return(ret)
    }

    if(is.numeric(K) & is.numeric(D)) {
        ret <- 1/D * diag(nrow(S)) - 1/D^2 * S %*% solve(1/K * diag(ncol(S)) + 1/D * crossprod(S),t(S))
        if(!missing(b)) ret <- ret %*% b
        return(ret)
    }
}


SWsolve2 <- function(Zlist,clist,b) {
    ## Invert a matrix of the form sum [ clist[i] tcrossprod(Zlist[[ii]])) ] using Sherman Morrison Woodbury Identities
    if(length(Zlist)!=(length(clist)-1)) stop()
    k <- length(Zlist)
    D <- clist[1] * tcrossprod(Zlist[[1]])
    diag(D) <- diag(D) + clist[k+1]
    Dinv <- SWsolve(Zlist[[1]],clist[1],clist[k+1])
    if(k==1) {
        if(!missing(b)) Dinv <- Dinv %*% b
        return(Dinv)
    }
    for(ii in 2:k) {
        Dinv <- SWsolve(Zlist[[ii]],clist[ii],D,Dinv)
        D <- D + clist[ii]*tcrossprod(Zlist[[ii]])
    }
    if(!missing(b)) Dinv <- Dinv %*% b
    return(Dinv)
}
