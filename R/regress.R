regress <- function(formula, Vformula, identity=TRUE, start=NULL, fraction=1, pos, verbose=0, gamVals=NULL, maxcyc=50, tol=1e-4, data, print.level=NULL){

  ## Vformula can just be something like ~ V0 + V1
  ## or leave it out or Vformula=NULL
  ## assume its in the form ~ V1 + V2 + ... + Vn or missing or Vformula=NULL
  ## for random effects and random interactions for factors A and B include
  ## ~ A + B + I(A:B)

  if(!is.null(print.level)) {
    cat("\nWarning: print.level has been replaced by verbose and has been deprecated.\nIt will be removed in the next version of regress\n\n")
    verbose <- print.level
  }

  if(missing(data)) data <- environment(formula)
  mf <- model.frame(formula,data=data,na.action=na.pass)
  mf <- eval(mf,parent.frame())
  y <- model.response(mf)

  X <- model.matrix(formula,data=data)

  model <- list()
  ##model$formula <- formula
  ##model$Vformula <- Vformula
  model <- c(model,mf)

  if(missing(Vformula)) Vformula <- NULL

  if(!is.null(Vformula))
  {
      V <- model.frame(Vformula,data=data,na.action=na.pass)
      V <- eval(V, parent.frame())
      Vcoef.names <- names(V)
      V <- as.list(V)
      k <- length(V)
  } else {
      V <- NULL
      k <- 0
      Vcoef.names=NULL
  }

  ## Remove missing values
  isNA <- is.na(y)
  ##y <- y[isNA==F]
  y <- na.omit(y)
  n <- length(y)
  Xcolnames <- dimnames(X)[[2]]
  if(is.null(Xcolnames)) {
      Xcolnames <- paste("X.column",c(1:dim(as.matrix(X))[2]),sep="")
  }

  X <- X[isNA==F,]
  X <- matrix(X, n, length(X)/n)
  qr <- qr(X)
  rankQ <- n-qr$rank
  if(qr$rank) {
      X <- matrix(X[, qr$pivot[1:qr$rank]],n,qr$rank)
      Xcolnames <- Xcolnames[qr$pivot[1:qr$rank]]
  } else {
      cat("\nERROR: X has rank 0\n\n")
  }

  if(missing(maxcyc)) maxcyc <- 50
  if(missing(tol)) tol <- 1e-4
  delta <- 1

  ## remove missing values
  for(i in 1:k)
  {
      if(is.matrix(V[[i]]))
      {
          V[[i]] <- V[[i]][isNA==F,]
          V[[i]] <- V[[i]][,isNA==F]
      }
      if(is.factor(V[[i]]))
      {
          V[[i]] <- V[[i]][isNA==F]
      }
  }

  In <- diag(rep(1,n),n,n)

  if(identity) {
      ##if(k) for(i in k:1) V[[i+1]] <- V[[i]]
      ##V[[1]] <- In
      V[[k+1]] <- as.factor(1:n)
      names(V)[k+1] <- "In"
      k <- k+1

      ##Vcoef.names <- c("Id",Vcoef.names)
      Vcoef.names <- c(Vcoef.names,"In")
      Vformula <- as.character(Vformula)
      Vformula[2] <- paste(Vformula[2],"+In")
      Vformula <- as.formula(Vformula)
  }

  model <- c(model,V)
  model$formula <- formula
  model$Vformula <- Vformula

  ## specify which parameters are positive and which are negative
  ## pos = c(1,1,0) means first two parameters are positive, third is either
  if(!missing(pos)) pos <- as.logical(pos)
  if(missing(pos)) pos <- rep(FALSE,k)
  pos <- c(pos,rep(FALSE,k))
  pos <- pos[1:k]

  ## Sherman Woodley identities for matrix inverses can be brought to bear here
  SWsolveINDICATOR <- FALSE
  if(all(sapply(V,is.factor))) {
      SWsolveINDICATOR <- TRUE
      Z <- list() ## see use below -
      for(i in 1:length(V))
      {
          if(is.factor(V[[i]])){
              Vi <- model.matrix(~V[[i]]-1)
              Z[[i]] <- Vi ## list of model matrices
              V[[i]] <- tcrossprod(Vi)
          }
      }
  } else {
      for(i in 1:length(V))
      {
          if(is.factor(V[[i]])){
              Vi <- model.matrix(~V[[i]]-1)
              Vi <- tcrossprod(Vi)
              V[[i]] <- Vi
          }
      }
  }
  ## So V is always a list of variance coavriance matrices, Z contains
  ## the model matrices of factors when we need to invoke the Sherman
  ## Woodley identities

  ## Expected Fisher Information
  A <- matrix(rep(0, k^2), k, k)
  entries <- expand.grid(1:k,1:k)

  x <- rep(0,k)
  sigma <- c(1,rep(0, k-1))
  stats <- rep(0, 0)

  ## START ALGORITHM

  if(!is.null(start)) {
      ## pad start with zeros if required
      start <- c(start, rep(1,k))
      start <- start[1:k]
  }

  if(k>2 && is.null(start)) start <- rep(var(y,na.rm=TRUE),k)
  if(k==1 && is.null(start)) start <- var(y,na.rm=TRUE)

  if(is.null(start) && k==2) {
      if(missing(gamVals)) {
          gamVals <- seq(0.01,0.02,length=3)^2
          gamVals <- sort(c(gamVals,seq(0.1,0.9,length=3),1-gamVals))
          gamVals <- 0.5
      }
      if(length(gamVals)>1) {
          if(verbose>=1) cat("Evaluating the llik at gamma = \n")
          if(verbose>=1) cat(gamVals)
          if(verbose>=1) cat("\n")
          reg.obj <- reml(gamVals,y,X,V[[1]],V[[2]],verbose=verbose)
          llik <- reg.obj$llik
          llik <- as.real(llik)
          if(verbose>=2) cat(llik,"\n")
          gam <- gamVals[llik==max(llik)]
          gam <- gam[1]
          if(verbose>=2) cat("MLE is near",gam,"and llik =",max(llik),"there\n")
      }
      if(length(gamVals)==1) {
          ## go straight to the Newton Raphson at gamVals
          gam <- gamVals[1]
          reg.obj <- list(rms=var(y))
      }
      start <- c(1-gam,gam)*reg.obj$rms
      ## it tends to take huge steps when starting at gam=0.9999
      if(gam==0.9999) {
          fraction <- fraction/100
          maxcyc <- maxcyc*10
      }
      if(verbose>=1) cat(c("start algorithm at",round(start,4),"\n"))
  }

  if(is.null(start) & k>2) {
      ## Never gets here by default - but this could be implemented,
      ## though it does add on a few extra iterations at the
      ## start.... not necessary in basic examples

      LLvals <- NULL
      ## equal weights
      V2 <- V[[2]]
      for(ii in 3:k) V2 <- V2 + V[[ii]]
      LLvals <- c(LLvals,reml(0.5,y,X,V[[1]],V2)$llik)
      ## Most at one end
      V2 <- V[[1]] + V2 ## total of all Vs
      for(ii in 1:k) {
          V2 <- V2 - V[[ii]]
          LLvals <- c(LLvals,reml(0.75,y,X,V2,V[[ii]])$llik)
      }
      best <- which.max(LLvals)
      if(verbose) {
          cat("Checking starting points\n")
          cat("llik values of", LLvals, "\n")
      }
      if(best==1) {
          start <- rep(var(y,na.rm=TRUE),k)
      } else {
          start <- rep(0.25,k)
          start[best] <- 0.75
      }
  }

  sigma <- start
  ## reparameterise so everything will get into the correct spot after exp
  ind <- which(pos)
  if(length(ind)) sigma[ind] <- log(start[ind])

  ## Set the memory requirements beforehand
  T <- vector("list", length=k)
  for(ii in 1:k) T[[ii]] <- matrix(NA,n,n)

  for(cycle in 1:maxcyc){

      ## newton-raphson step
      ## at the end of the iteration x will be -l'(sigma)/l''(sigma)
      ## hence we add instead of subtract
      ## fraction controls the proportion of each step to take
      ## for some reason 0.5 works very well

      if(all(pos==1)) x <- sign(x) * pmin(abs(x),5) ## limit maximum shift we can take in one step to 5 units on log scale
      sigma <- sigma+fraction*x

      ind <- which(pos)
      if(length(ind)) {
          sigma[ind] <- pmin(sigma[ind],20)
          sigma[ind] <- pmax(sigma[ind],-20) ## so on regular scale everything is between exp(-20) and exp(20)
      }

      if(verbose>=1) {
          cat(cycle, " ")
          sigmaTemp <- sigma
          ind <- which(pos)
          if(length(ind)) sigmaTemp[ind] <- exp(sigmaTemp[ind])
          cat(sigmaTemp)
      }

      ##Sigma <- matrix(0,dim(V[[1]])[1],dim(V[[1]])[2])
      if(!SWsolveINDICATOR) {
          Sigma <- 0
          ## can we get rid of this loop?
          for(i in 1:k)
          {
              if(pos[i]) {
                  Sigma <- Sigma + V[[i]]*exp(sigma[i])
              } else {
                  Sigma <- Sigma + V[[i]]*sigma[i]

              }
          }

          cholesky <- try(chol(Sigma, pivot=FALSE), silent=TRUE)
          if(class(cholesky) == "try-error" || min(diag(cholesky)) < 1e-9) return(1e+32)
          W <- chol2inv(cholesky)
          WX <- W %*% X
      } else {
          ind <- which(pos)
          sigmaTemp <- sigma
          if(length(ind)) sigmaTemp[ind] <- exp(sigmaTemp[ind])
          W <- SWsolve2(Z[1:(k-1)],sigmaTemp)
          WX <- W %*% X
      }

      ##WX <- solve(Sigma,cbind(X,In))
      ##W <- WX[,(dim(X)[2]+1):dim(WX)[2]]
      ##WX <- WX[,1:dim(X)[2]]
      XtWX <- t(X)%*%WX
      WQ <- W - WX%*%solve(XtWX,t(WX))
      rss <- as.numeric(t(y) %*% WQ %*% y)
      ldet <- sum(log(eigen(WQ,symmetric=TRUE,only.values=TRUE)$values[1:rankQ]))

      ## REML LLIK
      rllik1 <- ldet/2 - rss/2

      ## MARGINAL LLIK BASED ON RESIDUAL CONFIGURATION STATISTIC
      rllik2 <- ldet/2 - rankQ * log(rss)/2

      if(verbose) cat(" REML =",rllik1,"\n")

      ## change in REML llik
      if(cycle > 1) delta <- rllik1 - rllik0
      rllik0 <- rllik1

      ## now the fun starts, derivative and expected fisher info
      ## the 0.5 multiple is ignored, it is in both and they cancel

      ##T <- list(NULL)
      x <- NULL

      ## derivatives are now D[[i]] = var.components[i]*V[[i]]
      var.components <- rep(1,k)
      ind <- which(pos)
      if(length(ind)) var.components[ind] <- exp(sigma[ind])

      ## Slow part - order k n-squared
      if(!SWsolveINDICATOR) {
          if(identity) {
              T[[k]] <- WQ
              if(k>1) {
                  for(ii in (k-1):1) T[[ii]] <- WQ %*% V[[ii]]
              }
          } else {
              for(ii in 1:k) T[[ii]] <- WQ %*% V[[ii]]
          }
      } else {
          if(identity) {
              T[[k]] <- WQ
              if(k>1) {
                  for(ii in (k-1):1) T[[ii]] <- tcrossprod(WQ %*% Z[[ii]],Z[[ii]])
              }
          } else {
              for(ii in 1:k) T[[ii]] <- tcrossprod(WQ %*% Z[[ii]],Z[[ii]])
          }

          ##if(k>=6) {
          ##    ## One line to do all - may be memory inefficient though
          ##    T <- lapply(Z,function(x) tcrossprod(WQ %*% x, x))
          ##} else {
          ##    if(k>=1) T[[1]] <- tcrossprod(WQ %*% Z[[1]], Z[[1]])
          ##    if(k>=2) T[[2]] <- tcrossprod(WQ %*% Z[[2]], Z[[2]])
          ##    if(k>=3) T[[3]] <- tcrossprod(WQ %*% Z[[3]], Z[[3]])
          ##    if(k>=4) T[[4]] <- tcrossprod(WQ %*% Z[[4]], Z[[4]])
          ##    if(k>=5) T[[5]] <- tcrossprod(WQ %*% Z[[5]], Z[[5]])
          ##}
      }

      x <- sapply(T,function(x) as.numeric(t(y) %*% x %*% WQ %*% y - sum(diag(x))))
      x <- x * var.components


      ## See nested for loops commented out below - evaluating the Expected Fisher Information, A
      ff <- function(x) sum(T[[x[1]]] * t(T[[x[2]]])) * var.components[x[1]] * var.components[x[2]]
      aa <- apply(entries,1,ff)
      A[as.matrix(entries)] <- aa

      ##for(i in 1:k)
      ##{
      ##if(identity && i==1 && pos[i]) {
      ##  T[[i]] <- WQ
      ##} else

      ##if(SWsolveINDICATOR==FALSE) {
      ##    T[[i]] <- WQ %*% V[[i]]
      ##} else {
      ##    T[[i]] <- WQ %*% tcrossprod(V[[i]])
      ##}
      ##x[i] <- as.numeric(t(y) %*% T[[i]] %*% WQ %*% y - sum(diag(T[[i]])))
      ##x[i] <- x[i]*var.components[i]

      ##for(j in c(1:i))
      ##  {
      ## expected fisher information
      ## ##A[j,i] <- Tr(T[[j]] %*% T[[i]])
      ## ##the Ts are not symmetric, hence the transpose below

      ##A[j,i] <- sum(T[[j]] * t(T[[i]]))
      ##A[j,i] <- A[j,i]*var.components[i]*var.components[j]
      ##A[i,j] <- A[j,i]
      ##}
      ##}

      stats <- c(stats, rllik1, rllik2, sigma[1:k], x[1:k])
      if(verbose==-1) {
          ##cat(c(rllik1, rllik2, sigma[1:k], x[1:k]),"\n")
      }

      A.svd <- ginv(A)
      x <- A.svd %*% x

      if(qr(A)$rank < k){
          if(cycle==1) {
              if(verbose) {
                  cat("Warning: Non identifiable dispersion model\n")
                  ##print(round(A,6))
                  sigmaTemp <- sigma
                  ind <- which(pos)
                  if(length(ind)) sigmaTemp[ind] <- exp(sigmaTemp[ind])
                  cat(sigmaTemp)
                  cat("\n")
              }
          }
      }

      ## check the change in llik is small
      if(abs(delta) < tol) break
  }

  if(cycle==maxcyc)
  {
      ## issue a warning
      if(verbose) cat("WARNING:  maximum number of cycles reached before convergence\n")
  }
  y

  llik <- as.numeric(stats)
  llik <- t(matrix(llik, 2*k+2, length(llik)/(2*k+2)))

  cov <- XtWX
  cov <- solve(cov, cbind(t(WX),diag(1,dim(XtWX)[1])))
  beta.cov <- matrix(cov[,(dim(t(WX))[2]+1):dim(cov)[2]],dim(X)[2],dim(X)[2])
  cov <- matrix(cov[,1:dim(t(WX))[2]],dim(X)[2],dim(X)[1])

  beta <- cov %*% y
  beta <- matrix(beta,length(beta),1)
  row.names(beta) <- Xcolnames
  beta.se <- sqrt(abs(diag(beta.cov)))
  pos.cov <- (diag(beta.cov)<0)
  beta.se[pos.cov] <- NA
  beta.se <- matrix(beta.se,length(beta.se),1)
  row.names(beta.se) <- Xcolnames
  rms <- rss/rankQ

  fitted.values <- X%*%beta
  Q <- In -  X %*% cov

  predicted <- NULL
  if(identity) {
      gam <- sigma[k]  ## coefficient of identity, last variance term
      if(pos[k]) gam <- exp(gam)
      if(SWsolveINDICATOR) {
          ## Sigma is undefined
          Sigma <- 0
          for(i in 1:k)
          {
              if(pos[i]) {
                  Sigma <- Sigma + V[[i]]*exp(sigma[i])
              } else {
                  Sigma <- Sigma + V[[i]]*sigma[i]
              }
          }
      }
      predicted <- fitted.values + (Sigma - gam*In) %*% W%*%(y - fitted.values)
  }

  ## scale dictated by pos
  sigma.cov <- (A.svd[1:k, 1:k] * 2)

  ## Last Step:  June 17th 2005
  ## Convert the estimates for the variance parameters, their standard
  ## errors etc to the usual scale

  ind <- which(pos)
  if(length(ind)) sigma[ind] <- exp(sigma[ind])

  FI <- A/2

  ## convert FI using pos
  FI.c <- matrix(0,dim(FI)[1],dim(FI)[2])
  FI.c <- FI / tcrossprod((sigma-1)*pos+1)
  ##for(i in 1:dim(FI)[1])
  ##  for(j in 1:dim(FI)[2])
  ##    FI.c[i,j] <- FI[i,j]/(((sigma[i]-1)*pos[i]+1)*((sigma[j]-1)*pos[j]+1))

  sigma.cov <- ginv(FI.c)


  result <- list(trace=llik, llik=llik[cycle, 1], cycle=cycle,
                 rdf=rankQ, beta=beta, beta.cov=beta.cov, beta.se=beta.se,
                 sigma=sigma[1:k], sigma.cov=sigma.cov[1:k,1:k], W=W, Q=Q,
                 fitted=fitted.values, predicted=predicted, pos=pos,
                 Vnames=Vcoef.names, formula=formula, Vformula=Vformula, model=model)
  class(result) <- "regress"
  result
}

## generalised inverse

ginv <- function (X, tol = sqrt(.Machine$double.eps))
{
  ## taken from library MASS
  if (length(dim(X)) > 2 || !(is.numeric(X) || is.complex(X)))
    stop("X must be a numeric or complex matrix")
  if (!is.matrix(X))
    X <- as.matrix(X)
  Xsvd <- svd(X)
  if (is.complex(X))
    Xsvd$u <- Conj(Xsvd$u)
  Positive <- Xsvd$d > max(tol * Xsvd$d[1], 0)
  if (all(Positive))
    Xsvd$v %*% (1/Xsvd$d * t(Xsvd$u))
  else if (!any(Positive))
    array(0, dim(X)[2:1])
  else Xsvd$v[, Positive, drop = FALSE] %*% ((1/Xsvd$d[Positive]) *
                            t(Xsvd$u[, Positive, drop = FALSE]))
}

summary.regress <- function(object, ...) object

print.regress <- function(x, digits=3, fixed.effects=T, ...)
  {
    cat("\nMaximised Residual Log Likelihood is",round(x$llik,digits),"\n",sep=" ")
    indent.lin <- max(nchar(dimnames(x$beta)[[1]]))
    indent.var <- max(nchar(x$Vnames))
    indent <- max(indent.lin,indent.var)

    extra.space <- ""
    space.var <- extra.space
    for(i in 0:(indent-indent.var)) space.var <- paste(space.var," ",sep="")
    space.lin <- extra.space
    for(i in 0:(indent-indent.lin)) space.lin <- paste(space.lin," ",sep="")

    coefficients <- cbind(x$beta,x$beta.se)
    dimnames(coefficients)[[2]] <- c("Estimate","Std. Error")
    coefficients <- round(coefficients,digits)
    if(fixed.effects) {
      cat("\nLinear Coefficients:\n")
      row.names(coefficients) <- paste(space.lin,dimnames(x$beta)[[1]],sep="")
      print(coefficients)
      cat("\n")
    } else {
      cat("\nLinear Coefficients: not shown\n\n")
    }

    ## New version of regress automatically converts to the linear
    ## scale - as if pos was a vector of zeroes

    var.coefficients <- cbind(x$sigma,sqrt(diag(as.matrix(x$sigma.cov))))
    row.names(var.coefficients) <- paste(space.var,x$Vnames,sep="")
    dimnames(var.coefficients)[[2]] <- c("Estimate","Std. Error")
    var.coefficients <- round(var.coefficients,digits)
    cat("Variance Coefficients:\n")
    print(var.coefficients)
    cat("\n")
  }

## when two matrices are passed to regress this is also called
## to evaluate the REML at certain values of gamma and find a
## good place to start the regress algorithm

reml <- function(lambda, y, X, V0, V1,verbose=0){

  if(is.null(dim(y)))
    {
      isNA <- is.na(y)
      y <- y[isNA==F]
    } else {
      isNA <- apply(y,1,is.na)

      if(is.matrix(isNA))  {
        isNA <- as.logical(apply(isNA,2,sum))
      }
      y <- y[isNA==F,]
    }
  V0 <- V0[isNA==F,isNA==F]
  V1 <- V1[isNA==F,isNA==F]
  X <- X[isNA==F,]
  X <- as.matrix(X)

  qr <- qr(X)
  ##print(qr$rank)
  n <- dim(as.matrix(y))[1]
  In <- diag(1,n)

  X <- matrix(X[, qr$pivot[1:qr$rank]],n,qr$rank)
  llik <- rep(0, length(lambda))
  if(is.null(dim(y))) q <- 1 else q <- dim(y)[2]

  n <- dim(X)[1]
  if(missing(V0)) V0 <- diag(rep(1, n), n, n)
  rank <- n - qr$rank
  ##if(verbose==1) cat("n-p =",n,"-",qr$rank,"=",rank,"\n")
  for(i in 1:length(lambda))
    {
      if(verbose>=2) cat(lambda[i],"\n")

      Sigma <- (1-lambda[i])*V0 + lambda[i] * V1
       cholesky <- try(chol(Sigma, pivot=FALSE), silent=TRUE)
      if(class(cholesky) == "try-error" || min(diag(cholesky)) < 1e-9) return(1e+32)
      W <- chol2inv(cholesky)
      WX <- W %*% X
      ##WX <- solve(Sigma,cbind(X,In))
      ##W <- WX[,(dim(X)[2]+1):dim(WX)[2]]
      ##WX <- WX[,1:dim(X)[2]]
      XtWX <- t(X)%*%WX
      WQ <- W - WX%*%solve(XtWX,t(WX))
      rss <- t(y) %*% WQ %*% y
      logdetrss <- sum(log(eigen(rss)$values[1:q]))
      eVals <- eigen(WQ,symmetric=TRUE,only.values=TRUE)$values[1:rank]
      ldet <- sum(log(eVals))
      llik[i] <- Re(ldet*q/2 - rank*logdetrss/2)
    }
  imax <- sort.list(-llik)[1]
  lambdamax <- lambda[imax]
  curv <- 0
  if(imax > 1 && imax < length(lambda)){
    delta <- (lambda[imax+1] - lambda[imax-1])/2
    slope <-  (llik[imax+1] - llik[imax-1])/2
    curv <- llik[imax-1] -2*llik[imax] + llik[imax+1]
    lambdamax <- lambdamax - slope/curv * delta
    curv <- -curv/delta^2
  }
  lamMax <- lambdamax
  Sigma <- (1-lamMax)*V0 + lamMax * V1
  cholesky <- try(chol(Sigma, pivot=FALSE), silent=TRUE)
  if(class(cholesky) == "try-error" || min(diag(cholesky)) < 1e-9) return(1e+32)
  W <- chol2inv(cholesky)
  WX <- W %*% X
  ##WX <- solve(Sigma,cbind(X,In))
  ##W <- WX[,(dim(X)[2]+1):dim(WX)[2]]
  ##WX <- WX[,1:dim(X)[2]]
  XtWX <- t(X)%*%WX
  FItWX <- solve(XtWX,t(WX))
  WQ <- W - WX%*%FItWX
  rss <- t(y) %*% WQ %*% y
  beta <- FItWX %*% y

  list(llik=as.numeric(llik),rms=rss/rank, beta=beta, gamma=lambda, gamMax=lambdamax,W=W)
}

remlOptimize <- function(y, X, V0, V1,verbose=0,...){

  if(is.null(dim(y)))
    {
      isNA <- is.na(y)
      y <- y[isNA==F]
    } else {
      isNA <- apply(y,1,is.na)

      if(is.matrix(isNA))  {
        isNA <- as.logical(apply(isNA,2,sum))
      }
      y <- y[isNA==F,]
    }
  V0 <- V0[isNA==F,isNA==F]
  V1 <- V1[isNA==F,isNA==F]
  X <- X[isNA==F,]
  X <- as.matrix(X)

  qr <- qr(X)
  ##print(qr$rank)
  n <- dim(as.matrix(y))[1]
  In <- diag(1,n)

  X <- matrix(X[, qr$pivot[1:qr$rank]],n,qr$rank)
  if(is.null(dim(y))) q <- 1 else q <- dim(y)[2]

  n <- dim(X)[1]
  if(missing(V0)) V0 <- diag(rep(1, n), n, n)
  rank <- n - qr$rank
  ##if(verbose==1) cat("n-p =",n,"-",qr$rank,"=",rank,"\n")

  f <- function(lambda,verbose=verbose) {
    if(verbose>=2) cat(lambda,"\n")
    Sigma <- (1-lambda)*V0 + lambda * V1
    cholesky <- try(chol(Sigma, pivot=FALSE), silent=TRUE)
    if(class(cholesky) == "try-error" || min(diag(cholesky)) < 1e-9) return(1e+32)
    W <- chol2inv(cholesky)
    WX <- W %*% X
    ##WX <- solve(Sigma,cbind(X,In))
    ##W <- WX[,(dim(X)[2]+1):dim(WX)[2]]
    ##WX <- WX[,1:dim(X)[2]]
    XtWX <- t(X)%*%WX
    WQ <- W - WX%*%solve(XtWX,t(WX))
    rss <- t(y) %*% WQ %*% y
    logdetrss <- sum(log(eigen(rss)$values[1:q]))
    eVals <- eigen(WQ,symmetric=TRUE,only.values=TRUE)$values[1:rank]
    ldet <- sum(log(eVals))
    llik <- Re(ldet*q/2 - rank*logdetrss/2)
    llik
  }

  res <- optimize(f,interval=c(0,1),maximum=TRUE,verbose=verbose,...)
  lamMax <- res$maximum
  llikMax <- res$objective

  Sigma <- (1-lamMax)*V0 + lamMax * V1
  cholesky <- try(chol(Sigma, pivot=FALSE), silent=TRUE)
  if(class(cholesky) == "try-error" || min(diag(cholesky)) < 1e-9) return(1e+32)
  W <- chol2inv(cholesky)
  WX <- W %*% X
  ##WX <- solve(Sigma,cbind(X,In))
  ##W <- WX[,(dim(X)[2]+1):dim(WX)[2]]
  ##WX <- WX[,1:dim(X)[2]]
  XtWX <- t(X)%*%WX
  FItWX <- solve(XtWX,t(WX))
  WQ <- W - WX%*%FItWX
  rss <- t(y) %*% WQ %*% y
  beta <- FItWX %*% y

  list(llik=as.numeric(llikMax),rms=rss/rank, beta=beta, gamMax=lamMax,W=W)
}

SWsolve <- function(S,K,D,Dinv=NULL,b) {
    ## solve(a,b) where a has the form SKS' + D using the Sherman Woodley identities

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
