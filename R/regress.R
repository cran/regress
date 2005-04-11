regress <- function(formula, Vformula, identity=TRUE, start=NULL, fraction=1, pos, print.level=0, gamVals=NULL, maxcyc=50, tol=1e-4, data){

  ## Vformula can just be something like ~ V0 + V1
  ## or leave it out or Vformula=NULL
  ## assume its in the form ~ V1 + V2 + ... + Vn or missing or Vformula=NULL
  ## for random effects and random interactions for factors A and B include
  ## ~ A + B + I(A:B)

  if(missing(data)) data <- environment(formula)
  mf <- model.frame(formula,data=data,na.action=na.pass)
  mf <- eval(mf,parent.frame())
  y <- model.response(mf)

  X <- model.matrix(formula,data=data)
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

  In <- diag(rep(1,length(y)),length(y) ,length(y))

  if(identity) {
    ##if(k) for(i in k:1) V[[i+1]] <- V[[i]]
    ##V[[1]] <- In
    V[[k+1]] <- In
    k <- k+1
    ##Vcoef.names <- c("Id",Vcoef.names)
    Vcoef.names <- c(Vcoef.names,"I")
  }

  ## specify which parameters are positive and which are negative
  ## pos = c(1,1,0) means first two parameters are positive, third is either
  if(missing(pos)) pos <- rep(0,k)
  pos <- c(pos,rep(0,k))
  pos <- pos[1:k]
  
  for(i in 1:length(V))
    {
      if(is.factor(V[[i]])){
        Vi <- model.matrix(~V[[i]]-1)
        Vi <- Vi %*% t(Vi)
        V[[i]] <- Vi
      }
    } 

  A <- matrix(rep(0, k^2), k, k)

  x <- rep(0,k)
  sigma <- c(1,rep(0, k-1))
  stats <- rep(0, 0)

  In <- diag(1,n)

  ## START ALGORITHM
  
  if(!is.null(start)) {
    start <- c(start, rep(0,k))
    start <- start[1:k]
  }

  if(k>2 && is.null(start)) start <- rep(var(y,na.rm=T),k)
  if(k==1 && is.null(start)) start <- var(y,na.rm=T)
  
  if(is.null(start) && k==2) {
    if(missing(gamVals)) {
      gamVals <- seq(0.01,0.02,length=3)^2
      gamVals <- sort(c(gamVals,seq(0.1,0.9,length=3),1-gamVals))
      gamVals <- 0.5
    }
    if(length(gamVals)>1) {
      if(print.level>=1) cat("Evaluating the llik at gamma = \n")
      if(print.level>=1) cat(gamVals)
      if(print.level>=1) cat("\n")
      reg.obj <- reml(gamVals,y,X,V[[1]],V[[2]],print.level=print.level)
      llik <- reg.obj$llik
      llik <- as.real(llik)
      if(print.level>=2) cat(llik,"\n")
      gam <- gamVals[llik==max(llik)]
      gam <- gam[1]
      if(print.level>=2) cat("MLE is near",gam,"and llik =",max(llik),"there\n")
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
    if(print.level>=1) cat(c("start algorithm at",round(start,4),"\n"))
  }

  sigma <- start
  ## reparameterise so everything will get into the correct spot after exp
  for(i in 1:k)
    {
      if(pos[i]) sigma[i] <- log(start[i])
    }
  
  for(cycle in 1:maxcyc){

    ## newton-raphson step
    ## at the end of the iteration x will be -l'(sigma)/l''(sigma)
    ## hence we add instead of subtract
    ## fraction controls the proportion of each step to take
    ## for some reason 0.5 works very well
    
    sigma <- sigma+fraction*x

    if(print.level>=1) {
      cat(cycle, " ")
      for(i in 1:k) {
        if(pos[i]) {
          cat(exp(sigma[i])," ")
        } else {
          cat(sigma[i]," ")
        }
      }
    }
    
    ##Sigma <- matrix(0,dim(V[[1]])[1],dim(V[[1]])[2])
    Sigma <- 0
    for(i in 1:k)
      {
        if(pos[i]) {
          Sigma <- Sigma + V[[i]]*exp(sigma[i])
        } else {
          Sigma <- Sigma + V[[i]]*sigma[i]
        
	}
      }

    WX <- solve(Sigma,cbind(X,In))
    W <- WX[,(dim(X)[2]+1):dim(WX)[2]]
    WX <- WX[,1:dim(X)[2]]
    XtWX <- t(X)%*%WX
    WQ <- W - WX%*%solve(XtWX,t(WX))
    rss <- as.numeric(t(y) %*% WQ %*% y)
    ldet <- sum(log(eigen(WQ,symmetric=TRUE,only.values=TRUE)$values[1:rankQ]))
    
    ## REML LLIK
    rllik1 <- ldet/2 - rss/2
    
    ## MARGINAL LLIK BASED ON RESIDUAL CONFIGURATION STATISTIC
    rllik2 <- ldet/2 - rankQ * log(rss)/2
    
    if(print.level) cat("REML =",rllik1,"\n")

    ## change in REML llik
    if(cycle > 1) delta <- rllik1 - rllik0
    rllik0 <- rllik1
    
    ## now the fun starts, derivative and expected fisher info
    ## the 0.5 multiple is ignored, it is in both and they cancel
    
    T <- list(NULL)
    x <- NULL

    ## derivatives are now D[[i]] = var.components[i]*V[[i]]
    var.components <- rep(1,k)
    for(i in 1:k)
      {
        if(pos[i]) var.components[i] <- exp(sigma[i])
      }

    for(i in 1:k)
      {
        ##if(identity && i==1 && pos[i]) {
        ##  T[[i]] <- WQ
        ##} else
        T[[i]] <- WQ %*% V[[i]]
        x[i] <- as.numeric(t(y) %*% T[[i]] %*% WQ %*% y - sum(diag(T[[i]])))
        x[i] <- x[i]*var.components[i]
        
        for(j in c(1:i))
          {
            ## expected fisher information
            ##A[j,i] <- Tr(T[[j]] %*% T[[i]])
            ## the Ts are not symmetric, hence the transpose below
            A[j,i] <- sum(T[[j]] * t(T[[i]]))
            A[j,i] <- A[j,i]*var.components[i]*var.components[j]
            A[i,j] <- A[j,i]
          }
      }
    
    stats <- c(stats, rllik1, rllik2, sigma[1:k], x[1:k])
    if(print.level==-1) {
      cat(c(rllik1, rllik2, sigma[1:k], x[1:k]),"\n")
    }

    A.svd <- ginv(A)
    x <- A.svd %*% x

    if(qr(A)$rank < k){
      if(cycle==1) {
        if(print.level) {
          cat("Warning: Non identifiable dispersion model\n")
          ##print(round(A,6))
          for(i in 1:k) {
            if(pos[i]) {
              cat(exp(sigma[i])," ")
            } else {
              cat(sigma[i], " ")
            }
          }
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
      if(print.level) cat("WARNING:  maximum number of cycles reached before convergence\n")
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
    
    predicted <- fitted.values + (Sigma - gam*In) %*% W%*%(y - fitted.values)
  }

  result <- list(trace=llik, llik=llik[cycle, 1], cycle=cycle,
                 rdf=rankQ, beta=beta, beta.cov=beta.cov, beta.se=beta.se,
                 sigma=sigma[1:k], sigma.cov=(A.svd[1:k, 1:k] * 2), W=W, Q=Q,
                 fitted=fitted.values, predicted=predicted, pos=pos,
                 Vnames=Vcoef.names)
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

summary.regress <- function(object, digits=3, fixed.effects=T,...)
  {
    cat("\nMaximised Residual Log Likelihood is",round(object$llik,digits),"\n",sep=" ")
    indent.lin <- max(nchar(dimnames(object$beta)[[1]]))
    indent.var <- max(nchar(object$Vnames))
    indent <- max(indent.lin,indent.var)

    extra.space <- ""
    space.var <- extra.space
    for(i in 0:(indent-indent.var)) space.var <- paste(space.var," ",sep="")
    space.lin <- extra.space
    for(i in 0:(indent-indent.lin)) space.lin <- paste(space.lin," ",sep="")
    
    coefficients <- cbind(object$beta,object$beta.se)
    dimnames(coefficients)[[2]] <- c("Estimate","Std. Error")
    coefficients <- round(coefficients,digits)
    if(fixed.effects) {
      cat("\nLinear Coefficients:\n")
      row.names(coefficients) <- paste(space.lin,dimnames(object$beta)[[1]],sep="")
      print(coefficients)
      cat("\n")
    } else {
      cat("\nLinear Coefficients: not shown\n\n")
    }
    
    if(sum(object$pos)==0) {
      var.coefficients <- cbind(object$sigma,sqrt(diag(as.matrix(object$sigma.cov))))
      ##row.names(var.coefficients) <- paste(space.var,"coef.",c(1:length(object$sigma)),sep="")
      row.names(var.coefficients) <- paste(space.var,object$Vnames,sep="")
      dimnames(var.coefficients)[[2]] <- c("Estimate","Std. Error")
      var.coefficients <- round(var.coefficients,digits)
      cat("Variance Coefficients:\n")
      print(var.coefficients)
      cat("\n")
    } else {
      var.coefficients <- cbind(object$sigma, object$sigma - 2*sqrt(diag(as.matrix(object$sigma.cov))), object$sigma + 2*sqrt(diag(as.matrix(object$sigma.cov))))
      for(i in 1:length(object$pos))
        {
          if(object$pos[i]) var.coefficients[i,] <- exp(var.coefficients[i,])
        }
      ##row.names(var.coefficients) <- paste(space.var,c(1:length(object$sigma)),sep="")
      row.names(var.coefficients) <- paste(space.var,object$Vnames,sep="")
      dimnames(var.coefficients)[[2]] <- c("Estimate","Conf. Low.", "Conf. Upp.")
      var.coefficients <- round(var.coefficients,digits)
      cat("Variance Coefficients:\n")
      print(var.coefficients)
      cat("\n")
    }
  }

## when two matrices are passed to regress this is also called
## to evaluate the REML at certain values of gamma and find a
## good place to start the regress algorithm

reml <- function(lambda, y, X, V0, V1,print.level=0){

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
  ##if(print.level==1) cat("n-p =",n,"-",qr$rank,"=",rank,"\n")
  for(i in 1:length(lambda))
    {
      if(print.level>=2) cat(lambda[i],"\n")
      
      Sigma <- (1-lambda[i])*V0 + lambda[i] * V1
      WX <- solve(Sigma,cbind(X,In))
      W <- WX[,(dim(X)[2]+1):dim(WX)[2]]
      WX <- WX[,1:dim(X)[2]]
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
  WX <- solve(Sigma,cbind(X,In))
  W <- WX[,(dim(X)[2]+1):dim(WX)[2]]
  WX <- WX[,1:dim(X)[2]]
  XtWX <- t(X)%*%WX
  FItWX <- solve(XtWX,t(WX))
  WQ <- W - WX%*%FItWX
  rss <- t(y) %*% WQ %*% y
  beta <- FItWX %*% y

  list(llik=as.numeric(llik),rms=rss/rank, beta=beta, gamma=lambda, gamMax=lambdamax)
}

