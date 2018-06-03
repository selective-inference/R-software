# Cross-validation

#--------------------------------------
# Functions for creating cv folds
# -------------------------------------

cvMakeFolds <- function(x, nfolds = 5) {
    inds <- sample(1:nrow(x), replace=FALSE)
    #inds <- 1:nrow(x)
    foldsize <- floor(nrow(x)/nfolds)
    folds <- lapply(1:nfolds, function(f) return(inds[1:foldsize+(f-1)*foldsize]))
    if (nfolds*foldsize < nrow(x)) {
      # remainder observations added to first several folds
      for (i in 1:(nrow(x) - nfolds*foldsize)) {
        folds[[i]] <- c(folds[[i]], inds[nfolds*foldsize + i])
      }
    }
    return(folds)
}

# To interface with glmnet
foldidglmnet <- function(folds) {
  n <- sum(sapply(folds, length))
  glmnetfoldid <- rep(0, n)
  for (ind in 1:length(folds)) {
    glmnetfoldid[folds[[ind]]] <- ind
  }
  glmnetfoldid
}

# cv.glmnet and estimateSigma mashup
cvglmnetlar <- function(x, y, foldid) {
  cvfit <- cv.glmnet(x, y, intercept = FALSE, foldid = foldid)
  lamhat <- cvfit$lambda.min
  yhat <- predict(cvfit, x, s = lamhat)
  nz <- sum(coef(cvfit, s = lamhat) !=0)
  cvfit$sigma <- sqrt(sum((y-yhat)^2)/(length(y)-nz-1))
  cvfit$df <- nz
  return(cvfit)
}


#--------------------------------------
# Functions for computing quadratic form for cv-error
#--------------------------------------

# There seems to be a problem here, picking overly conservative models

# Can this be optimized using svdu_thresh? 
cvHatMatrix <- function(x, folds, active.sets) {
    nfolds <- length(folds)
    lapply(1:nfolds, function(f) {
        fold <- folds[[f]]
        active <- active.sets[[f]]
        x_tr <- x[ -fold, active]
        x_te <- x[fold, active]
        hatm <- matrix(0, nrow=length(fold), ncol=nrow(x))
        svdtr <- svd(x_tr)
        inds <- svdtr$d > .Machine$double.eps * svdtr$d[1]
        xtrinv <- svdtr$v[, inds, drop = FALSE] %*% ((1/svdtr$d[inds]) * t(svdtr$u[, inds, drop = FALSE]))
        hatm[, -fold] <- x_te %*% xtrinv
        return(hatm)
    })
}

cvProductHat <- function(folds, inds, finds, ginds, hat_matrices) {
    nfolds <- length(folds)
    terms <- lapply(inds, function(h) {
        t(hat_matrices[[h]][, finds]) %*% hat_matrices[[h]][, ginds]
    })
    return(Reduce('+', terms))
}

# This is too "clever," I can't easily understand it
# simpler code is preferable for maintenance and forking etc
cvRSSquad <- function(x, folds, active.sets) {
    hat_matrices <- cvHatMatrix(x, folds, active.sets)
    nfolds <- length(folds)
    rows <- lapply(1:nfolds, function(f) {
        do.call(cbind, lapply(1:nfolds, function(g) {
            ginds <- folds[[g]]
            finds <- folds[[f]]
            if (f == g) {
                return(cvProductHat(folds, setdiff(1:nfolds, f), finds, ginds, hat_matrices))
            } else {
                return(
                    cvProductHat(folds, setdiff(1:nfolds, c(f,g)), finds, ginds, hat_matrices) - hat_matrices[[f]][, ginds] - t(hat_matrices[[g]][, finds]))
            }
        }))
    })
    Q <- do.call(rbind, rows)
    return(Q)
}

cvopt <- function(x, y, maxsteps, folds, active.sets) {
  yperm <- y[order(unlist(folds))]
  RSSquads <- list()
  # Can this loop be optimized with smart updating of each model along each path?
  for (s in 1:maxsteps) {
    initial.active <- lapply(active.sets, function(a) a[1:s])
    RSSquads[[s]] <- cvRSSquad(x, folds, initial.active)
  }
  
  RSSs <- lapply(RSSquads, function(Q) t(y) %*% Q %*% y)
  sstar <- which.max(RSSs)
  quadstar <- RSSquads[sstar][[1]]
  
  RSSquads <- lapply(RSSquads, function(quad) quad - quadstar)
  RSSquads[[sstar]] <- NULL # remove the all zeroes case
  return(list(sstar = sstar, RSSquads = RSSquads))
}


#--------------------------------------
# Functions for forward stepwise
# broke this while making cvlar
#--------------------------------------

cvfs <- function(x, y, index = 1:ncol(x), maxsteps, sigma = NULL, intercept = TRUE, center = TRUE, normalize = TRUE, nfolds = 5) {

    n <- nrow(x)
    if (maxsteps >= n*(1-1/nfolds)) {
        maxsteps <- floor(n*(1-1/nfolds))
        warning(paste("maxsteps too large for training fold size, set to", maxsteps))
    }

    folds <- cvMakeFolds(x, nfolds)
    nfolds <- length(folds)
    projections <- list(1:nfolds)
    maxprojs <- list(1:nfolds)
    active.sets <- list(1:nfolds)
    cvobj <- list(1:nfolds)
    cv_perm <- sample(1:n)
    Y <- y[cv_perm]
    X <- x[cv_perm, ]

    # Initialize copies of data for loop
    by <- mean(Y)
    if (intercept) Y <- Y - by

    # Center and scale design matrix
    xscaled <- scaleGroups(X, index, center, normalize)
    xm <- xscaled$xm
    xs <- xscaled$xs
    X <- xscaled$x

    # Flatten list or something?
    for (f in 1:nfolds) {
        fold <- folds[[f]]
        fit <- groupfs(X[-fold,], Y[-fold], index=index, maxsteps=maxsteps, sigma=sigma, intercept=FALSE, center=FALSE, normalize=FALSE)
        fit$fold <- fold
        # Why is this commented out?
        ## projections[[f]] <- lapply(fit$projections, function(step.projs) {
        ##     lapply(step.projs, function(proj) {
        ##         # Reduce from n by n matrix to svdu_thresh
        ##         expanded.proj <- matrix(0, n, ncol(proj))
        ##         expanded.proj[-fold, ] <- proj
        ##         return(expanded.proj)
        ##     })
        ## })
        active.sets[[f]] <- fit$action
        cvobj[[f]] <- fit
    }
    #projections <- do.call(c, projections)



    fit <- groupfs(X, Y, index=index, maxsteps=sstar, sigma=sigma, intercept=intercept, center=center, normalize=normalize)
    fit$cvobj <- cvobj
    fit$cvquad <- RSSquads

    fit$cvperm <- cv_perm

    invisible(fit)
}


#--------------------------------------
# Functions for lar
#--------------------------------------

cvlar <- function(x, y, maxsteps, folds = NULL) { # other args
  this.call = match.call()
  if (is.null(folds)) folds <- cvMakeFolds(x)
  models <- lapply(folds, function(fold) {
    x.train <- x
    y.train <- y
    x.train[fold,] <- 0
    y.train[fold] <- 0
    x.test <- x[fold,]
    y.test <- y[fold]
    larpath.train <- lar(x.train, y.train, maxsteps = maxsteps, intercept = F, normalize = F)
    return(larpath.train)
  })
  
  active.sets <- lapply(models, function(model) model$action)
  #lambdas <- lapply(models, function(model) model$lambda)
  #lmin <- min(unlist(lambdas))
  cvmin <- cvopt(x, y, maxsteps, folds, active.sets)
  sstar <- cvmin$sstar
  fit <- lar(x, y, maxsteps=sstar, intercept = F, normalize = F)
  fit$ols <- lsfit(x[, fit$action, drop = F], y, intercept = F)
  names(fit$ols$coefficients) <- fit$action
  fit$sigma <- sqrt(sum((fit$ols$residuals)^2)/(length(y)-length(fit$action)-1))
  fit$RSSquads <- cvmin$RSSquads
  # tall Gamma encoding all cv-model paths
  fit$tallGamma <- do.call(rbind, lapply(models, function(model) return(model$Gamma)))
  fit$khat <- sstar
  fit$folds <- folds
  fit$call <- this.call
  class(fit) <- "cvlar"
  # more to do here?
  return(fit)
}

# cvlarInf <- function(obj, ...) {
#   pv.unadj <- larInf(obj, type = "all", k = obj$khat, verbose = T, ...)
#   obj$Gamma <- rbind(obj$Gamma, obj$tallGamma)
#   pv.adj <- larInf(obj, type = "all", k = obj$khat, verbose = T, ...)
#   return(list(pv.unadj = pv.unadj, pv.adj = pv.adj))
# }

cvlarInf <- function (obj, sigma, alpha = 0.1,
                    k = NULL,
                    gridrange = c(-100, 100),
                    bits = NULL, mult = 2, 
                    ntimes = 2, verbose = FALSE) {
  this.call = match.call()
  #checkargs.misc(sigma = sigma, alpha = alpha, k = k, gridrange = gridrange, mult = mult, ntimes = ntimes)
  if (class(obj) != "cvlar") 
    stop("obj must be an object of class cvlar")
  if (!is.null(bits) && !requireNamespace("Rmpfr", quietly = TRUE)) {
    warning("Package Rmpfr is not installed, reverting to standard precision")
    bits = NULL
  }
  x = obj$x
  y = obj$y
  p = ncol(x)
  n = nrow(x)
  G = obj$Gamma
  #nk = obj$nk
  sx = obj$sx
  k = obj$khat
  sigma = obj$sigma
  # may the gods of OOP have mercy on us
  class(obj) <- "lar"
  pv.unadj <- larInf(obj, type = "all", sigma = sigma, k = obj$khat)
  class(obj) <- "cvlar"
  #pv.spacing = pv.modspac = pv.covtest = khat = NULL
  
  G = rbind(obj$Gamma, obj$tallGamma) #G[1:nk[k], ]
  u = rep(0, nrow(G))
  kk = k
  pv = vlo = vup = numeric(kk)
  vmat = matrix(0, kk, n)
  ci = tailarea = matrix(0, kk, 2)
  sign = numeric(kk)
  vars = obj$action[1:kk]
  xa = x[, vars]
  M = pinv(crossprod(xa)) %*% t(xa)
  for (j in 1:kk) {
    if (verbose) 
      cat(sprintf("Inference for variable %i ...\n", 
                  vars[j]))
    vj = M[j, ]
    mj = sqrt(sum(vj^2))
    vj = vj/mj
    sign[j] = sign(sum(vj * y))
    vj = sign[j] * vj
    Gj = rbind(G, vj)
    uj = c(u, 0)
    a = poly.pval(y, Gj, uj, vj, sigma, bits)
    pv[j] = a$pv
    sxj = sx[vars[j]]
    vlo[j] = a$vlo * mj/sxj
    vup[j] = a$vup * mj/sxj
    vmat[j, ] = vj * mj/sxj
    
    #a = poly.int(y, Gj, uj, vj, sigma, alpha, gridrange = gridrange, flip = (sign[j] == -1), bits = bits)
    #ci[j, ] = a$int * mj/sxj
    #tailarea[j, ] = a$tailarea
  }
  out = list(type = type, k = k, khat = khat, pv = pv, 
             pv.unadj = pv.unadj, vlo = vlo, vup = vup, vmat = vmat, 
             y = y, vars = vars, sign = sign, sigma = sigma, 
             alpha = alpha, call = this.call)
  class(out) = "cvlarInf"
  return(out)
}



poly.pval <- function(y, G, u, v, sigma, bits=NULL) {
  z = sum(v*y)
  vv = sum(v^2)
  sd = sigma*sqrt(vv)
  
  rho = G %*% v / vv
  vec = (u - G %*% y + rho*z) / rho
  vlo = suppressWarnings(max(vec[rho>0]))
  vup = suppressWarnings(min(vec[rho<0]))
  
  pv = tnorm.surv(z,0,sd,vlo,vup,bits)
  return(list(pv=pv,vlo=vlo,vup=vup))
}

pinv <- function(A, tol=.Machine$double.eps) {
  e = eigen(A)
  v = Re(e$vec)
  d = Re(e$val)
  d[d > tol] = 1/d[d > tol]
  d[d < tol] = 0
  if (length(d)==1) return(v*d*v)
  else return(v %*% diag(d) %*% t(v))
}

tnorm.surv <- function(z, mean, sd, a, b, bits=NULL) {
  z = max(min(z,b),a)
  
  # Check silly boundary cases
  p = numeric(length(mean))
  p[mean==-Inf] = 0
  p[mean==Inf] = 1
  
  # Try the multi precision floating point calculation first
  o = is.finite(mean)
  mm = mean[o]
  pp = mpfr.tnorm.surv(z,mm,sd,a,b,bits) 
  
  # If there are any NAs, then settle for an approximation
  oo = is.na(pp)
  if (any(oo)) pp[oo] = bryc.tnorm.surv(z,mm[oo],sd,a,b)
  
  p[o] = pp
  return(p)
}

mpfr.tnorm.surv <- function(z, mean=0, sd=1, a, b, bits=NULL) {
  # If bits is not NULL, then we are supposed to be using Rmpf
  # (note that this was fail if Rmpfr is not installed; but
  # by the time this function is being executed, this should
  # have been properly checked at a higher level; and if Rmpfr
  # is not installed, bits would have been previously set to NULL)
  if (!is.null(bits)) {
    z = Rmpfr::mpfr((z-mean)/sd, precBits=bits)
    a = Rmpfr::mpfr((a-mean)/sd, precBits=bits)
    b = Rmpfr::mpfr((b-mean)/sd, precBits=bits)
    return(as.numeric((Rmpfr::pnorm(b)-Rmpfr::pnorm(z))/
                        (Rmpfr::pnorm(b)-Rmpfr::pnorm(a))))
  }
  
  # Else, just use standard floating point calculations
  z = (z-mean)/sd
  a = (a-mean)/sd
  b = (b-mean)/sd
  return((pnorm(b)-pnorm(z))/(pnorm(b)-pnorm(a)))
}

bryc.tnorm.surv <- function(z, mean=0, sd=1, a, b) {
  z = (z-mean)/sd
  a = (a-mean)/sd
  b = (b-mean)/sd
  n = length(mean)
  
  term1 = exp(z*z)
  o = a > -Inf
  term1[o] = ff(a[o])*exp(-(a[o]^2-z[o]^2)/2)
  term2 = rep(0,n)
  oo = b < Inf
  term2[oo] = ff(b[oo])*exp(-(b[oo]^2-z[oo]^2)/2)
  p = (ff(z)-term2)/(term1-term2)
  
  # Sometimes the approximation can give wacky p-values,
  # outside of [0,1] ..
  #p[p<0 | p>1] = NA
  p = pmin(1,pmax(0,p))
  return(p)
}

ff <- function(z) {
  return((z^2+5.575192695*z+12.7743632)/
           (z^3*sqrt(2*pi)+14.38718147*z*z+31.53531977*z+2*12.77436324))
}



print.cvlar <- function(x, ...) {
  cat("\nCall:\n")
  dput(x$call)
  
  cat("\nSequence of LAR moves:\n")
  nsteps = length(x$action)
  tab = cbind(1:nsteps,x$action,x$sign)
  colnames(tab) = c("Step","Var","Sign")
  rownames(tab) = rep("",nrow(tab))
  print(tab)
  invisible()
}

print.cvlarInf <- function(x, ...) {
  cat("\nCall:\n")
  dput(x$call)
  
  cat(sprintf("\nStandard deviation of noise (specified or estimated) sigma = %0.3f\n",
              x$sigma))
  
  
  cat(sprintf("\nTesting results at step = %i, with alpha = %0.3f\n",x$k,x$alpha))
  cat("",fill=T)
  tab = cbind(x$vars,
              round(x$sign*x$vmat%*%x$y,3),
              round(x$sign*x$vmat%*%x$y/(x$sigma*sqrt(rowSums(x$vmat^2))),3),
              round(x$pv,3))
  colnames(tab) = c("Var", "Coef", "Z-score", "P-value")
  rownames(tab) = rep("",nrow(tab))
  print(tab)
  invisible()
}