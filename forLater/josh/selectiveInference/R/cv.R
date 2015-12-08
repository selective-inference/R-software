# ------------------------------------------------
# Cross-validation, preliminary

cvMakeFolds <- function(x, nfolds = 10) {
    #inds <- sample(1:nrow(x), replace=FALSE)
    inds <- 1:nrow(x)
    foldsize <- floor(nrow(x)/nfolds)
    lapply(1:nfolds, function(f) return(inds[1:foldsize+(f-1)*foldsize]))
}

############################################
# Can this be optimized using svdu_thresh? #
############################################
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

    RSSquads <- list()
    for (s in 1:maxsteps) {
        initial.active <- lapply(active.sets, function(a) a[1:s])
        RSSquads[[s]] <- cvRSSquad(X, folds, initial.active)
    }

    RSSs <- lapply(RSSquads, function(Q) t(Y) %*% Q %*% Y)
    sstar <- which.min(RSSs)
    quadstar <- RSSquads[sstar][[1]]

    RSSquads <- lapply(RSSquads, function(quad) quad - quadstar)
    RSSquads[[sstar]] <- NULL # remove the all zeroes case

    fit <- groupfs(X, Y, index=index, maxsteps=sstar, sigma=sigma, intercept=intercept, center=center, normalize=normalize)
    fit$cvobj <- cvobj
    fit$cvquad <- RSSquads

    fit$cvperm <- cv_perm

    invisible(fit)
}

