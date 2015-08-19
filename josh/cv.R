source("groupfs.R")

# ------------------------------------------------
# Cross-validation, preliminary


cv_make_folds <- function(x, nfolds = 10) {
    #inds <- sample(1:nrow(x), replace=FALSE)
    inds <- 1:nrow(x)
    foldsize <- floor(nrow(x)/nfolds)
    lapply(1:nfolds, function(f) return(inds[1:foldsize+(f-1)*foldsize]))
}

cv_hat_matrix <- function(x, folds, active.sets) {
    nfolds <- length(folds)
    lapply(1:nfolds, function(f) {
        fold <- folds[[f]]
        active <- active.sets[[f]]
        x_tr <- x[ -fold, active]
        x_te <- x[fold, active]
        hatm <- matrix(0, nrow=length(fold), ncol=nrow(x))
        hatm[, -fold] <- x_te %*% ginv(x_tr)
        return(hatm)
    })
}

product_cv_hat <- function(folds, inds, finds, ginds, hat_matrices) {
    nfolds <- length(folds)
    terms <- lapply(inds, function(h) {
        t(hat_matrices[[h]][, finds]) %*% hat_matrices[[h]][, ginds]
    })
    return(Reduce('+', terms))
}

cv_RSSquad <- function(x, folds, active.sets) {
    hat_matrices <- cv_hat_matrix(x, folds, active.sets)
    nfolds <- length(folds)
    rows <- lapply(1:nfolds, function(f) {
        do.call(cbind, lapply(1:nfolds, function(g) {
            ginds <- folds[[g]]
            finds <- folds[[f]]                
            if (f == g) {
                return(product_cv_hat(folds, setdiff(1:nfolds, f), finds, ginds, hat_matrices))
            } else {
                return(
                    product_cv_hat(folds, setdiff(1:nfolds, c(f,g)), finds, ginds, hat_matrices) - hat_matrices[[f]][, ginds] - t(hat_matrices[[g]][, finds]))
            }
        }))
    })
    Q <- do.call(rbind, rows)
    return(Q)
}

cv_fs <- function(x, y, steps, nfolds = 10) {
    
    n <- nrow(x)
    if (steps >= n*(1-1/nfolds)) stop("Too many steps")
    
    folds <- cv_make_folds(x, nfolds)
    nfolds <- length(folds)
    projections <- list()
    active.sets <- list()
    cv_perm <- sample(1:n)
    Y <- y[cv_perm]
    X <- x[cv_perm, ]
    
    for (f in 1:nfolds) {
        fold <- folds[[f]]
        fit <- groupfs(X[-fold,], Y[-fold], steps=steps)
        path.projs <- fit$projections
        path.projs <- lapply(path.projs, function(step.projs) {
            lapply(step.projs, function(proj) {
                expanded.proj <- matrix(0, n, n)
                expanded.proj[-fold, -fold] <- proj
                return(expanded.proj)
            })
        })
        projections[[f]] <- path.projs
        active.sets[[f]] <- fit$variable
    }
    projections <- do.call(c, projections)

    RSSquads <- list()
    for (s in 1:steps) {
        initial.active <- lapply(active.sets, function(a) a[1:s])
        RSSquads[[s]] <- cv_RSSquad(X, folds, initial.active)
    }

    RSSs <- lapply(RSSquads, function(Q) t(y) %*% Q %*% y)
    sstar <- which.min(RSSs)
    quadstar <- RSSquads[sstar][[1]]

    RSSquads <- lapply(RSSquads, function(quad) quadstar - quad)

    fit <- groupfs(X, Y, steps=sstar)
    fit$projections <- c(fit$projections, projections, RSSquads)
    fit$cvperm <- cv_perm

    invisible(fit)
}

