# Functions for simulation/testing

randomGroupSizes <- function(G, lambda = 2) return(2 + rpois(G, lambda))

randomGroups <- function(G, lambda = 2) {
    rles <- randomGroupSizes(G, lambda)
    return(rep(1:G, rles))
}

randomIndexFixedP <- function(p, G) sort(c(sample(1:G), sample(1:G, size = p-G, replace=T)))

randomFactorDesign <- function(n, G, lambda = 2) {
    if (n < (1+lambda)*G) stop("Larger n required to avoid duplicate columns")
    rles <- randomGroupSizes(G, lambda)
    print(rles)
    df <- data.frame(do.call(cbind, lapply(rles, function(g) {
        sample(LETTERS[1:g], n, replace = TRUE, prob = runif(g))
    })), stringsAsFactors = TRUE)
    if (any(apply(df, 2, function(col) length(unique(col))) == 1)) return(randomFactorDesign(n, G, lambda))
    fd <- factorDesign(df)
    if (any(duplicated(fd$x, MARGIN = 2))) return(randomFactorDesign(n, G, lambda))
    return(list(df=df, fd=fd))
}

randomFactorsFixedP <- function(p, G) {
#    index <-
}

randomGaussianFixedP <- function(n, p, G = p, sparsity = 0, snr = 0, sigma = 1, rho = 0) {
    index <- 1:p
    if (G < p) index <- randomIndexFixedP(p, G)
    x <- matrix(rnorm(n*p), nrow=n)
    if (rho != 0) {
        z <- matrix(rep(t(rnorm(n)), p), nrow = n)
        x <- sqrt(1-rho)*x + sqrt(rho)*z
    }
    beta <- rep(0, p)
    if (sparsity > 0 && snr > 0) {
        for (j in 1:sparsity) {
            inds <- which(index == j)
            beta[inds] <- snr * sqrt(2*log(G)/(n*length(inds))) * sample(c(-1,1), length(inds), replace=T)
        }
    }
    y <- x %*% beta + sigma * rnorm(n)
    return(list(x=x, y=y, beta = beta, index=index, sigma = sigma))
}
