library(selectiveInference)


## Approximates inverse covariance matrix theta
InverseLinfty <- function(sigma, n, resol=1.5, mu=NULL, maxiter=50, threshold=1e-10, verbose = TRUE) {
  isgiven <- 1;
  if (is.null(mu)){
    isgiven <- 0;
  }
  
  p <- nrow(sigma);
  M <- matrix(0, p, p);
  xperc = 0;
  xp = round(p/10);
  for (i in 1:p) {
    if ((i %% xp)==0){
      xperc = xperc+10;
      if (verbose) {
        print(paste(xperc,"% done",sep="")); }
    }
    if (isgiven==0){
      mu <- (1/sqrt(n)) * qnorm(1-(0.1/(p^2)));
    }
    mu.stop <- 0;
    try.no <- 1;
    incr <- 0;
    while ((mu.stop != 1)&&(try.no<10)){
      last.beta <- beta
      output <- InverseLinftyOneRow(sigma, i, mu, maxiter=maxiter, threshold=threshold)
      beta <- output$optsol
      iter <- output$iter
      if (isgiven==1){
        mu.stop <- 1
      }
      else{
        if (try.no==1){
          if (iter == (maxiter+1)){
            incr <- 1;
            mu <- mu*resol;
          } else {
            incr <- 0;
            mu <- mu/resol;
          }
        }
        if (try.no > 1){
          if ((incr == 1)&&(iter == (maxiter+1))){
            mu <- mu*resol;
          }
          if ((incr == 1)&&(iter < (maxiter+1))){
            mu.stop <- 1;
          }
          if ((incr == 0)&&(iter < (maxiter+1))){
            mu <- mu/resol;
          }
          if ((incr == 0)&&(iter == (maxiter+1))){
            mu <- mu*resol;
            beta <- last.beta;
            mu.stop <- 1;
          }
        }
      }
      try.no <- try.no+1
    }
    M[i,] <- beta;
  }
  return(M)
}

InverseLinftyOneRow <- function ( sigma, i, mu, maxiter=50, threshold=1e-10) {
  p <- nrow(sigma);
  rho <- max(abs(sigma[i,-i])) / sigma[i,i];
  mu0 <- rho/(1+rho);
  beta <- rep(0,p);
  
  #if (mu >= mu0){
  #  beta[i] <- (1-mu0)/sigma[i,i];
  #  returnlist <- list("optsol" = beta, "iter" = 0);
  #  return(returnlist);
  #}
  
  diff.norm2 <- 1;
  last.norm2 <- 1;
  iter <- 1;
  iter.old <- 1;
  beta[i] <- (1-mu0)/sigma[i,i];
  beta.old <- beta;
  sigma.tilde <- sigma;
  diag(sigma.tilde) <- 0;
  vs <- -sigma.tilde%*%beta;
  
  while ((iter <= maxiter) && (diff.norm2 >= threshold*last.norm2)){
    
    for (j in 1:p){
      oldval <- beta[j];
      v <- vs[j];
      if (j==i)
        v <- v+1;
      beta[j] <- SoftThreshold(v,mu)/sigma[j,j];
      if (oldval != beta[j]){
        vs <- vs + (oldval-beta[j])*sigma.tilde[,j];
      }
    }
    
    iter <- iter + 1;
    if (iter==2*iter.old){
      d <- beta - beta.old;
      diff.norm2 <- sqrt(sum(d*d));
      last.norm2 <-sqrt(sum(beta*beta));
      iter.old <- iter;
      beta.old <- beta;
      #if (iter>10)
      #  vs <- -sigma.tilde%*%beta;
    }

    # print(c(iter, maxiter, diff.norm2, threshold * last.norm2, threshold, mu))

  }
  
  returnlist <- list("optsol" = beta, "iter" = iter)
  return(returnlist)
}

SoftThreshold <- function( x, lambda ) {
  #
  # Standard soft thresholding
  #
  if (x>lambda){
    return (x-lambda);}
  else {
    if (x< (-lambda)){
      return (x+lambda);}
    else {
      return (0); }
  }
}


### Test

n = 100; p = 50

X = matrix(rnorm(n * p), n, p)
S = t(X) %*% X / n

mu = 7.791408e-02

tol = 1.e-12

rows = c(1:2)
A1 = debiasingMatrix(S, FALSE, n, rows, mu=mu, max_iter=1000, kkt_tol=tol, objective_tol=tol, parameter_tol=tol)
A2 = debiasingMatrix(S / n, FALSE, n, rows, mu=mu, max_iter=1000, kkt_tol=tol, objective_tol=tol, parameter_tol=tol)

B1 = debiasingMatrix(X, TRUE, n, rows, mu=mu, max_iter=1000, kkt_tol=tol, objective_tol=tol, parameter_tol=tol)
B2 = debiasingMatrix(X / sqrt(n), TRUE, n, rows, mu=mu, max_iter=1000, kkt_tol=tol, objective_tol=tol, parameter_tol=tol)

C1 = InverseLinfty(S, n, mu=mu, maxiter=1000)[rows,]
C2 = InverseLinfty(S / n, n, mu=mu, maxiter=1000)[rows,]

par(mfrow=c(2,3))

plot(A1[1,], C1[1,])
plot(A1[1,], B1[1,])
plot(B1[1,], C1[1,])

plot(A1[1,], A2[1,])
plot(B1[1,], B2[1,])
plot(C1[1,], C2[1,])

print(c('A', sum(A1[1,] == 0)))
print(c('B', sum(B1[1,] == 0)))
print(c('C', sum(C1[1,] == 0)))

## Are our points feasible

feasibility = function(S, soln, j, mu) {
     p = nrow(S)
     E = rep(0, p)
     E[j] = 1
     G = S %*% soln - E
     return(c(max(abs(G)), mu))
}

print(c('feasibility A', feasibility(S, A1[1,], 1, mu)))
print(c('feasibility B', feasibility(S, B1[1,], 1, mu)))
print(c('feasibility C', feasibility(S, C1[1,], 1, mu)))

active_KKT = function(S, soln, j, mu) {
     p = nrow(S)
     E = rep(0, p)
     E[j] = 1
     G = S %*% soln - E
     print(which(soln != 0))
     print(G[j])
     return(c(G[soln != 0] * sign(soln)[soln != 0], mu))
}

print(c('active_KKT A', active_KKT(S, A1[1,], 1, mu)))
print(c('active_KKT B', active_KKT(S, B1[1,], 1, mu)))
print(c('active_KKT C', active_KKT(S, C1[1,], 1, mu)))


