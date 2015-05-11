#########
#
# Sampling multivariate normal under affine constraints
#
#
# Original code in Pyx -  Jonathan Taylor (sample_truncnorm)
# Ported from Pyx to R -  Yuval Benjamini
#
#
#
#
##########



sample_truncnorm = function(A, b, initial, eta,
                                    Sigma = diag(ncol(A)),
                                    mu = rep(0,ncol(A)),
                                    how_often=1000,
                                    burnin=500,
                                    ndraw=1000,
                                    thinning=1,
                                    use_A=TRUE) {
   # A the linear constraints of Ax <= b
   # b The truncations of Ax <= b
   # mean, sigma, the parameters of the multivariate Gaussian to be truncated
   # 
  
  if (is.null(ncol(Sigma))) { 
    sigma = sqrt(Sigma)
    new_b = b - A * mu
    new_initial = initial - mu
    Z = sample_truncnorm_white(A = A,  b = new_b,
                                           initial = initial, 
                                           eta = eta,
                                           how_often=how_often,
                                           ndraw=ndraw, 
                                           burnin=burnin,
                                           sigma= sigma, # the sd, not the variance
                                           thinning = thinning,
                                           use_A=use_A)
    Z = Z + mu
    
  } else {
  
    whitened = whiten(A = A, b = b, Sigma = Sigma, mu = mu)
    initial_new = whitened$forward_map(initial)
    eta_new = whitened$forward_map(eta)

    white_samples = sample_truncnorm_white(A = whitened$A,
                               b = whitened$b,
                               initial = initial_new, 
                               eta = eta_new,
                               how_often=how_often,
                               ndraw=ndraw, 
                               burnin=burnin,
                               sigma=1.,
                               thinning = thinning,
                               use_A=use_A)
    
    
    Z = t(whitened$inverse_map(t(white_samples)))
  }
  return(Z)  
}

thresh2constraints = function(d, lower = rep(-Inf, d), upper = rep(Inf,d)){
  # Converts a marginal constraint on vector x in R^d 
  #               lower_i < x_i <upper_i 
  # into a linear-constraints representation Ax <= b. 
  #
  # lower and upper can have -Inf or Inf coordinates 
  
  stopifnot(is.element(length(lower),c(1,d)))
  stopifnot(is.element(length(upper),c(1,d)))
  
  if (length(lower) == 1){
    lower = rep(lower, d)
  }
  if (length(upper) == 1){
    upper = rep(upper, d)
  }
  
  
  A = matrix(nc = d, nr = 0)
  b = numeric(0)
  lower_constraints = which(lower > -Inf)
  for (l in lower_constraints){
    new_vec = rep(0,d)
    new_vec[l] = -1
    A = rbind(A, new_vec)
    b = c(b, -lower[l])
  }
  upper_constraints = which(upper < Inf)
  for (u in upper_constraints){
    new_vec = rep(0,d)
    new_vec[u] = 1
    A = rbind(A, new_vec)
    b = c(b, upper[u])
  }
  
  constraints = list(A = A, b = b)
  return(constraints)
}



whiten = function(A, b, Sigma, mu){
  
  #     Return a whitened version of constraints in a different
  #      basis, and a change of basis matrix.
  
  #If `self.covariance` is rank deficient, the change-of
  #basis matrix will not be square.

  require("Matrix")
  rank = rankMatrix(Sigma)
  eigSigma = eigen(Sigma)
  D = sqrt(eigSigma$values[1:rank])
  U = eigSigma$vectors[,1:rank]
  
  sqrt_cov = U 
  sqrt_inv = U
  for (i in 1:rank){ 
    sqrt_cov[,i] = sqrt_cov[,i] * D[i]
    sqrt_inv[,i] = sqrt_inv[,i] / D[i]
  }
  sqrt_inv = t(sqrt_inv)
   
  white = list()
  white$A = A %*% sqrt_cov
  white$b = b - A %*% mu
  
  white$inverse_map = function(x) {sqrt_cov %*% x + mu}
  white$forward_map = function(x) {sqrt_inv %*% (x - mu) } 
  return (white)
}


sample_truncnorm_white =  function(A, b, initial, eta,
                                   how_often=1000,
                                   sigma=1.,
                                   burnin=500,
                                   ndraw=1000,
                                   thinning=1,
                                   use_A=TRUE) {
  
  #     Sample from a truncated normal with covariance
  #     equal to sigma**2 I.
  #     
  #     Constraint is $Ax \leq b$ where `A` has shape
  #     `(q,n)` with `q` the number of constraints and
  #     `n` the number of random variables.
  # 
  # 
  #     Parameters
  #     ----------
  # 
  #     A : np.float((q,n))
  #        Linear part of affine constraints.
  
  #     b : np.float(q)
  #         Offset part of affine constraints.
  # 
  #     initial : np.float(n)
  #         Initial point for Gibbs draws.
  #         Assumed to satisfy the constraints.
  # 
  #     bias_direction : np.float (optional)
  #         Which projection is of most interest?
  # 
  #     how_often : int (optional)
  #         How often should the sampler make a move along `direction_of_interest`?
  #         If negative, defaults to ndraw+burnin (so it will never be used).
  # 
  #     sigma : float
  #         Variance parameter.
  # 
  #     burnin : int
  #         How many iterations until we start
  #         recording samples?
  # 
  #     ndraw : int
  #         How many samples should we return?
  # 
  #     Returns
  #     -------
  # 
  #     trunc_sample : np.float((ndraw, n))
  
  
  nvar = ncol(A)
  nconstraint = nrow(A)
  
  trunc_sample = matrix(NA, nr = ndraw, nc = nvar)
  state = initial
  
  irow = 1
  ivar = 1
  
  lower_bound = 0
  upper_bound = 0
  V = 0
  val = 0
  alpha = 0
  
  tol = 1.e-7
  
  U = A %*% state - b
#  usample = runif(burnin + ndraw)
  
  # directions not parallel to coordinate axes
  
  rand_dirs = ceiling(nvar/5)
  if (use_A){
    directions = rbind(A,matrix(rnorm(rand_dirs*nvar,mean = 0, sd =1),nr = rand_dirs, nc = nvar))
  } else { #use_A
    directions = matrix(rnorm(nvar^2,mean = 0, sd =1),nc = nvar, nr = nvar)
  }
  
  ndir = nrow(directions)
  directions[ndir,] = eta  # Set eta as the last direction 
  
  # Normalize the rows 
  for (i in 1:ndir) {
    directions[i,] = directions[i,] / sqrt(sum(directions[i,]^2))
  }
  
  
  alphas_dir = A %*% t(directions)
  alphas_coord = A
  alphas_max_dir = apply(abs(alphas_dir), 2, max) * tol    
  alphas_max_coord = apply(abs(alphas_coord), 2, max) * tol
  
  # for switching between coordinate updates and
  # other directions
  
  invperiod = 13
  docoord = FALSE
  iperiod = 0
  ibias = 0
  dobias = FALSE
  
  samp_idx = 1
  
  for (iter_count in 1:((ndraw+1)*(thinning) + burnin)){
    
    docoord = 1
    iperiod = iperiod + 1
    ibias = ibias + 1
    
    if (iperiod == invperiod){ 
      docoord = FALSE
      iperiod = 1
      dobias = FALSE
    }
    
    if (ibias == how_often){
      docoord = FALSE
      ibias = 1
      dobias = TRUE
    }
    
    if (docoord == 1) {
      idx = sample(nvar, size = 1)
      V = state[idx]
    } else { # docoord == 1
      if (!dobias){
        idx = sample(ndir,size = 1)
      } else { # !dobias 
        idx = ndir # last row of directions is bias_direction
      }
      V = sum(directions[idx,] * state)
    }
    
    lower_bound = -1e12
    upper_bound = 1e12
    for (irow in 1:nconstraint){ #anyway to makke this faster?
      if (docoord) { 
        alpha = alphas_coord[irow,idx]
        val = -U[irow] / alpha + V
        if ((alpha > alphas_max_coord[idx]) & (val < upper_bound)) {
          upper_bound = val
        } else if ((alpha < -alphas_max_coord[idx]) & (val > lower_bound)){
          lower_bound = val
        }
      } else { #docoord
        alpha = alphas_dir[irow,idx]
        val = -U[irow] / alpha + V
        if ((alpha > alphas_max_dir[idx]) & (val < upper_bound)) {
          upper_bound = val
        } else if ((alpha < -alphas_max_dir[idx]) & (val > lower_bound)){
          lower_bound = val
        }
      }
    }
    if (lower_bound > V){
      lower_bound = V - tol * sigma
    } 
    if (upper_bound < V) { 
      upper_bound = V + tol * sigma
    }
    
    lower_bound = lower_bound / sigma
    upper_bound = upper_bound / sigma
    
    stopifnot(lower_bound <= upper_bound ) #  "bound violation" 
    
    if (upper_bound < -10){ # use Exp approximation
      unif = runif(1) * (1 - exp((lower_bound - upper_bound) * upper_bound))
      tnorm = upper_bound - log(1 - unif) / upper_bound 
    } else if (lower_bound > 10) { 
      unif = runif(1) * (1 - exp(-(upper_bound - lower_bound) * lower_bound))
      tnorm = -log(1 - unif) / lower_bound + lower_bound
    } else if (lower_bound < 0) { 
      cdfL = pnorm(lower_bound)
      cdfU = pnorm(upper_bound)
      unif = runif(1) * (cdfU - cdfL) + cdfL
      if (unif < 0.5){
        tnorm = qnorm(unif) * sigma
      } else { # unif < 0.5
        tnorm = -qnorm(1-unif) * sigma
      }
    } else {  # upper_bound < -10, lower_bound > 10, lower_bound < 0 
      cdfL = pnorm(-lower_bound)
      cdfU = pnorm(-upper_bound)
      unif = runif(1) * (cdfL - cdfU) + cdfU
      if (unif < 0.5) {
        tnorm = -qnorm(unif) * sigma
      } else { #unif < 0.5
        tnorm = qnorm(1-unif) * sigma
      }
    }
    
    if (docoord == 1) { 
      state[idx] = tnorm
      tnorm = tnorm - V
      U = U + tnorm * A[, idx] 
    } else {
      tnorm = tnorm - V
      for (ivar in 1:nvar) {
        state[ivar] = state[ivar] + tnorm * directions[idx,ivar]
        U = (U + A[, ivar] * 
                       tnorm * directions[idx,ivar])
      }
    }
    
    if ((iter_count >= burnin) & ((iter_count %% thinning) ==  0)){
      trunc_sample[samp_idx, ] = state
      samp_idx = samp_idx + 1
      if (samp_idx > ndraw) {
        break
      }
    }
  }
  return (trunc_sample)
}

