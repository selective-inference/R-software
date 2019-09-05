# solves full lasso problem via glmnet

solve_problem_glmnet = function(X, y, lambda_glmnet, penalty_factor, family){
  if (is.null(lambda_glmnet)){
    cv = cv.glmnet(x=X, 
                   y=y, 
                   family=family,
                   penalty.factor=penalty_factor, 
                   intercept=FALSE,  
                   thresh = 1e-12)
    beta_hat = coef(cv, s="lambda.min")
  }
  else {
    lasso = glmnet(x=X, 
                   y=y, 
                   family=family,
                   penalty.factor=penalty_factor,
                   alpha=1,
                   standardize=FALSE, 
                   intercept=FALSE, 
                   thresh=1e-20)
    beta_hat = coef(lasso, s=lambda_glmnet)
  }
  return(beta_hat[-1])
}

# solves the restricted problem
solve_restricted_problem = function(X, y, var, lambda_glmnet, penalty_factor, loss){
  restricted_soln=rep(0, ncol(X))
  restricted_soln[-var] = solve_problem_glmnet(X[,-var], 
                                               y, 
                                               lambda_glmnet, 
                                               penalty_factor[-var], 
                                               family=family_label(loss))
  return(restricted_soln)
}

solve_problem_Q = function(Xdesign, 
                           Qbeta_bar, 
                           lambda_glmnet, 
                           penalty_factor,
                           max_iter=50,
                           kkt_tol=1.e-4, 
                           objective_tol=1.e-4, 
                           parameter_tol=1.e-4,
                           kkt_stop=TRUE,
                           objective_stop=TRUE,	
                           parameter_stop=TRUE){
  n=nrow(Xdesign)
  p=ncol(Xdesign)

  linear_func = -as.numeric(Qbeta_bar)
  soln = as.numeric(rep(0., p))
  ever_active = as.integer(rep(0, p))
  ever_active[1] = 1
  ever_active = as.integer(ever_active)
  nactive = as.integer(1)
  Xsoln = as.numeric(rep(0, n))
  gradient = 1. * linear_func 
  max_active = as.integer(p)
  
  linear_func = linear_func/n
  gradient = gradient/n
  
  #solve_QP_wide solves n*linear_func^T\beta+1/2(X\beta)^T (X\beta)+\sum\lambda_i|\beta_i|

  result = solve_QP_wide(Xdesign,                                       # this is a design matrix
                         as.numeric(penalty_factor*lambda_glmnet),      # vector of Lagrange multipliers
                         0,                                             # ridge_term 
                         max_iter, 
                         soln, 
                         linear_func, 
                         gradient, 
                         Xsoln,
                         ever_active, 
                         nactive, 
                         kkt_tol, 
                         objective_tol, 
                         parameter_tol,
                         max_active,
                         kkt_stop,
                         objective_stop,	
                         parameter_stop)
  
  return(result$soln)
}


# the selection event is |sigma_est^2*(target_cov)^{-1}Z+center|>radius
# var is one of 1..p variables that we are buliding truncation set for

# JT: this function can be written to not depend on sigma_est -- see python code

truncation_set = function(X, 
                          y, 
                          Qbeta_bar, 
                          QE, 
                          Xdesign, 
                          target_stat, 
			  QiE,
                          var, 
                          active_set,
                          lambda_glmnet, 
                          penalty_factor, 
                          loss, 
                          solver){
  
  if (solver=="QP"){
    penalty_factor_rest = rep(penalty_factor)
    penalty_factor_rest[var] = 10^10
    restricted_soln = solve_problem_Q(Xdesign, 
                                      Qbeta_bar, 
                                      lambda_glmnet, 
                                      penalty_factor=penalty_factor_rest)
    #print('restrict')
    #print(restricted_soln)
  } else {
    restricted_soln = solve_restricted_problem(X, 
                                               y, 
                                               var, 
                                               lambda_glmnet, 
                                               penalty_factor=penalty_factor, 
                                               loss=loss)
  }

  n = nrow(X)
  idx = match(var, active_set)                    # active_set[idx]=var
  nuisance_res = (Qbeta_bar[var] -                # nuisance stat restricted to active vars
                  solve(QiE) %*% target_stat)/n 
  center = nuisance_res - (QE[idx,] %*% restricted_soln/n)
  radius = penalty_factor[var]*lambda_glmnet
  return(list(center=center*n, radius=radius*n))
}


linear_contrast = function(i, 
                           target_stat, 
                           target_cov, 
                           sigma_est,
                           center, 
                           radius) {

  target_cov_inv = solve(target_cov)
  target_stat_scaled = target_cov_inv %*% target_stat
  I=diag(length(target_stat)) * sigma_est^2
  nuisance = target_stat_scaled - I[,i]*target_stat[i]/target_cov[i,i]
  new_center = center[i]+ nuisance[i]
  new_radius = sqrt(max(radius^2 - sum((nuisance[-i])^2),0))
  return(list(target_cov=target_cov[i,i], center=new_center, radius=new_radius))

}

# returns P(Z > z | z in union of intervals)
tnorm.union.surv = function(z, 
                            mean, 
                            sd, 
                            intervals, 
                            bits=NULL){
  # intervals is a I x 2 matrix of disjoint intervals where the first column contains the lower endpoint
  
  pval = matrix(NA, nrow = dim(intervals)[1], ncol = length(mean))

  for(jj in 1:dim(intervals)[1]){
    if(z <= intervals[jj,1]){
      pval[jj,] = 1
    } else if(z >= intervals[jj,2]){
      pval[jj,] = 0
    }else{
      pval[jj,] = tnorm.surv(z, mean, sd, intervals[jj,1], intervals[jj,2], bits=bits)
    }
  }
  
  ww = matrix(NA, nrow=dim(intervals)[1], ncol=length(as.vector(mean)))

  for(jj in 1:dim(intervals)[1]){
    ww[jj,] = (pnorm(intervals[jj,2], mean=mean, sd=sd) - 
               pnorm(intervals[jj,1], mean=mean, sd=sd))
  }
  
  ww = ww%*%diag(as.vector(1/apply(ww,2,sum)), nrow=ncol(ww))
  
  pval = apply(pval*ww,2,sum)
  
  return(as.numeric(pval))
}

create_tnorm_interval = function(z, 
                                 sd, 
                                 alpha, 
                                 intervals, 
                                 gridrange=c(-20,20), 
                                 gridpts = 10000, 
                                 griddepth = 2, 
                                 bits = NULL){
  
  grid = seq(gridrange[1]*sd,gridrange[2]*sd,length=gridpts)
  fun = function(x) { 
        pv = tnorm.union.surv(z, x, sd, intervals, bits)
        return(pv)
  }
  
  conf_interval = grid.search(grid, fun, alpha/2, 1-alpha/2, gridpts, griddepth)
  
  return(conf_interval)
}


# GLM based functions

gradient = function(X, 
                    y, 
                    beta,
                    loss) {
  fit = X %*% beta
  if (loss=="logit"){
    fit = exp(fit) / (1 + exp(fit))
  }
  return(-t(X)%*%(y-fit))
}

hessian = function(X, 
                   beta, 
                   loss) {
  if (loss=="logit"){
    fit = X%*%beta
    W=diag(as.vector(exp(fit)/((1+exp(fit))^2)))
  } else if (loss=="ls"){
    W=diag(nrow(X))
  }
  return(t(X) %*% W %*% X)
}

hessian_active = function(X, 
                          beta, 
                          loss, 
                          active_set) {
  if (loss=="logit"){
    fit = X%*%beta
    W=diag(as.vector(exp(fit)/((1+exp(fit))^2)))
  } else if (loss=="ls"){
    W=diag(nrow(X))
  }
  return(t(X[, active_set]) %*% W %*% X)
}


mle = function(X,y,loss){
  reg = glm(y~X-1, family=family_label(loss))
  return(reg$coefficients)
}

# Functions to compute different debiasing matrices

approximate_JM = function(X, active_set){
  n=nrow(X)
  p=ncol(X)
  nactive=length(active_set)
  inactive_set=setdiff(1:p, active_set)
 
  is_wide = n < (2 * p)

  if (!is_wide) {
    hsigma = 1/n*(t(X) %*% X)
    htheta = debiasingMatrix(hsigma, is_wide, n, active_set)
  } else {
    htheta = debiasingMatrix(X, is_wide, n, active_set)
  }
  M_active <- htheta / n # the n is so that we approximate 
                         # the inverse of (X^TX)^{-1} instead of (X^TX/n)^{-1}.
  return(M_active)
}

approximate_BN = function(X, active_set){
	       # inspired by (6) of https://arxiv.org/pdf/1703.03282.pdf
               # approximate inverse is a scaled pseudo-inverse
	       # in the paper above, they use certain rows scaled X^{\dagger}
               # here we compute some scaled rows of (X^TX)^{\dagger} in such a way that
	       # our final answer, when multiplied by X^T agrees with (6).
  n=nrow(X)
  p=ncol(X)
  nactive=length(active_set)
  
  svdX = svd(X)
  
  approx_rank = sum(svdX$d > max(svdX$d) * 1.e-9)
  inv_d = 1. / svdX$d
  if (approx_rank < length(svdX$d)) { 
     inv_d[(approx_rank+1):length(svdX$d)] = 0.
  }
  inv = svdX$u[,1:approx_rank,drop=FALSE] %*% diag(inv_d[1:approx_rank]^2) %*% t(svdX$u)[1:approx_rank,,drop=FALSE]

  D = rep(0, nactive)

  for (i in 1:nactive){
    var = active_set[i]
    D[i] = 1/(t(X[,var]) %*% inv %*% X[,var])
  }

  pseudo_XTX = svdX$v[active_set,1:approx_rank,drop=FALSE] %*% diag(inv_d[1:approx_rank]^2) %*% t(svdX$v)[1:approx_rank,,drop=FALSE]

  M_active = diag(D) %*% pseudo_XTX # last two terms: projection onto row(X)
  return(M_active)
}


setup_Qbeta = function(X, 
                       y, 
                       soln, 
                       active_set, 
                       loss, 
                       debiasing_method,
                       use_debiased=FALSE){
  n=nrow(X)
  p=ncol(X)
  
  fit = X%*%soln
  if (loss=="ls"){
    W = rep(1,n)
  }
  if (loss=="logit"){
    W = exp(fit/2)/(1+exp(fit))  ## sqrt(pi*(1-pi))
  } 
  W_root=diag(as.vector(W))
  
  if ((n>p) & (!use_debiased)) {

    Q = hessian(X, soln, loss=loss)  ## X^TWX
    QE = Q[active_set,]
    Qi = solve(Q)   ## (X^TWX)^{-1}
    QiE = Qi[active_set,][, active_set]

    Xdesign = W_root %*% X
    beta_bar = soln - Qi %*% gradient(X, y, soln, loss=loss)
    Qbeta_bar = Q%*%soln - gradient(X, y, soln, loss=loss)
    beta_barE = beta_bar[active_set]
    M_active = Qi[active_set,]
    
  } else {
    
    if (debiasing_method == "JM") {
      ## this should be the active rows of \hat{\Sigma}(W^{1/2} X)^T/n, so size |E|\times p 
      M_active = approximate_JM(W_root %*% X, active_set) 
    }  else if (debiasing_method == "BN") {
      M_active = approximate_BN(W_root %*% X, active_set)
    }
      
    G = gradient(X, y, soln, loss=loss)
    beta_barE = soln[active_set] - M_active %*% G    

    M2 = M_active %*% t(W_root %*% X)
    QiE = M2 %*% t(M2) # size |E|\times |E|
    QE = hessian_active(X, soln, loss, active_set)
    Xdesign = W_root %*% X
    Qbeta_bar = t(QE)%*%soln[active_set] - G
  }
  
  return(list(QE=QE, 
              Xdesign=Xdesign, 
              Qbeta_bar=Qbeta_bar, 
              QiE=QiE, 
              beta_barE=beta_barE,
              M_active=M_active))  
}

ROSI = function(X, 
                y, 
                soln, 
                lambda, 
                penalty_factor=NULL, 
                dispersion=1,
                family=c('gaussian', 'binomial'),
                solver=c('QP', 'glmnet'),
                construct_ci=TRUE, 
                debiasing_method=c("JM", "BN"),
                verbose=FALSE,
		level=0.9,
                use_debiased=TRUE) {
  
  this.call = match.call()

  family = match.arg(family)
  solver = match.arg(solver)
  debiasing_method = match.arg(debiasing_method)

  active_set = which(soln!=0)
  nactive_set = length(active_set)
  
  if (is.null(penalty_factor)) { 
     penalty_factor = rep(1, ncol(X))
  }

  if (verbose){
    print(c("nactive", nactive_set))
  }
  
  begin_setup = Sys.time()
  setup_params = setup_Qbeta(X=X, 
                             y=y, 
                             soln=soln, 
                             active_set=active_set, 
                             loss=loss_label(family), 
                             debiasing_method=debiasing_method,
                             use_debiased=use_debiased)

  QE = as.matrix(setup_params$QE)
  Xdesign = setup_params$Xdesign
  QiE = as.matrix(setup_params$QiE)
  beta_barE = setup_params$beta_barE
  Qbeta_bar = setup_params$Qbeta_bar
  end_setup = Sys.time()
  
  if (verbose) {
    cat("setup time", end_setup-begin_setup, "\n")
  }
  
  pvalues = NULL
  selected_vars = NULL
  sel_intervals = NULL
  lower_trunc = NULL
  upper_trunc = NULL

  sigma_est = sqrt(dispersion)

  for (i in 1:nactive_set){
    
    target_stat = beta_barE[i]
    target_cov = as.matrix(QiE)[i,i]
    
    begin_TS = Sys.time()

    n = nrow(X)
    lambda_glmnet = lambda / n
    TS = truncation_set(X=X, 
                        y=y, 
                        Qbeta_bar=Qbeta_bar, 
                        QE=QE, 
                        Xdesign=Xdesign, 
                        QiE=QiE[i,i,drop=FALSE], # this is part of inverse Hessian, i.e. without dispersion
                        target_stat=target_stat, 
                        var=active_set[i], 
                        active_set=active_set,
                        lambda_glmnet=lambda_glmnet, 
                        penalty_factor=penalty_factor, 
                        loss=loss_label(family), 
                        solver=solver)
    end_TS = Sys.time()

    if (verbose) {
      print(c("current var", active_set[i]))
      cat("TS time", end_TS-begin_TS, "\n")
      cat("variable", active_set[i], "\n")
    }
    
    center = TS$center
    radius = TS$radius
    
    lower = target_cov * (-center - radius)
    upper = target_cov * (-center + radius)
    lower_trunc = c(lower_trunc, lower)
    upper_trunc = c(upper_trunc, upper)  

    if (target_stat<=lower | target_stat>=upper)
      {
        intervals = matrix(nrow=2, ncol=2)
        intervals[1,] = c(-Inf, lower)
        intervals[2,] = c(upper, Inf)
        pval = tnorm.union.surv(target_stat, 
                                mean=0, 
                                sd=sqrt(target_cov) * sigma_est, 
                                intervals)
        pval = 2*min(pval, 1-pval)

        if (verbose==TRUE){
          print(c("target stat", target_stat))
          print(c("target cov", target_cov))
          print(c("pval", pval))
        }

        pvalues = c(pvalues, pval)
       
        if (construct_ci) {
          sel_int = create_tnorm_interval(z=target_stat, 
                                          sd=sqrt(target_cov) * sigma_est, 
                                          alpha=1-level, 
                                          intervals=intervals)
          if (verbose==TRUE){
            cat("sel interval", sel_int, "\n")
          }
          sel_intervals = rbind(sel_intervals, sel_int)
        }
      } else {
           pvalues = c(pvalues, NA)
	   sel_intervals = rbind(sel_intervals, c(NA, NA))
           warning("observation not within the truncation limits!")
           print("observation not within the truncation limits!")
      }
  }
  
  out = list(pvalues=pvalues, 
             active_set=active_set,
             intervals=sel_intervals,
             estimate=as.numeric(beta_barE),
             std_err=sqrt(diag(QiE) * dispersion),
             dispersion=dispersion,
             lower_trunc=as.numeric(lower_trunc),
             upper_trunc=as.numeric(upper_trunc),
             lambda=lambda,
             penalty_factor=penalty_factor,
             level=level,
             call=this.call)
   class(out) = "ROSI"
   return(out)
}

print.ROSI <- function(x, ...) {
  cat("\nCall:\n")
  dput(x$call)

  cat(sprintf("\nDispersion taken to be dispersion = %0.3f\n",
              x$dispersion))

  cat(sprintf("\nTesting results at lambda = %0.3f, with level = %0.2f\n",x$lambda,x$level))
  cat("",fill=T)
  tab = cbind(signif(x$std_err^2,3),
              round(x$estimate,3),
              round(x$estimate / x$std_err,3),
              round(x$pvalues,3),
              round(x$intervals,3))
  colnames(tab) = c("Var", "Coef", "Z-score", "P-value", "LowConfPt", "UpConfPt")
  rownames(tab) = rep("",nrow(tab))
  print(tab)

  cat("\nNote: coefficients shown are full regression coefficients.\n")
  invisible()
}


# Some little used functions -- not exported

compute_coverage = function(ci, beta){

  nactive=length(beta)
  coverage_vector = rep(0, nactive)
  for (i in 1:nactive){
  if(!is.na(ci[i,1]) & !is.na(ci[i,2])) {
    if (beta[i]>=ci[i,1] && beta[i]<=ci[i,2]){
      coverage_vector[i]=1
    } 
  } else {
    coverage_vector[i] = NA
    }
  }
  return(coverage_vector)
}

loss_label = function(family) {
  if (family=="gaussian"){
    return("ls")
  } else if (family=="binomial"){
    return("logit")
  }
}

family_label = function(loss){
  if (loss=="ls"){
    return("gaussian")
  } else if (loss=="logit"){
    return("binomial")
  }
}


