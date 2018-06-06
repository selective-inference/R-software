family_label = function(loss){
  if (loss=="ls"){
    return("gaussian")
  } else if (loss=="logit"){
    return("binomial")
  }
}

# solves full lasso problem via glmnet

solve_problem_glmnet = function(X, y, lambda, penalty_factor, loss){
  if (is.null(lambda)){
    cv = cv.glmnet(x=X, 
                   y=y, 
                   family=family_label(loss), 
                   penalty.factor=penalty_factor, 
                   intercept=FALSE,  
                   thresh = 1e-12)
    beta_hat = coef(cv, s="lambda.min")
  }
  else {
    lasso = glmnet(x=X, 
                   y=y, 
                   family=family_label(loss), 
                   penalty.factor=penalty_factor,
                   alpha=1,
                   standardize=FALSE, 
                   intercept=FALSE, 
                   thresh=1e-20)
    beta_hat = coef(lasso, s=lambda)
  }
  return(beta_hat[-1])
}

# solves full group lasso problem via gglasso
solve_problem_gglasso = function(X, y, groups, lambda, penalty_factor, loss){
  if (is.null(lambda)){
    cv <- cv.gglasso(x=X, 
                     y=y, 
                     group=groups, 
                     loss=loss, 
                     pf=penalty_factor, 
                     intercept=FALSE, 
                     eps=1e-12)
    beta_hat = coef(cv, s="lambda.min")
  }
  else {
      # gglasso for logit loss needs the response to be in {-1,1}
      if (loss=="logit"){
        y_pm1 = rep(y)
        y_pm1[which(y==0)]=-1
      } else if (loss=="ls"){
       y_pm1 = rep(y)
      }
      m <- gglasso(x=X, 
                   y=y_pm1, 
                   group=groups, 
                   loss=loss, 
                   pf=penalty_factor, 
                   intercept=FALSE, 
                   eps=1e-20)
      beta_hat = coef(m, s=lambda)
  }
  return(beta_hat[-1])
}

# solves the restricted problem
solve_restricted_problem = function(X, y, var, lambda, penalty_factor, loss, algo){
  if (algo=="glmnet"){
    restricted_soln=rep(0, ncol(X))
    restricted_soln[-var] = solve_problem_glmnet(X[,-var], 
                                                 y, 
                                                 lambda, 
                                                 penalty_factor[-var], 
                                                 loss=loss)
  } else if (algo=="gglasso"){
    penalty_factor_rest = rep(penalty_factor)
    penalty_factor_rest[var] = 10^10
    restricted_soln = solve_problem_gglasso(X,
                                            y, 
                                            1:ncol(X), 
                                            lambda, 
                                            penalty_factor=penalty_factor_rest, 
                                            loss=loss)
  }
  return(restricted_soln)
}

solve_problem_Q = function(Q_sq, Qbeta_bar, lambda, penalty_factor,
                           max_iter=50,
                           kkt_tol=1.e-4, 
                           objective_tol=1.e-4, 
                           parameter_tol=1.e-4,
                           kkt_stop=TRUE,
                           objective_stop=TRUE,	
                           parameter_stop=TRUE){
  n=nrow(Q_sq)
  p=ncol(Q_sq)

  Xinfo = Q_sq
  linear_func = -as.numeric(Qbeta_bar)
  soln = as.numeric(rep(0., p))
  ever_active = as.integer(rep(0, p))
  ever_active[1] = 1
  ever_active = as.integer(ever_active)
  nactive = as.integer(1)
  Xsoln = as.numeric(rep(0, nrow(Xinfo)))
  gradient = 1. * linear_func 
  max_active=as.integer(p)
  
  linear_func = linear_func/n
  gradient = gradient/n
  
  #solve_QP_wide solves n*slinear_func^T\beta+\beta^T Xinfo\beta+\sum\lambda_i|\beta_i|
  result = solve_QP_wide(Xinfo,         # this is a design matrix
                         as.numeric(penalty_factor*lambda),  # vector of Lagrange multipliers
                         0,                          # ridge_term 
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
                          Q_sq, 
                          target_stat, 
                          target_cov,
                          var, 
                          active_vars,
                          lambda, 
                          penalty_factor, 
                          loss, 
                          algo){
  
  if (algo=="Q"){
    penalty_factor_rest = rep(penalty_factor)
    penalty_factor_rest[var] = 10^10
    restricted_soln = solve_problem_Q(Q_sq, 
                                      Qbeta_bar, 
                                      lambda, 
                                      penalty_factor=penalty_factor_rest)
  } else {
    restricted_soln = solve_restricted_problem(X, 
                                               y, 
                                               var, 
                                               lambda, 
                                               penalty_factor=penalty_factor, 
                                               loss=loss, 
                                               algo=algo)
  }

  n = nrow(X)
  idx = match(var, active_vars)                    # active_vars[idx]=var
  nuisance_res = (Qbeta_bar[var] -                 # nuisance stat restricted to active vars
                  solve(target_cov) %*% target_stat)/n 
  center = nuisance_res - (QE[idx,] %*% restricted_soln/n)
  radius = penalty_factor[var]*lambda
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

# the pvalue if prob(Z>obs given  |sigma_est^2/variance * Z+center|>radius)
# where Z~N(param, variance)
test_TG = function(param, 
                   observed, 
                   variance, 
                   sigma_est, 
                   center, 
                   radius, 
                   alt) {
  st.error = sqrt(variance)
  lower = variance*(-center-radius)/(sigma_est^2)
  upper = variance*(-center+radius)/(sigma_est^2)
  if (observed<=lower){
     case=1
     num = (pnorm(upper, 
                  mean=param, 
                  sd=st.error, 
                  lower.tail=FALSE) + 
            pnorm(lower, 
                  mean=param, 
                  sd=st.error) - 
            pnorm(observed, mean=param, sd=st.error))
  } else if (observed>=upper){
    case=2
    num = pnorm(observed, 
                mean=param, 
                sd=st.error, 
                lower.tail=FALSE)
  } else{
    case=3
    num = pnorm(upper, 
                mean=param, 
                sd=st.error, 
                lower.tail=FALSE)
    return(NULL)
  }
  den = (pnorm(upper, 
               mean=param, 
               sd=st.error, 
               lower.tail=FALSE) + 
         pnorm(lower, 
               mean=param, 
               sd=st.error))
         pivot = num/den

  if (alt=="two-sided"){
    return(2*pmin(pivot, 1-pivot))
  } else if (alt=="upper") {
    return(pivot)
  } else if (alt=="lower"){
    return(1-pivot)
  }
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
    }else if(z >= intervals[jj,2]){
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

selective_CI = function(observed, 
                        variance, 
                        sigma_est, 
                        center, 
                        radius,
                        alpha=0.1, 
                        gridrange=c(-20,20), 
                        gridpts=10000, 
                        griddepth=2) {
  
  pivot = function(param){
    return(test_TG(param, 
                   observed, 
                   variance, 
                   sigma_est, 
                   center, 
                   radius, 
                   alt="upper"))
  }
  st.error = sqrt(variance)
  param_grid = seq(observed+gridrange[1]*st.error, 
                   observed+gridrange[2]*st.error, 
                   length=gridpts)
  interval = grid.search(param_grid, 
                         pivot, 
                         alpha/2, 
                         1-alpha/2, 
                         gridpts, 
                         griddepth)
  return(interval)
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
  n=nrow(X)
  p=ncol(X)
  nactive=length(active_set)
  
  svdX = svd(X)
  inv = solve(svdX$v %*% diag(svdX$d^2) %*% t(svdX$v))
  D = matrix(rep(0, nactive*p), nrow=nactive, ncol=p)
  for (i in 1:nactive){
    var = active_set[i]
    D[i, var] = 1/(t(X[,var]) %*% inv %*% X[,var])
  }
  M_active = D %*% svdX$v %*% t(svdX$v) # last two terms: projection onto row(X)
  return(M_active)
}


setup_Qbeta = function(X, 
                       y, 
                       soln, 
                       active_set, 
                       loss, 
                       debias_mat,
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

    Q_sq = W_root %*% X
    beta_bar = soln - Qi %*% gradient(X, y, soln, loss=loss)
    Qbeta_bar = Q%*%soln - gradient(X, y, soln, loss=loss)
    beta_barE = beta_bar[active_set]
    M_active = Qi[active_set,]
    
  } else {
    
    if (debias_mat == "JM"){
      ## this should be the active rows of \hat{\Sigma}(W^{1/2} X)^T/n, so size |E|\times p 
      M_active = approximate_JM(W_root %*% X, active_set) 
    }  else if (debias_mat == "BN"){
      M_active = approximate_BN(W_root %*% X, active_set)
    }
      
    G = gradient(X, y, soln, loss=loss)
    beta_barE = soln[active_set] - M_active %*% G    

    M2 = M_active %*% t(W_root %*% X)
    QiE = M2 %*% t(M2) # size |E|\times |E|
    QE = hessian_active(X, soln, loss, active_set)
    Q_sq = W_root %*% X
    Qbeta_bar = t(QE)%*%soln[active_set] - G
  }
  
  return(list(QE=QE, 
              Q_sq=Q_sq, 
              Qbeta_bar=Qbeta_bar, 
              QiE=QiE, 
              beta_barE=beta_barE,
              M_active=M_active))  
}

inference_debiased_full = function(X, 
                                   y, 
                                   soln, 
                                   lambda, 
                                   penalty_factor, 
                                   sigma_est,
                                   loss, 
                                   algo, 
                                   construct_ci, 
                                   debias_mat="JM", 
                                   verbose=FALSE,
                                   use_debiased=FALSE){
  
  active_vars = which(soln!=0)
  nactive_vars = length(active_vars)
  
  if (nactive_vars==0){
    return(list(pvalues=NULL, naive_pvalues=NULL))
  }
  
  if (verbose){
    print(c("nactive", nactive_vars))
  }
  
  begin_setup = Sys.time()
  setup_params = setup_Qbeta(X=X, 
                             y=y, 
                             soln=soln, 
                             active_set=active_vars, 
                             loss=loss, 
                             debias_mat=debias_mat,
                             use_debiased=use_debiased)

  QE = as.matrix(setup_params$QE)
  Q_sq = setup_params$Q_sq
  QiE = as.matrix(setup_params$QiE)
  beta_barE = setup_params$beta_barE
  Qbeta_bar = setup_params$Qbeta_bar
  end_setup = Sys.time()
  
  if (verbose) {
    cat("setup time", end_setup-begin_setup, "\n")
  }
  
  pvalues = NULL
  naive_pvalues = NULL
  selected_vars = NULL
  sel_intervals = NULL
  naive_intervals= NULL
  
  for (i in 1:nactive_vars){
    
    target_stat = beta_barE[i]
    target_cov = as.matrix(QiE)[i,i]
    
    begin_TS = Sys.time()
    # JT would be good 
    TS =  truncation_set(X=X, 
                         y=y, 
                         Qbeta_bar=Qbeta_bar, 
                         QE=QE, 
                         Q_sq=Q_sq, 
                         target_cov=target_cov, # this Hessian, i.e. without dispersion
                         target_stat=target_stat, 
                         var=active_vars[i], 
                         active_vars=active_vars,
                         lambda=lambda, 
                         penalty_factor=penalty_factor, 
                         loss=loss, 
                         algo=algo)
    end_TS = Sys.time()

    if (verbose) {
      print(c("current var", active_vars[i]))
      cat("TS time", end_TS-begin_TS, "\n")
      cat("variable", active_vars[i], "\n")
    }
    
    center = TS$center
    radius = TS$radius
    
    lower = target_cov * (-center - radius)
    upper = target_cov * (-center + radius)
      
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
        
        selected_vars = c(selected_vars, active_vars[i])
        
        # JT: why are we reporting naive ones -- which we know are bad?
        naive_pval =  pvalue_naive_linear(target_stat, target_cov * sigma_est^2)
        naive_pvalues = c(naive_pvalues, naive_pval)
        
        if (construct_ci) {
          sel_int = create_tnorm_interval(z=target_stat, 
                                          sd=sqrt(target_cov) * sigma_est, 
                                          alpha=0.1, 
                                          intervals=intervals)
          naive_int = naive_CI(target_stat, target_cov * sigma_est^2)
          if (verbose==TRUE){
            cat("sel interval", sel_int, "\n")
            cat("naive interval", naive_int, "\n")
          }
          sel_intervals = cbind(sel_intervals, sel_int)
          naive_intervals = cbind(naive_intervals, naive_int)
        }
      } else {
           warning("observation not within the truncation limits!")
      }
  }
  
  return(list(pvalues=pvalues, 
              naive_pvalues=naive_pvalues, 
              active_vars=selected_vars,
              sel_intervals=sel_intervals, 
              naive_intervals=naive_intervals,
              M_active=setup_params$M_active))
}

# Some little used functions -- not exported

pvalue_naive_linear = function(observed, variance){
  pval = pnorm(observed, 
               mean=0, 
               sd=sqrt(variance), 
               lower.tail=FALSE)
  return(2*min(pval, 1-pval))
}

naive_CI = function(observed, variance, alpha=0.1){
  quantile = abs(qnorm(alpha/2, mean=0, sd=1))
  st.error = sqrt(variance)
  interval = c(observed - st.error*quantile, observed + st.error*quantile)
  return(interval)
}


compute_coverage = function(ci, beta){
  nactive=length(beta)
  coverage_vector = rep(0, nactive)
  for (i in 1:nactive){
    if (beta[i]>=ci[1,i] && beta[i]<=ci[2,i]){
      coverage_vector[i]=1
    } 
  }
  return(coverage_vector)
}
