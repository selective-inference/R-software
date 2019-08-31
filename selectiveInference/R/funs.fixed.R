# Lasso inference function (for fixed lambda). Note: here we are providing inference
# for the solution of
# min 1/2 || y - \beta_0 - X \beta ||_2^2 + \lambda || \beta ||_1

fixedLassoInf <- function(x, y, beta, 
                          lambda, 
                          family=c("gaussian","binomial","cox"),
                          intercept=TRUE, 
                          add.targets=NULL, 
                          status=NULL,
                          sigma=NULL, 
                          alpha=0.1,
                          type=c("partial", "full"), 
                          tol.beta=1e-5, 
                          tol.kkt=0.1,
                          gridrange=c(-100,100), 
			  bits=NULL, 
			  verbose=FALSE, 
                          linesearch.try=10) {

  family = match.arg(family)
  this.call = match.call()
  type = match.arg(type)
  if(family=="binomial")  {
    if(type!="partial") stop("Only type= partial allowed with binomial family")
    out=fixedLogitLassoInf(x,y,beta,lambda,alpha=alpha, type=type, tol.beta=tol.beta, tol.kkt=tol.kkt,
                           gridrange=gridrange, bits=bits, verbose=verbose,
                           linesearch.try=linesearch.try,
                           this.call=this.call)
    return(out)
  }
  else if(family=="cox")  {
    if(type!="partial") stop("Only type= partial allowed with Cox family")
    out=fixedCoxLassoInf(x,y,status,beta,lambda,alpha=alpha, type="partial",tol.beta=tol.beta,
                         tol.kkt=tol.kkt, gridrange=gridrange, bits=bits, verbose=verbose,this.call=this.call)
    return(out)
  }
  
  else{
    
    checkargs.xy(x,y)
    if (missing(beta) || is.null(beta)) stop("Must supply the solution beta")
    if (missing(lambda) || is.null(lambda)) stop("Must supply the tuning parameter value lambda")
    
    n = nrow(x)
    p = ncol(x)
    beta = as.numeric(beta)
    if (type == "full") {
      if (p > n) {
        # need intercept (if there is one) for debiased lasso
        hbeta = beta
        if (intercept == T) {
          if (length(beta) != p + 1) {
            stop("Since type='full', p > n, and intercept=TRUE, beta must have length equal to ncol(x)+1")
          }
          # remove intercept if included
          beta = beta[-1]
        } else if (length(beta) != p) {
          stop("Since family='gaussian', type='full' and intercept=FALSE, beta must have length equal to ncol(x)")
        }
      }
    } else if (length(beta) != p) {
      stop("Since family='gaussian' and type='partial', beta must have length equal to ncol(x)")
    }
    
    checkargs.misc(beta=beta,lambda=lambda,sigma=sigma,alpha=alpha,
                   gridrange=gridrange,tol.beta=tol.beta,tol.kkt=tol.kkt)
    if (!is.null(bits) && !requireNamespace("Rmpfr",quietly=TRUE)) {
      warning("Package Rmpfr is not installed, reverting to standard precision")
      bits = NULL
    }
    
    if (!is.null(add.targets) && (!is.vector(add.targets)
                                  || !all(is.numeric(add.targets)) || !all(add.targets==floor(add.targets))
                                  || !all(add.targets >= 1 && add.targets <= p))) {
      stop("'add.targets' must be a vector of integers between 1 and p")
    }
    
    # If glmnet was run with an intercept term, center x and y
    if (intercept==TRUE) {
      obj = standardize(x,y,TRUE,FALSE)
      x = obj$x
      y = obj$y
    }
    
    # Check the KKT conditions
    g = t(x)%*%(y-x%*%beta) / lambda
    if (any(abs(g) > 1+tol.kkt * sqrt(sum(y^2))))
      warning(paste("Solution beta does not satisfy the KKT conditions",
                    "(to within specified tolerances)"))
    
    tol.coef = tol.beta * sqrt(n / colSums(x^2))
    vars = which(abs(beta) > tol.coef)
    sign_vars = sign(beta[vars])

    if(sum(vars)==0){
      cat("Empty model",fill=T)
      return()
    }

    if (any(sign(g[vars]) != sign_vars)) {
      warning(paste("Solution beta does not satisfy the KKT conditions",
                    "(to within specified tolerances). You might try rerunning",
                    "glmnet with a lower setting of the",
                    "'thresh' parameter, for a more accurate convergence."))
    }

    # Get lasso polyhedral region, of form Gy >= u

    logical.vars=rep(FALSE,p)
    logical.vars[vars]=TRUE
    
    if (type == 'full') {
       out = fixedLassoPoly(x, y, lambda, beta, logical.vars, inactive=TRUE)
    } 
    else {
       out = fixedLassoPoly(x, y, lambda, beta, logical.vars)
    }
    
    A = out$A
    b = out$b
    
    # Check polyhedral region
    tol.poly = 0.01
    if (max(A %*% y - b) > tol.poly * sqrt(sum(y^2)))
      stop(paste("Polyhedral constraints not satisfied; you must recompute beta",
                 "more accurately. With glmnet, make sure to use exact=TRUE in coef(),",
                 "and check whether the specified value of lambda is too small",
                 "(beyond the grid of values visited by glmnet).",
                 "You might also try rerunning glmnet with a lower setting of the",
                 "'thresh' parameter, for a more accurate convergence."))
    
    # Estimate sigma
    if (is.null(sigma)) {
      if (n >= 2*p) {
        oo = intercept
        sigma = sqrt(sum(lsfit(x,y,intercept=oo)$res^2)/(n-p-oo))
      }
      else {
        sigma = sd(y)
        warning(paste(sprintf("p > n/2, and sd(y) = %0.3f used as an estimate of sigma;",sigma),
                      "you may want to use the estimateSigma function"))
      }
    }
    
    # add additional targets for inference if provided
    if (!is.null(add.targets)) {
       # vars is boolean...
       old_vars = vars & TRUE
       vars[add.targets] = TRUE
       sign_vars = sign(beta[vars]) 
       sign_vars[!old_vars] = NA
       stop("`add.targets` not fully implemented yet")
    }

    k = length(vars)
    pv = vlo = vup = numeric(k)
    vmat = matrix(0,k,n)
    ci = tailarea = matrix(0,k,2)
      
    if (type=="full" & p > n) {
      if (intercept == TRUE) {
        pp=p+1
        Xint <- cbind(rep(1,n),x)
        # indices of selected predictors
        S = c(1,vars + 1)
      } else {
        pp=p
        Xint <- x
        # indices of selected predictors
        S = vars
        # notS = which(abs(beta) <= tol.coef)
      }
      
      notS = setdiff(1:pp,S)
      
      XS = Xint[,S]
      hbetaS = hbeta[S]
      
      # Reorder so that active set S is first
      Xordered = Xint[,c(S,notS,recursive=T)]
      hsigmaS = 1/n*(t(XS)%*%XS) # hsigma[S,S]
      hsigmaSinv = solve(hsigmaS) # pinv(hsigmaS)
      
      FS = rbind(diag(length(S)),matrix(0,pp-length(S),length(S)))
      GS = cbind(diag(length(S)),matrix(0,length(S),pp-length(S)))

      is_wide = n < (2 * p) # somewhat arbitrary decision -- it is really for when we don't want to form with pxp matrices

      # Approximate inverse covariance matrix for when (n < p) from lasso_Inference.R
      if (!is_wide) {
           hsigma = 1/n*(t(Xordered)%*%Xordered)
           htheta = debiasingMatrix(hsigma, 
                                    is_wide, 
                                    n, 
                                    1:length(S), 
                                    verbose=FALSE, 
                                    max_try=linesearch.try, 
                                    warn_kkt=TRUE)
           ithetasigma = (GS-(htheta%*%hsigma))
      } else {
           htheta = debiasingMatrix(Xordered, 
                                    is_wide, 
                                    n, 
                                    1:length(S), 
                                    verbose=FALSE, 
                                    max_try=linesearch.try, 
                                    warn_kkt=TRUE)
           ithetasigma = (GS-((htheta%*%t(Xordered)) %*% Xordered)/n)
      }

      M <- (((htheta%*%t(Xordered))+ithetasigma%*%FS%*%hsigmaSinv%*%t(XS))/n)

      # vector which is offset for testing debiased beta's
      null_value <- (((ithetasigma%*%FS%*%hsigmaSinv)%*%sign(hbetaS))*lambda/n)

      if (intercept == T) {
        M = M[-1,] # remove intercept row
        null_value = null_value[-1] # remove intercept element
      }
    } else if (type=="partial" || p > n) {
      xa = x[,vars,drop=F]
      M = pinv(crossprod(xa)) %*% t(xa)
      null_value = rep(0,k)
    } else {
      M = pinv(crossprod(x)) %*% t(x)
      M = M[vars,,drop=F]
      null_value = rep(0,k)
    }

  for (j in 1:k) {
    if (verbose) cat(sprintf("Inference for variable %i ...\n",vars[j]))

    vj = M[j,]
    mj = sqrt(sum(vj^2))
    vj = vj / mj        # Standardize (divide by norm of vj)

    if (!is.na(sign_vars[j])) {
        vj = sign_vars[j] * vj
    }

    limits.info = TG.limits(y, A, b, vj, Sigma=diag(rep(sigma^2, n)))
    a = TG.pvalue.base(limits.info, null_value=null_value[j], bits=bits)
    pv[j] = a$pv
    
    if (is.na(sign_vars[j])) { # for variables not in the active set, report 2-sided pvalue
       pv[j] = 2 * min(pv[j], 1 - pv[j])
    }
    
    vlo[j] = a$vlo * mj # Unstandardize (mult by norm of vj)
    vup[j] = a$vup * mj # Unstandardize (mult by norm of vj)
    if (!is.na(sign_vars[j])) { 
        vmat[j,] = vj * mj * sign_vars[j]  # Unstandardize (mult by norm of vj) and fix sign
    } else {
        vmat[j,] = vj * mj # Unstandardize (mult by norm of vj)
    }
    a = TG.interval.base(limits.info, 
                         alpha=alpha,
                         gridrange=gridrange,
			 flip=(sign_vars[j]==-1),
                         bits=bits)
    ci[j,] = (a$int-null_value[j]) * mj # Unstandardize (mult by norm of vj)
    tailarea[j,] = a$tailarea
  }

  out = list(type=type,
             lambda=lambda,
             pv=pv,
             ci=ci,
             tailarea=tailarea,
             vlo=vlo,
             vup=vup,
             vmat=vmat,
             y=y,
             vars=vars,
             sign=sign_vars,
             sigma=sigma,
             alpha=alpha,
             sd=sigma*sqrt(rowSums(vmat^2)),
             coef0=vmat%*%y,
             call=this.call)

  class(out) = "fixedLassoInf"
  return(out)
}
}

#############################


fixedLassoPoly =
  function(X, y, lambda, beta, active, inactive = FALSE) {

    XA = X[, active, drop=FALSE]
    XI = X[, !active, drop=FALSE]
    XAi = pinv(crossprod(XA))
    XAp = XAi %*% t(XA)
    Ir = t(XI) %*% t(XAp)  # matrix in the "irrepresentable" condition

    if(length(lambda)>1) {
       lambdaA= lambda[active]
       lambdaI = lambda[!active]
    } else {
       lambdaA = rep(lambda, sum(active))
       lambdaI = rep(lambda, sum(!active))
    }

    penalized = lambdaA != 0
    signA = sign(beta[active])
    active_subgrad = signA * lambdaA
    if (length(signA)>1) sign_diag = diag(signA)
    if (length(signA)==1) sign_diag = matrix(signA, 1, 1)
    
    if (inactive) { # should we include the inactive constraints?
      RA = diag(rep(1, nrow(XA))) - XA %*% XAp # RA is residual forming matrix of selected model
      
      A = rbind(
        t(XI) %*% RA,
        -t(XI) %*% RA,
        -(sign_diag %*% XAp)[penalized,] # no constraints for unpenalized
      )

      b = c(
        lambdaI - Ir %*% active_subgrad,
        lambdaI + Ir %*% active_subgrad,
        -(sign_diag %*% XAi %*% active_subgrad)[penalized])
    } else {
      A = -(sign_diag %*% XAp)[penalized,]  # no constraints for unpenalized
      b = -(sign_diag %*% XAi %*% active_subgrad)[penalized]
    }
    
    return(list(A=A, b=b))
  }

##############################

## Approximates inverse covariance matrix theta
## using coordinate descent 

debiasingMatrix = function(Xinfo,               # could be X or t(X) %*% X / n 
                                                # depending on is_wide
                           is_wide,
                           nsample, 
                           rows, 
                           verbose=FALSE, 
                           bound=NULL,          # starting value of bound
                           linesearch=TRUE,     # do a linesearch?
                           scaling_factor=1.5,  # multiplicative factor for linesearch
                           max_active=NULL,     # how big can active set get?
                           max_try=10,          # how many steps in linesearch?
                           warn_kkt=FALSE,      # warn if KKT does not seem to be satisfied?
                           max_iter=50,         # how many iterations for each optimization problem
                           kkt_stop=TRUE,       # stop based on KKT conditions?
                           parameter_stop=TRUE, # stop based on relative convergence of parameter?
                           objective_stop=TRUE, # stop based on relative decrease in objective?
                           kkt_tol=1.e-4,       # tolerance for the KKT conditions
                           parameter_tol=1.e-4, # tolerance for relative convergence of parameter
                           objective_tol=1.e-4  # tolerance for relative decrease in objective
                           ) {


  if (is.null(max_active)) {
     max_active = max(50, 0.3 * nsample)
  } 

  p = ncol(Xinfo);
  M = matrix(0, length(rows), p);

  if (is.null(bound)) {
      bound = (1/sqrt(nsample)) * qnorm(1-(0.1/(p^2)))
  }
 
  xperc = 0;
  xp = round(p/10);
  idx = 1;
  for (row in rows) {
    if ((idx %% xp)==0){
      xperc = xperc+10;
      if (verbose) {
        print(paste(xperc,"% done",sep="")); }
    }

    output = debiasingRow(Xinfo,               # could be X or t(X) %*% X / n 
                                               # depending on is_wide
                          is_wide,
                          row,
                          bound,
                          linesearch=linesearch,
                          scaling_factor=scaling_factor,
                          max_active=max_active,
			  max_try=max_try,
			  warn_kkt=FALSE,
			  max_iter=max_iter,
			  kkt_stop=kkt_stop,
			  parameter_stop=parameter_stop,
			  objective_stop=objective_stop,
			  kkt_tol=kkt_tol,
			  parameter_tol=parameter_tol,
			  objective_tol=objective_tol)

    if (warn_kkt && (!output$kkt_check)) {
       warning("Solution for row of M does not seem to be feasible")
    } 
  
    if (!is.null(output$soln)) {
        M[idx,] = output$soln;
    } else {
        stop(paste("Unable to approximate inverse row ", row));
    }

    idx = idx + 1;
  }
  return(M)
}

# Find one row of the debiasing matrix -- assuming X^TX/n is not too large -- i.e. X is tall

debiasingRow = function (Xinfo,               # could be X or t(X) %*% X / n depending on is_wide
                         is_wide, 
                         row, 
                         bound, 
	                 linesearch=TRUE,     # do a linesearch?
		         scaling_factor=1.5,  # multiplicative factor for linesearch
		         max_active=NULL,     # how big can active set get?
			 max_try=10,          # how many steps in linesearch?
			 warn_kkt=FALSE,      # warn if KKT does not seem to be satisfied?
			 max_iter=50,         # how many iterations for each optimization problem
                         kkt_stop=TRUE,       # stop based on KKT conditions?
                         parameter_stop=TRUE, # stop based on relative convergence of parameter?
                         objective_stop=TRUE, # stop based on relative decrease in objective?
                         kkt_tol=1.e-4,       # tolerance for the KKT conditions
			 parameter_tol=1.e-4, # tolerance for relative convergence of parameter
			 objective_tol=1.e-4  # tolerance for relative decrease in objective
                         ) {

  p = ncol(Xinfo)

  if (is.null(max_active)) {
      max_active = min(nrow(Xinfo), ncol(Xinfo))
  }

   
  # Initialize variables 

  soln = rep(0, p)
  soln = as.numeric(soln)
  ever_active = rep(0, p)
  ever_active[1] = row      # 1-based
  ever_active = as.integer(ever_active)
  nactive = as.integer(1)

  linear_func = rep(0, p)
  linear_func[row] = -1
  linear_func = as.numeric(linear_func)
  gradient = 1. * linear_func 

  counter_idx = 1;
  incr = 0;

  last_output = NULL

  if (is_wide) {
     Xsoln = as.numeric(rep(0, nrow(Xinfo)))
  }

  while (counter_idx < max_try) {

      if (!is_wide) {
          result = solve_QP(Xinfo, # this is non-neg-def matrix
                            bound, 
                            max_iter, 
                            soln, 
                            linear_func, 
                            gradient, 
                            ever_active, 
                            nactive, 
                            kkt_tol, 
                            objective_tol, 
			    parameter_tol,
                            max_active,
			    kkt_stop,
			    objective_stop,
			    parameter_stop)
      } else {
          result = solve_QP_wide(Xinfo,                      # this is a design matrix
                                 as.numeric(rep(bound, p)),  # vector of Lagrange multipliers
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

      }

      iter = result$iter

      # Logic for whether we should continue the line search

      if (!linesearch) {
        break
      }

      if (counter_idx == 1){
        if (iter == (max_iter+1)){
           incr = 1; # was the original problem feasible? 1 if not
         } else {
           incr = 0; # original problem was feasible
         }
      } 

      if (incr == 1) { # trying to find a feasible point
         if ((iter < (max_iter+1)) && (counter_idx > 1)) { 
           break;      # we've found a feasible point and solved the problem            
         }
         bound = bound * scaling_factor;
      } else {         # trying to drop the bound parameter further
         if ((iter == (max_iter + 1)) && (counter_idx > 1)) {
            result = last_output; # problem seems infeasible because we didn't solve it
   	    break;                # so we revert to previously found solution
         }
         bound = bound / scaling_factor;
      }

      # If the active set has grown to a certain size
      # then we stop, presuming problem has become
      # infeasible.

      # We revert to the previous solution
	
      if (result$max_active_check) {
	  result = last_output;
	  break;
      }
      
      counter_idx = counter_idx + 1
      last_output = list(soln=result$soln,
                         kkt_check=result$kkt_check)
    }


  # Check feasibility

  if (warn_kkt && (!result$kkt_check)) {
     warning("Solution for row of M does not seem to be feasible")
  } 

  return(list(soln=result$soln,
              kkt_check=result$kkt_check,
	      gradient=result$gradient))

}


##############################

print.fixedLassoInf <- function(x, tailarea=TRUE, ...) {
  cat("\nCall:\n")
  dput(x$call)

  cat(sprintf("\nStandard deviation of noise (specified or estimated) sigma = %0.3f\n",
              x$sigma))

  cat(sprintf("\nTesting results at lambda = %0.3f, with alpha = %0.3f\n",x$lambda,x$alpha))
  cat("",fill=T)
  tab = cbind(x$vars,
    round(x$coef0,3),
    round(x$coef0 / x$sd,3),
    round(x$pv,3),round(x$ci,3))
  colnames(tab) = c("Var", "Coef", "Z-score", "P-value", "LowConfPt", "UpConfPt")
  if (tailarea) {
    tab = cbind(tab,round(x$tailarea,3))
    colnames(tab)[(ncol(tab)-1):ncol(tab)] = c("LowTailArea","UpTailArea")
  }
  rownames(tab) = rep("",nrow(tab))
  print(tab)

  cat(sprintf("\nNote: coefficients shown are %s regression coefficients\n",
              ifelse(x$type=="partial","partial","full")))
  invisible()
}

#estimateLambda <- function(x, sigma, nsamp=1000){
#  checkargs.xy(x,rep(0,nrow(x)))
#  if(nsamp < 10) stop("More Monte Carlo samples required for estimation")
#  if (length(sigma)!=1) stop("sigma should be a number > 0")
 # if (sigma<=0) stop("sigma should be a number > 0")

 # n = nrow(x)
 # eps = sigma*matrix(rnorm(nsamp*n),n,nsamp)
 # lambda = 2*mean(apply(t(x)%*%eps,2,max))
 # return(lambda)
#}
