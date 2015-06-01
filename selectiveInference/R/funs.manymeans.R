### functions for computing the many means estimates- Stephen Reid
### returns 
###   i) selected indices
###   ii) selection adjusted point estimates
###   iii) selection adjusted interval estimates
###   iv) selection adjusted p-value of hypothesis testing whether underlying signal is 0




#########################
##### MAIN FUNCTION #####
#########################

#### user-facing function for computing
#### selected set
#### point and interval estimates
#### p-values
#### input:
####    - y = Vector of observations
####    - alpha = Significance level used in CI construction
####    - bh.q = q parameter for BH(q) procedure (default: NULL)
####    - k = Number of largest elements to consider (default: NULL)
####    - sigma = Estimate of standard deviation of one of the components
#### output:
####    * A list (of class "mm") with the following components :
####        - muhat = Vector of length length(y) containing the estimated signal size. If a sample element is not selected, then its signal size estimate is 0
####        - selected.set = Indices into the vector y of the sample elements that were selected by our procedure (either BH(q) or top-K)
####        - CIs = Matrix with two columns and number of rows equal to number of elements in selected.set. Provides the post-selection CI bounds for the estimated signal sizes of selected elements. CIs given is rows in the same order as encountered in selected.set
####        - p.vals = Vector of p-values for the test of nullity of the signals of the selected sample elemetns. P-values given in the same order as selected.set

manyMeans = function(y, alpha=0.05, bh.q=NULL, k=NULL, sigma=1){
    this.call=match.call()
  #### error checks
    checkargs(y,alpha=alpha,bh.q=bh.q,k=k,sigma=sigma)
    if(!is.wholenumber(k) | k<1 | k>length(y)) stop("k must be an integer between 1 and length of y")
  if (is.null(bh.q) & is.null(k)){
    print ("Please set either the BH control parameter (bh.q) or the top-k parameter (k). They cannot both be NULL.")
    return ()
  }
  
  
  n = length(y)
  if (!is.null(bh.q)){ # use BH selection procedure
    if (bh.q > 1 | bh.q < 0){
      print ("Please specify a bh.q value between 0 and 1.")
      return ()
    }
    if (!is.null(k)){
      warning ("Both bh.q and k have been specified. k ignored.")
    }
    
    ### find the selected set and threshold
    p.vals = 2*pnorm (abs(y)/sigma, 0, 1, lower.tail=FALSE)
    order.p.vals = order(p.vals)
    sorted.p.vals = p.vals[order.p.vals]
    
    options (warn=-1) # ignore warning if max is over empty set
    last.reject = max(which (sorted.p.vals <= bh.q*(1:n)/n))
    options (warn=0) # reinstitute warnings
    
    if (last.reject == -Inf){ # none rejected
      print ("No sample elements selected.")
      
      out = list(muhat=rep(0, n), selected.set=NULL, CIs=NULL, p.vals=NULL, method="BH(q)", q=bh.q, k=NULL, threshold=NULL)
      class (out) = "manyMeans"
      return (out)
    }
    
    selected.set = order.p.vals[1:n <= last.reject]
    threshold = sigma*qnorm (last.reject/2/n, lower.tail=FALSE)
  }else{ # use top-k selection procedure
    ### find the selected set and threshold
    if (k == n){ # make no changes - return MLE
      z.alpha = qnorm (alpha/2, lower.tail=FALSE)
      cis = cbind(y - z.alpha*sigma, y + z.alpha*sigma)
      p.vals = 2*pnorm (abs(y), 0, sigma, lower.tail=FALSE)
      
      out = list(muhat=y, selected.set=1:n, CIs=cis, p.vals=round(p.vals, 4), method="top-K", q=NULL, k=k, threshold=NULL)
      class(out) = "manyMeans"
      return (out)
    }
    
    order.abs.y = order (-abs(y))
    sorted.abs.y = y[order.abs.y]
    
    selected.set = order.abs.y[1:k]
    threshold = abs(sorted.abs.y[k+1])
  }
  
  ### estimate their underlying signal sizes
  mu.hat = sapply (selected.set, function(s){
    uniroot (f=function(mu){tn.mean(mu, -threshold, threshold, sigma=sigma) - y[s]}, lower=-10000*sigma, upper=10000*sigma)$root
  })
  
  ### and CIs
  right.ci = sapply (selected.set, function(s){
    uniroot (f=function(mu){tn.cdf (y[s], mu, -threshold, threshold, sigma=sigma) - (alpha/2)}, lower=-10000*sigma, upper=10000*sigma)$root
  })
  left.ci = sapply (selected.set, function(s){
    uniroot (f=function(mu){tn.cdf (y[s], mu, -threshold, threshold, sigma=sigma) - (1-alpha/2)}, lower=-10000*sigma, upper=10000*sigma)$root
  })
  
  ### and p-values
  p.val = sapply (selected.set, function(s){tn.cdf (y[s], 0, -threshold, threshold, sigma=sigma)})
  p.val = 2*pmin(p.val, 1-p.val)
  
  ### arrange
  order.selected.set = order (selected.set)
  selected.set = selected.set[order.selected.set]
  mu.hat = mu.hat[order.selected.set]
  left.ci = left.ci[order.selected.set]
  right.ci = right.ci[order.selected.set]
  p.val = p.val[order.selected.set]
  
  mu.hat.final = rep(0, n)
  mu.hat.final[selected.set] = mu.hat
  
  out = list(muhat=mu.hat.final, selected.set=selected.set,selint=cbind(left.ci, right.ci), p.vals=round(p.val, 4), method=ifelse(is.null(bh.q), "top-K", "BH(q)"), q=bh.q, k=k, threshold=threshold, call=this.call)
  class (out) = "manyMeans"
  return (out)
}



###############################
##### AUXILIARY FUNCTIONS #####
###############################

#### function returning the cumulative distribution function value
#### of a truncated Gaussian RV, truncated to interval (-Inf, a) \union (b, Inf)
#### with underlying Gaussian having mean parameter mu and standard deviation sigma
#### at value y
tn.cdf = function(y, mu, a, b, sigma=1){
  ## denominator
  d_right = pnorm (b, mu, sigma, lower.tail=FALSE, log.p=TRUE)
  d_left = pnorm (a, mu, sigma, lower.tail=TRUE, log.p=TRUE)
  d_max = max(d_right, d_left)
  d_log = d_max + log(exp(d_left - d_max) + exp(d_right - d_max))
  
  
  # numerator
  if (y > a & y < b){
    n_log = d_left
    return (exp(n_log-d_log))
  }else{
    if (y > b){
      # b and y
      n_y_tilde = pnorm (y, mu, sigma, lower.tail=FALSE, log.p=TRUE)
      n_b_tilde = pnorm (b, mu, sigma, lower.tail=FALSE, log.p=TRUE)
      n_yb = n_b_tilde + log(1 - exp(n_y_tilde-n_b_tilde))
      
      # a
      n_a = d_left
      
      # combine
      return(exp(n_yb-d_log) + exp(n_a-d_log))
    }else{
      n_log = pnorm (y, mu, sigma, lower.tail=TRUE, log.p=TRUE)
      return (exp(n_log-d_log))
    }
  }
}


##### function for computing the mean of an N(mu, 1) RV
##### truncated to be on the interval (-Inf, a) \union (b, Inf)
tn.mean = function(mu, a, b, sigma=1){
  # denominator
  d_left = pnorm (a, mu, sigma, lower.tail=TRUE, log.p=TRUE)
  d_right = pnorm (b, mu, sigma, lower.tail=FALSE, log.p=TRUE)
  d_max = max(d_left, d_right)
  d_log = d_max + log(exp(d_left - d_max) + exp(d_right - d_max))
  
  # numerator
  n_left = dnorm (b, mu, sigma, log=TRUE)
  n_right = dnorm (a, mu, sigma, log=TRUE)
  
  if (n_left > n_right){
    mu + exp(n_left + log(1 - exp(n_right-n_left)) - d_log)
  }else{
    mu - exp(n_right + log(1 - exp(n_left-n_right)) - d_log)
  }
}




#### returns a pretty data frame summarising the information of an object of the mm class
#### columns for index, signal size estimate, left and right CI bounds and p values
#### only for those sample elements selected by the selection procedure associated with the mmObj
summary.manyMeans = function (object, ...){
  selected.set = object$selected.set
  mu.hat = object$muhat[selected.set]
  selint = object$selint # note that we need not apply the selected set
  p.vals=object$p.vals # note that we need not apply the selected set
  
  pretty.out = data.frame(index=selected.set, muhat=mu.hat,lowerConfPt=selint[,1], upperConfPt=selint[,2], p.value=p.vals)
  
  return (pretty.out)
}

#### function for printing an manyMeans object
#### calls the summary function, which makes a pretty data frame
#### then prints the data frame
print.manyMeans = function(x, ...){
  ## pretty print
  pretty.out = summary(x)
  print (pretty.out)
}








