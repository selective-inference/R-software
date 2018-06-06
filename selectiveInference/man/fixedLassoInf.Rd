\name{fixedLassoInf}
\alias{fixedLassoInf}

\title{
Inference for the lasso, with a fixed lambda 
}
\description{
Compute p-values and confidence intervals for the lasso estimate, at a 
fixed value of the tuning parameter lambda
}
\usage{
fixedLassoInf(x, 
              y, 
              beta, 
              lambda, 
              family = c("gaussian", "binomial", "cox"),
              intercept=TRUE, 
              add.targets=NULL, 
              status=NULL, 
              sigma=NULL, 
              alpha=0.1,
              type=c("partial","full"), 
              tol.beta=1e-5, 
              tol.kkt=0.1,
              gridrange=c(-100,100), 
              bits=NULL, 
              verbose=FALSE, 
              linesearch.try=10) 
}
\arguments{
  \item{x}{
Matrix of predictors (n by p); 
}
  \item{y}{
Vector of outcomes (length n)
}
  \item{beta}{
Estimated lasso coefficients (e.g., from glmnet). This is of length p
(so the intercept is not included as the first component).

 Be careful! This function uses the "standard" lasso objective
  \deqn{
    1/2 \|y - x \beta\|_2^2 + \lambda \|\beta\|_1.
  }
 In contrast, glmnet multiplies the first term by a factor of 1/n.
 So after running glmnet, to extract the beta corresponding to a value lambda, 
 you need to use \code{beta = coef(obj, s=lambda/n)[-1]},
 where obj is the object returned by glmnet (and [-1] removes the intercept,
 which glmnet always puts in the first component)
}      
  \item{lambda}{
Value of lambda used to compute beta. See the above warning
} 

\item{family}{Response type: "gaussian" (default), "binomial", or 
    "cox" (for censored survival data) }

\item{sigma}{
Estimate of error standard deviation. If NULL (default), this is estimated 
using the mean squared residual of the full least squares fit when n >= 2p, and 
using the standard deviation of y when n < 2p. In the latter case, the user 
should use \code{\link{estimateSigma}} function for a more accurate estimate.
Not used for family= "binomial", or  "cox"
}
\item{alpha}{
Significance level for confidence intervals (target is miscoverage alpha/2 in each tail)
}
\item{intercept}{
Was the lasso problem solved (e.g., by glmnet) with an intercept in the model? 
Default is TRUE. Must be TRUE for "binomial" family. Not used for 'cox" family, where no intercept is assumed.
}
\item{add.targets}{Optional vector of predictors to be included as targets of inference, regardless of whether or not they are selected by the lasso. Default is NULL.}
\item{status}{Censoring status for Cox model; 1=failurem 0=censored}
\item{type}{Contrast type for p-values and confidence intervals: default is
"partial"---meaning that the contrasts tested are the partial population 
regression coefficients, within the active set of predictors; the alternative is
"full"---meaning that the full population regression coefficients are tested.
The latter does not make sense when p > n.}
\item{tol.beta}{
Tolerance for determining if a coefficient is zero
}
\item{tol.kkt}{
Tolerance for determining if an entry of the subgradient is zero 
}
\item{gridrange}{
Grid range for constructing confidence intervals, on the standardized scale
}
\item{bits}{
Number of bits to be used for p-value and confidence interval calculations. Default is 
NULL, in which case standard floating point calculations are performed. When not NULL, 
multiple precision floating point calculations are performed with the specified number 
of bits, using the R package \code{Rmpfr} (if this package is not installed, then a 
warning is thrown, and standard floating point calculations are pursued).
Note: standard double precision uses 53 bits
so, e.g., a choice of 200 bits uses about 4 times double precision. The confidence
interval computation is sometimes numerically challenging, and the extra precision can be 
helpful (though computationally more costly). In particular, extra precision might be tried 
if the values in the output columns of \code{tailarea} differ noticeably from alpha/2.
}
\item{verbose}{
Print out progress along the way? Default is FALSE
}
\item{linesearch.try}{
When running type="full" (i.e. debiased LASSO) how many attempts in the line search?
}
}

\details{
This function computes selective p-values and confidence intervals for the lasso,
given a fixed value of the tuning parameter lambda. 
Three different response types are supported:  gaussian, binomial and Cox.
The confidence interval construction involves numerical search and can be fragile:
if the observed statistic is too close to either end of the truncation interval
(vlo and vup, see references), then one or possibly both endpoints of the interval of 
desired coverage cannot be computed, and default to +/- Inf. The output \code{tailarea} 
gives the achieved Gaussian tail areas for the reported intervals---these should be close 
to alpha/2, and can be used for error-checking purposes.

Important!: Before running glmnet (or some other lasso-solver)  x should be centered, that is x <- scale(X,TRUE,FALSE).
In addition, if standardization of the predictors is desired, x should be scaled as well:  x <- scale(x,TRUE,TRUE).
Then when running glmnet, set standardize=F. See example below.

The penalty.factor facility in glmmet-- allowing different penalties lambda for each predictor,
is not yet implemented in fixedLassoInf. However you can finesse this--- see the example below. One caveat- using this approach, a penalty factor of zero (forcing a predictor in)
is not allowed.

Note that the coefficients and standard errors reported are unregularized.
Eg for the Gaussian, they are the usual least squares estimates and standard errors
for the model fit to the active set from the lasso.
}
\value{  
\item{type}{Type of coefficients tested (partial or full)}
\item{lambda}{Value of tuning parameter lambda used}
\item{pv}{One-sided P-values for active variables, uses the fact we have conditioned on the sign.}
\item{ci}{Confidence intervals}
\item{tailarea}{Realized tail areas (lower and upper) for each confidence interval}
\item{vlo}{Lower truncation limits for statistics}
\item{vup}{Upper truncation limits for statistics}
\item{vmat}{Linear contrasts that define the observed statistics}
\item{y}{Vector of outcomes}
\item{vars}{Variables in active set}
\item{sign}{Signs of active coefficients}
\item{alpha}{Desired coverage (alpha/2 in each tail)}
\item{sigma}{Value of error standard deviation (sigma) used}
\item{call}{The call to fixedLassoInf}
}

\references{
Jason Lee, Dennis Sun, Yuekai Sun, and Jonathan Taylor (2013). 
Exact post-selection inference, with application to the lasso. arXiv:1311.6238.

 Jonathan Taylor and Robert Tibshirani (2016) Post-selection inference for L1-penalized likelihood models.
arXiv:1602.07358

}
\author{Ryan Tibshirani, Rob Tibshirani, Jonathan Taylor, Joshua Loftus, Stephen Reid}

\examples{
set.seed(43)
n = 50
p = 10
sigma = 1

x = matrix(rnorm(n*p),n,p)
x = scale(x,TRUE,TRUE)

beta = c(3,2,rep(0,p-2))
y = x\%*\%beta + sigma*rnorm(n)

# first run glmnet
gfit = glmnet(x,y,standardize=FALSE)

# extract coef for a given lambda; note the 1/n factor!
# (and we don't save the intercept term)
lambda = .8
beta = coef(gfit, x=x, y=y, s=lambda/n, exact=TRUE)[-1]

# compute fixed lambda p-values and selection intervals
out = fixedLassoInf(x,y,beta,lambda,sigma=sigma)
out


## as above, but use lar function instead to get initial 
## lasso fit (should get same results)
lfit = lar(x,y,normalize=FALSE)
beta = coef(lfit, s=lambda, mode="lambda")
out2 = fixedLassoInf(x, y, beta, lambda, sigma=sigma)
out2

## mimic different penalty factors by first scaling x
 set.seed(43)
n = 50
p = 10
sigma = 1

x = matrix(rnorm(n*p),n,p)
x=scale(x,TRUE,TRUE)

beta = c(3,2,rep(0,p-2))
y = x\%*\%beta + sigma*rnorm(n)
pf=c(rep(1,7),rep(.1,3))  #define penalty factors
pf=p*pf/sum(pf)   # penalty factors should be rescaled so they sum to p
xs=scale(x,FALSE,pf) #scale cols of x by penalty factors
# first run glmnet
gfit = glmnet(xs, y, standardize=FALSE)

# extract coef for a given lambda; note the 1/n factor!
# (and we don't save the intercept term)
lambda = .8
beta_hat = coef(gfit, x=xs, y=y, s=lambda/n, exact=TRUE)[-1]

# compute fixed lambda p-values and selection intervals
out = fixedLassoInf(xs,y,beta_hat,lambda,sigma=sigma)

#rescale conf points to undo the penalty factor
out$ci=t(scale(t(out$ci),FALSE,pf[out$vars]))
out

#logistic model
set.seed(43)

n = 50
p = 10
sigma = 1

x = matrix(rnorm(n*p),n,p)
x=scale(x,TRUE,TRUE)

beta = c(3,2,rep(0,p-2))
y = x\%*\%beta + sigma*rnorm(n)
y=1*(y>mean(y))
# first run glmnet
gfit = glmnet(x,y,standardize=FALSE,family="binomial")

# extract coef for a given lambda; note the 1/n factor!
# (and here  we DO  include the intercept term)
lambda = .8
beta_hat = coef(gfit, x=x, y=y, s=lambda/n, exact=TRUE)

# compute fixed lambda p-values and selection intervals
out = fixedLassoInf(x,y,beta_hat,lambda,family="binomial")
out


# Cox model

set.seed(43)
n = 50
p = 10
sigma = 1

x = matrix(rnorm(n*p), n, p)
x = scale(x, TRUE, TRUE)

beta = c(3,2,rep(0,p-2))
tim = as.vector(x\%*\%beta + sigma*rnorm(n))
tim= tim-min(tim)+1
status=sample(c(0,1),size=n,replace=TRUE)
# first run glmnet


y = Surv(tim,status)
gfit = glmnet(x, y, standardize=FALSE, family="cox")

# extract coef for a given lambda; note the 1/n factor!

lambda = 1.5
beta_hat = as.numeric(coef(gfit, x=x, y=y, s=lambda/n, exact=TRUE))

# compute fixed lambda p-values and selection intervals
out = fixedLassoInf(x, tim, beta_hat, lambda, status=status, family="cox")
out

# Debiased lasso or "full"

n = 50
p = 100
sigma = 1

x = matrix(rnorm(n*p),n,p)
x = scale(x,TRUE,TRUE)

beta = c(3,2,rep(0,p-2))
y = x\%*\%beta + sigma*rnorm(n)

# first run glmnet
gfit = glmnet(x, y, standardize=FALSE, intercept=FALSE)

# extract coef for a given lambda; note the 1/n factor!
# (and we don't save the intercept term)
lambda = 2.8
beta = coef(gfit, x=x, y=y, s=lambda/n, exact=TRUE)[-1]

# compute fixed lambda p-values and selection intervals
out = fixedLassoInf(x, y, beta, lambda, sigma=sigma, type='full', intercept=FALSE)
out

# When n > p and "full" we use the full inverse
# instead of Javanmard and Montanari's approximate inverse

n = 200
p = 50
sigma = 1

x = matrix(rnorm(n*p),n,p)
x = scale(x,TRUE,TRUE)

beta = c(3,2,rep(0,p-2))
y = x\%*\%beta + sigma*rnorm(n)

# first run glmnet
gfit = glmnet(x, y, standardize=FALSE, intercept=FALSE)

# extract coef for a given lambda; note the 1/n factor!
# (and we don't save the intercept term)
lambda = 2.8
beta = coef(gfit, x=x, y=y, s=lambda/n, exact=TRUE)[-1]

# compute fixed lambda p-values and selection intervals
out = fixedLassoInf(x, y, beta, lambda, sigma=sigma, type='full', intercept=FALSE)
out

}
 
