\name{randomizedLassoInf}
\alias{randomizedLassoInf}

\title{
Inference for the randomized lasso, with a fixed lambda 
}
\description{
Compute p-values and confidence intervals based on selecting
an active set with the randomized lasso, at a 
fixed value of the tuning parameter lambda and using Gaussian
randomization.
}
\usage{
randomizedLassoInf(X, 
                   y, 
                   lam, 
                   sigma=NULL, 
                   noise_scale=NULL, 
                   ridge_term=NULL, 
                   condition_subgrad=TRUE, 
                   level=0.9,
                   nsample=10000,
                   burnin=2000,
                   max_iter=100,       
                   kkt_tol=1.e-4,      
                   parameter_tol=1.e-8,
                   objective_tol=1.e-8,
                   objective_stop=FALSE,
                   kkt_stop=TRUE,
                   parameter_stop=TRUE)
}
\arguments{
  \item{X}{
Matrix of predictors (n by p); 
}
  \item{y}{
Vector of outcomes (length n)
}
  \item{lam}{
Value of lambda used to compute beta. See the above warning
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
\item{sigma}{
Estimate of error standard deviation. If NULL (default), this is estimated 
using the mean squared residual of the full least squares based on 
selected active set.
}
\item{noise_scale}{
Scale of Gaussian noise added to objective. Default is 
0.5 * sd(y) times the sqrt of the mean of the trace of X^TX.
}
\item{ridge_term}{
A small "elastic net" or ridge penalty is added to ensure
the randomized problem has a solution. 
0.5 * sd(y) times the sqrt of the mean of the trace of X^TX divided by
sqrt(n).
}
\item{condition_subgrad}{
In forming selective confidence intervals and p-values should we condition
on the inactive coordinates of the subgradient as well?
Default is TRUE. 
}
\item{level}{
Level for confidence intervals.
}
\item{nsample}{
Number of samples of optimization variables to sample.
}
\item{burnin}{
How many samples of optimization variable to discard (should be less than nsample).
}
\item{max_iter}{
How many rounds of updates used of coordinate descent in solving randomized
LASSO.
}
\item{kkt_tol}{
Tolerance for checking convergence based on KKT conditions.
}
\item{parameter_tol}{
Tolerance for checking convergence based on convergence
of parameters.
}
\item{objective_tol}{
Tolerance for checking convergence based on convergence
of objective value.
}
\item{kkt_stop}{
Should we use KKT check to determine when to stop?
}
\item{parameter_stop}{
Should we use convergence of parameters to determine when to stop?
}
\item{objective_stop}{
Should we use convergence of objective value to determine when to stop?
}
}

\details{
This function computes selective p-values and confidence intervals for a
randomized version of the lasso,
given a fixed value of the tuning parameter lambda. 

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
\item{call}{The call to lassoInf}
}

\references{
Xiaoying Tian, and Jonathan Taylor (2015).
Selective inference with a randomized response. arxiv.org:1507.06739

Xiaoying Tian, Snigdha Panigrahi, Jelena Markovic, Nan Bi and Jonathan Taylor (2016).
Selective inference after solving a convex problem. 
arxiv:1609.05609

}
\author{Jelena Markovic, Jonathan Taylor}

\examples{
set.seed(43)
n = 50
p = 10
sigma = 0.2
lam = 0.5

X = matrix(rnorm(n*p), n, p)
X = scale(X, TRUE, TRUE) / sqrt(n-1)

beta = c(3,2,rep(0,p-2))
y = X\%*\%beta + sigma*rnorm(n)

result = randomizedLassoInf(X, y, lam)

}
 
