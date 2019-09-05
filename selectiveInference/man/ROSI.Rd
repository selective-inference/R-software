\name{ROSI}
\alias{ROSI}
\title{
Relevant One-step Selective Inference for the LASSO
}
\description{
Compute p-values and confidence intervals for the lasso estimate, at a 
fixed value of the tuning parameter lambda using the "relevant"
conditioning event of arxiv.org/1801.09037. 
}
\usage{
ROSI(X, 
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
     use_debiased=TRUE) 
}
\arguments{
  \item{X}{
Matrix of predictors (n by p); 
}
  \item{y}{
Vector of outcomes (length n)
}
  \item{soln}{
Estimated lasso coefficients (e.g., from glmnet). This is of length p
(so the intercept is not included as the first component).

 Be careful! This function uses the "standard" lasso objective
  \deqn{
    1/2 \|y - X \beta\|_2^2 + \lambda \|\beta\|_1.
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
\item{penalty_factor}{
Penalty factor as used by glmnet.
Actual penalty used in solving the problem is 
\deqn{
\lambda \cdot \sum_{i=1}^p f_i |\beta_i|
}
with f being the penalty_factor. Defaults to vector of 1s.
}
\item{dispersion}{
Estimate of dispersion in the GLM. Can be taken to be
1 for logisitic and should be an estimate of the error variance
for the Gaussian.
}
\item{family}{
Family used for likelihood.
}
\item{solver}{
Solver used to solve restricted problems needed to find truncation set.
Each active variable requires solving a new LASSO problem obtained
by zeroing out one coordinate of original problem. The "QP" choice
uses coordinate descent for a specific value of lambda, rather than glmnet which
would solve for a new path each time.
}
\item{construct_ci}{
Report confidence intervals or just p-values?
}
\item{debiasing_method}{
Which method should be used for debiasing? Choices are "JM" (Javanmard, Montanari)
or "BN" (method described in arxiv.org/1703.03282).
}
\item{verbose}{
Print out progress along the way? Default is FALSE.
}
\item{level}{
Confidence level for intervals.
}
\item{use_debiased}{
Use the debiased estimate of the parameter or not. When FALSE, this is the
method desribed in arxiv.org/1801.09037. The default TRUE often
produces noticably shorter intervals and more powerful tests when
p is comparable to n. Ignored when n<p and set to TRUE.
Also note that with "BN" as debiasing method and n > p, this agrees with
method in arxiv.org/1801.09037.
}
}

\details{
???
}
\value{  
\item{active_set}{Active set of LASSO.}
\item{pvalues}{Two-sided P-values for active variables.}
\item{intervals}{Confidence intervals}
\item{estimate}{Relaxed (i.e. unshrunk) selected estimates.}
\item{std_err}{Standard error of relaxed estimates (pre-selection).}
\item{dispersion}{Dispersion parameter.}
\item{lower_trunc}{Lower truncation point. The estimates should be outside the interval formed by the lower and upper truncation poitns.}
\item{upper_trunc}{Lower truncation point. The estimates should be outside the interval formed by the lower and upper truncation poitns.}
\item{lambda}{Value of tuning parameter lambda used.}
\item{penalty_factor}{Penalty factor used for solving problem.}
\item{level}{Confidence level.}
\item{call}{The call to fixedLassoInf.}
}

\references{


Keli Liu, Jelena Markovic, Robert Tibshirani. More powerful post-selection 
inference, with application to the Lasso. arXiv:1801.09037

Tom Boot, Didier Nibbering. Inference in high-dimensional linear regression models.
arXiv:1703.03282

}
\author{Jelena Markovic, Jonathan Taylor}

\examples{

\donttest{
library(selectiveInference)
library(glmnet)
set.seed(43)

n = 100
p = 200
s = 2
sigma = 1

x = matrix(rnorm(n*p),n,p)
x = scale(x,TRUE,TRUE)

beta = c(rep(10, s), rep(0,p-s)) / sqrt(n)
y = x \%*\% beta + sigma*rnorm(n)

# first run glmnet
gfit = glmnet(x,y,standardize=FALSE)

# extract coef for a given lambda; note the 1/n factor!
# (and we don't save the intercept term)
lambda = 4 * sqrt(n)
lambda_glmnet = 4 / sqrt(n)
beta = selectiveInference:::solve_problem_glmnet(x, 
                                                 y, 
                                                 lambda_glmnet, 
                                                 penalty_factor=rep(1, p),
                                                 family="gaussian")
# compute fixed lambda p-values and selection intervals
out = ROSI(x,
           y,
           beta,
           lambda,
           dispersion=sigma^2)
out

# an alternate approximate inverse from Boot and Nibbering

out = ROSI(x,
           y,
           beta,
           lambda,
           dispersion=sigma^2,
           debiasing_method="BN")
out
}
}
 
