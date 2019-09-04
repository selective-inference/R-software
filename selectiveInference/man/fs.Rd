\name{fs}
\alias{fs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Forward stepwise regression
}
\description{
This function implements forward stepwise regression, for use in the 
selectiveInference package
}
\usage{
fs(x, y, maxsteps=2000, intercept=TRUE, normalize=TRUE, verbose=FALSE)
}
\arguments{
  \item{x}{
Matrix of predictors (n by p)
}
  \item{y}{
Vector of outcomes (length n)
}
  \item{maxsteps}{
Maximum number of steps to take 
}
\item{intercept}{Should an intercept be included on the model? Default is TRUE}
\item{normalize}{Should the predictors be normalized? Default is TRUE. (Note:
this argument has no real effect on model selection since forward stepwise is 
scale invariant already; however, it is included for completeness, and to match
the interface for the \code{lar} function)
}
\item{verbose}{Print out progress along the way? Default is FALSE}
}

\details{
This function implements forward stepwise regression, adding the predictor at each 
step that maximizes the absolute correlation between the predictors---once 
orthogonalized with respect to the current model---and the residual. This entry
criterion is standard, and is equivalent to choosing the variable that achieves 
the biggest drop in RSS at each step; it is used, e.g., by the \code{step} function
in R. Note that, for example, the \code{lars} package implements a stepwise option 
(with type="step"), but uses a (mildly) different entry criterion, based on maximal
absolute correlation between the original (non-orthogonalized) predictors and the 
residual. 
}
\value{
\item{action}{Vector of predictors in order of entry}
\item{sign}{Signs of coefficients of predictors, upon entry}
\item{df}{Degrees of freedom of each active model}
\item{beta}{Matrix of regression coefficients for each model along the path,
one column per model}
\item{completepath}{Was the complete stepwise path computed?}
\item{bls}{If completepath is TRUE, the full least squares coefficients}
\item{Gamma}{Matrix that captures the polyhedral selection at each step}
\item{nk}{Number of polyhedral constraints at each step in path}
\item{vreg}{Matrix of linear contrasts that gives coefficients of variables
to enter along the path}
\item{x}{Matrix of predictors used}
\item{y}{Vector of outcomes used}
\item{bx}{Vector of column means of original x}
\item{by}{Mean of original y}
\item{sx}{Norm of each column of original x}
\item{intercept}{Was an intercept included?}
\item{normalize}{Were the predictors normalized?}
\item{call}{The call to fs}
}

\author{Ryan Tibshirani, Rob Tibshirani, Jonathan Taylor, Joshua Loftus, Stephen Reid}

\seealso{
 \code{\link{fsInf}}, \code{\link{predict.fs}},\code{\link{coef.fs}}, \code{\link{plot.fs}}
}

\examples{
set.seed(33)
n = 50
p = 10
sigma = 1
x = matrix(rnorm(n*p),n,p)
beta = c(3,2,rep(0,p-2))
y = x\%*\%beta + sigma*rnorm(n)

# run forward stepwise, plot results
fsfit = fs(x,y)
plot(fsfit)

# compute sequential p-values and confidence intervals
# (sigma estimated from full model)
out = fsInf(fsfit)
out
}


