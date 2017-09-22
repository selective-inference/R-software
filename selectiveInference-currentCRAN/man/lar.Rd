\name{lar}
\alias{lar}
\title{
Least angle regression
}
\description{
This function implements least angle regression, for use in the 
selectiveInference package
}
\usage{
lar(x, y, maxsteps=2000, minlam=0, intercept=TRUE, normalize=TRUE,
    verbose=FALSE) 
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
\item{minlam}{
Minimum value of lambda to consider
}
\item{intercept}{Should an intercept be included on the model? Default is TRUE}
\item{normalize}{Should the predictors be normalized? Default is TRUE}
\item{verbose}{Print out progress along the way? Default is FALSE}
}

\details{
The least angle regression algorithm is described in detail by Efron et al. (2002). 
This function should match (in terms of its output) that from the \code{lars} package,
but returns additional information (namely, the polyhedral constraints) needed for the
selective inference calculations.
}

\value{
\item{lambda}{Values of lambda (knots) visited along the path}
\item{action}{Vector of predictors in order of entry}
\item{sign}{Signs of coefficients of predictors, upon entry}
\item{df}{Degrees of freedom of each active model}
\item{beta}{Matrix of regression coefficients for each model along the path,
one model per column}
\item{completepath}{Was the complete stepwise path computed?}
\item{bls}{If completepath is TRUE, the full least squares coefficients}
\item{Gamma}{Matrix that captures the polyhedral selection at each step}
\item{nk}{Number of polyhedral constraints at each step in path}
\item{vreg}{Matrix of linear contrasts that gives coefficients of variables
to enter along the path}     
\item{mp}{Value of M+ (for internal use with the spacing test)}	 
\item{x}{Matrix of predictors used}
\item{y}{Vector of outcomes used}
\item{bx}{Vector of column means of original x}
\item{by}{Mean of original y}
\item{sx}{Norm of each column of original x}
\item{intercept}{Was an intercept included?}
\item{normalize}{Were the predictors normalized?}
\item{call}{The call to lar}
}

\references{
Brad Efron, Trevor Hastie, Iain Johnstone, and Rob Tibshirani (2002).
Least angle regression. Annals of Statistics (with discussion).

See also the descriptions in Trevor Hastie, Rob Tibshirani, and 
Jerome Friedman (2002, 2009). Elements of Statistical Learning.
}

\author{Ryan Tibshirani, Rob Tibshirani, Jonathan Taylor, Max G'Sell, Joshua Loftus, Stephen Reid}

\seealso{
 \code{\link{larInf}}, \code{\link{predict.lar}}, \code{\link{coef.lar}}, \code{\link{plot.lar}}
}

\examples{
set.seed(43)
n = 50
p = 10
sigma = 1
x = matrix(rnorm(n*p),n,p)
beta = c(3,2,rep(0,p-2))
y = x\%*\%beta + sigma*rnorm(n)

# run LAR, plot results
larfit = lar(x,y)
plot(larfit)

# compute sequential p-values and confidence intervals
# (sigma estimated from full model)
out = larInf(larfit)
out                                    
}
