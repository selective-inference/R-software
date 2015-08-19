\name{lars.glm}
\alias{lars.glm}
\title{
Function to fit least angle regression path of solution for Ll-penalized (lasso) logistic regression
and the Cox proportional hazards model.
}
\description{
Function to fit least angle regression path of solution for Ll-penalized (lasso) logistic regression
and the Cox proportional hazards model.
}
\usage{
lars.glm(x, y, status = NULL, family = c("binomial", "cox"), standardize=TRUE,frac.arclength=.1)
}
\arguments{
 \item{x}{
N by p matrix of predictors
}
  \item{y}{
N-vector of outcome values
}
  \item{status}{
Optional N-vector of censoring indicators for Cox Proportioanl hazards model. 1=failed; 0=censored.
}
  \item{family}{
   "binomial" or "cox"
}
\item{standardize}{
Should predictor be standardized? Default TRUE
}
  \item{frac.arclength}{
   Step length parameter for \code{glmpath}. See help file for \code{glmpath}
   for details
}
}
\details{
This function estimates the  least angle regression path of solution for Ll-penalized (lasso) logistic regression
and the Cox proportional hazards model, using the R functions \code{glmpath} and \code{coxpath} written by Park and Hastie.
These latter functions use the predictor-corrector strategy devised by Park and Hastie (2007).
An additional L2 penalty can be added for stability.
}
\value{
  \item{beta }{Matrix of estimated coefficients, with LAR steps in the rows.}
  \item{a0}{Estimate of intercept}
  \item{lambda0}{Raw values of lambda used}
  \item{lambda}{Values of lambda multiplied by sdx, the standard deviation of each predictor}
  \item{lambda2}{Value of lambda2 (L2 penalty parameter)} 
  \item{act}{Actions (predictor added) at each step}
  \item{maxp}{Maximum number of predictors entered}
\item{family}{family used- "binomial" or "cox"}
\item{call}{Call to lars.glm}
\item{pathobj}{Result of call to  glmpath or coxpath}
}
\references{
Park, M.Y. and Hastie, T. (2007) 1l regularization path algorithm for generalized linear
models. JRSSB B 69(4), 659-677

}
\author{
Rob Tibshirani
}


\seealso{
covTest, predict.glm.Rd
}
\examples{
x=matrix(rnorm(100*10),ncol=10)
x=scale(x,TRUE,TRUE)/sqrt(99)

#logistic
y=5*x[,2]+rnorm(100)
y=1*(y>0)
a=lars.glm(x,y,family="binomial")


# Cox model
y=6*x[,2]+rnorm(100)+10
status=sample(c(0,1),size=length(y),replace=TRUE)
a=lars.glm(x,y,status=status,family="cox")
  }
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ LARS }
\keyword{ lasso }% __ONLY ONE__ keyword per line
\keyword{ logistic }% __ONLY ONE__ keyword per line
\keyword{ Cox model }% __ONLY ONE__ keyword per line
