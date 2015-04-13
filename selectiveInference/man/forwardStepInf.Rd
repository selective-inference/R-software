\name{forwardStepInf}
\alias{forwardStepInf}
\title{
Selective inference for forward stepwise regression
}
\description{
Computes exact p-values and selection intervals for forward stepwise regression
}
\usage{
forwardStepInf(fsfit, x, y, sigma, nsteps = NULL, alpha = 0.1, fixed.step = NULL, aic.stop = FALSE, trace = F, compute.ci = TRUE, one.sided = TRUE)
}
\arguments{
  \item{fsfit}{
Object returned by forwardStep
}
  \item{x}{
Matrix of predictors (n by p)
}
  \item{y}{
Vector of outcomes (length n)
}
\item{sigma}{
Estimate of error standard deviation; default is value returned by forwardStep
}
  \item{nsteps}{
Number of steps to analyze (default is number used in forwardStep)
}
  \item{alpha}{
significance level for  selection intervals (target is miscoverage alpha/2 in each tail)
}
  \item{fixed.step}{
If an integer k, results at step k are computed. If NULL and aic.stop is FALSE (default),
results for each FS step are computed. 
}
  \item{aic.stop}{
 Use fixed step prescibed by  AIC? Default FALSE
}
  \item{trace}{
Print out more details along the way.
}
  \item{compute.ci}{
Compute selection interals? Default TRUE. Can slow down the computation.
}
  \item{one.sided}{
One-sided p-values? Default TRUE.
}
}
\details{
}

\value{
  \item{pv}{p-values for each predictor}
  \item{vm}{Lower bounds for functionals}
 \item{vp}{Upper bounds for functionals}
\item{ci}{Selection intervals}
\item{miscov}{Realized tail area for each  selection interval}
\item{pred}{Predictors entered at each step of forward stepwise}
\item{which.steps }{Which FS steps were analyzed}
\item{step.ind}{Step indicator for each row of constraint matrix A (returned below)}
\item{alpha}{Desired coverage (alpha/2 in each tail)}
\item{sigma}{Value of error standard devarion (sigma) used}
\item{one.sided}{Are the pvalues one-sided}
\item{A}{Constraint matrix}
\item{b}{Right hand side vector for constraints}
\item{call}{The call to forwardStepInf}

}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Joshua Loftus, Jonathan Taylor, Ryan Tibshirani, Rob Tibshirani
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
