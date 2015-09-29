\name{plot.fs}
\alias{plot.fs}

\title{
Plot function for forward stepwise regression
}
\description{
Plot coefficient profiles along the forward stepwise  path
}

\usage{
\method{plot}{fs} (x, breaks=TRUE, omit.zeros=TRUE, var.labels=TRUE, ...)
}

\arguments{
\item{x}{
Object returned by a call to \code{fs} function
}
\item{breaks}{Should vertical lines be drawn at each break point in the piecewise
linear coefficient paths? Default is TRUE}
\item{omit.zeros}{Should segments of the coefficients paths that are equal to 
zero be omitted (to avoid clutter in the figure)? Default is TRUE}
\item{var.labels}{Should paths be labelled with corresponding variable numbers? 
Default is TRUE}
\item{...}{Additional arguments for plotting}
}

\author{Ryan Tibshirani, Rob Tibshirani, Jonathan Taylor, Joshua Loftus, Stephen Reid}

\examples{
#NOT RUN
#set.seed(33)
#n=200
#p=20
#sigma=1
#x=matrix(rnorm(n*p),n,p)
#x=scale(x,T,T)/sqrt(n-1)
#generate y
#beta=c(3,-2,rep(0,p-2))
#beta=c(rep(3,10),rep(0,p-10))
#y=x%*%beta+sigma*rnorm(n)
#y=y-mean(y)
#fsfit=fs(x,y)
#plot(fsfit)
}
