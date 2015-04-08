\name{forwardStep}
\alias{forwardStep}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Forward stepwise regression
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
forwardStep(x, y, nsteps = min(nrow(x), ncol(x)))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{y}{
%%     ~~Describe \code{y} here~~
}
  \item{nsteps}{
%%     ~~Describe \code{nsteps} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x, y, nsteps = min(nrow(x), ncol(x))) 
{
    this.call = match.call()
    BIG = 1e+10
    p = ncol(x)
    x = scale(x, T, F)
    y = y - mean(y)
    pred = s = scor = bhat = rep(NA, nsteps)
    ip = t(x) \%*\% y/sqrt(diag(t(x) \%*\% x))
    pred[1] = which.max(abs(ip))
    s[1] = sign(sum(x[, pred[1]] * y))
    scor[1] = ip[pred[1]]
    bhat[1] = ip[pred[1]]/sqrt(sum(x[, pred[1]]^2))
    r = lsfit(x[, pred[1]], y)$res
    for (j in 2:nsteps) {
        mod = pred[1:(j - 1)]
        r = lsfit(x[, mod], r)$res
        xr = lsfit(x[, mod], x)$res
        ip = t(xr) \%*\% r/sqrt(diag(t(xr) \%*\% xr))
        ip[mod] = 0
        pred[j] = which.max(abs(ip))
        scor[j] = ip[pred[j]]
        s[j] = sign(sum(xr[, pred[j]] * r))
        bhat[j] = ip[pred[j]]/sqrt(sum(xr[, pred[j]]^2))
    }
    aic = rss = val = rep(NA, nsteps + 1)
    rss[1] = sum((y - mean(y))^2)
    aic[1] = rss[1]
    for (j in 1:nsteps) {
        xx = x[, pred[1:j], drop = F]
        rss[j + 1] = sum(lsfit(xx, y)$res^2)
        aic[j + 1] = rss[j + 1] + 2 * j * sigma^2
    }
    d = diff(aic)
    o = which(d > 0)
    aichat = o[1] - 1
    out = list(pred = pred, s = s, scor = scor, bhat = bhat, 
        rss = rss, aic = aic, aichat = aichat, call = this.call)
    class(out) = "forwardStep"
    return(out)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
