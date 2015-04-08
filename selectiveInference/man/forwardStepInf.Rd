\name{forwardStepInf}
\alias{forwardStepInf}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Selective inference for forward stepwise regression
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
forwardStepInf(fsfit, x, y, sigma, nsteps = NULL, alpha = 0.1, fixed.step = NULL, aic.stop = FALSE, trace = F, compute.ci = TRUE, one.sided = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fsfit}{
%%     ~~Describe \code{fsfit} here~~
}
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{y}{
%%     ~~Describe \code{y} here~~
}
  \item{sigma}{
%%     ~~Describe \code{sigma} here~~
}
  \item{nsteps}{
%%     ~~Describe \code{nsteps} here~~
}
  \item{alpha}{
%%     ~~Describe \code{alpha} here~~
}
  \item{fixed.step}{
%%     ~~Describe \code{fixed.step} here~~
}
  \item{aic.stop}{
%%     ~~Describe \code{aic.stop} here~~
}
  \item{trace}{
%%     ~~Describe \code{trace} here~~
}
  \item{compute.ci}{
%%     ~~Describe \code{compute.ci} here~~
}
  \item{one.sided}{
%%     ~~Describe \code{one.sided} here~~
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
function (fsfit, x, y, sigma, nsteps = NULL, alpha = 0.1, fixed.step = NULL, 
    aic.stop = FALSE, trace = F, compute.ci = TRUE, one.sided = TRUE) 
{
    this.call = match.call()
    n = nrow(x)
    p = ncol(x)
    if (aic.stop) {
        fixed.step = fsfit$aichat
    }
    if (is.null(nsteps)) {
        nsteps = length(fsfit$pred)
    }
    if (!is.null(fixed.step)) 
        nsteps = fixed.step
    x = scale(x, TRUE, FALSE)
    y = y - mean(y)
    SMALL = 1e-07
    pred = fsfit$pred
    s = fsfit$s
    a = NULL
    proj = 0
    which.steps = 1:nsteps
    if (!is.null(fixed.step)) {
        which.steps = rep(fixed.step, fixed.step)
    }
    if (aic.stop & fsfit$aichat == 0) {
        out = list(pv = rep(NA, p), ci = NA, A = NA, b = NA)
        cat("Null model picked by AIC; no  results returned", 
            fill = T)
        return(out)
    }
    xx = x
    vv = 1/sum(x[, pred[1]]^2)
    stepind = NULL
    A = b = NULL
    for (j in 1:nsteps) {
        if (trace) 
            cat(c("step=", j), fill = T)
        notin = rep(T, p)
        if (j > 1) {
            mod = pred[1:(j - 1)]
            notin[mod] = F
            xx[, -mod] = lsfit(xx[, mod], xx[, -mod])$res
            vv = solve(t(xx[, mod, drop = F]) \%*\% xx[, mod, drop = F])
            proj = x[, mod] \%*\% vv \%*\% t(xx[, mod])
        }
        jj = pred[j]
        rest = notin
        rest[jj] = F
        w = sqrt(sum(xx[, jj]^2))
        for (k in which(rest)) {
            w2 = sqrt(sum(xx[, k]^2))
            A = rbind(A, -(s[j] * xx[, jj]/w - xx[, k]/w2) \%*\% 
                (diag(n) - proj))
            A = rbind(A, -(s[j] * xx[, jj]/w + xx[, k]/w2) \%*\% 
                (diag(n) - proj))
            b = c(b, 0, 0)
            stepind = c(stepind, j, j)
        }
        if (aic.stop) {
            mod2 = pred[1:j]
            vv = solve(t(x[, mod2, drop = F]) \%*\% x[, mod2, drop = F])
            temp = vv \%*\% t(x[, mod2, drop = F])
            temp2 = -s[j] * temp[j, , drop = F]/(sigma * sqrt(vv[j, 
                j]))
            A = rbind(A, temp2)
            b = c(b, -sqrt(2))
        }
    }
    if (is.null(fixed.step)) {
        mod = pred[1:nsteps]
        fitter = solve(t(xx[, mod]) \%*\% xx[, mod]) \%*\% t(xx[, 
            mod])
        aa = -1 * fitter[nsteps, ] * s[nsteps]
        A = rbind(A, aa)
        stepind = c(stepind, nsteps)
        b = c(b, 0)
    }
    if (aic.stop) {
        mod2 = pred[1:(nsteps + 1)]
        vv = solve(t(x[, mod2, drop = F]) \%*\% x[, mod2, drop = F])
        pp = length(mod2)
        temp = vv \%*\% t(x[, mod2, drop = F])
        snew = sign(sum(temp[pp, ] * y))
        temp2 = snew * temp[pp, , drop = F]/(sigma * sqrt(vv[pp, 
            pp]))
        A = rbind(A, temp2)
        b = c(b, sqrt(2))
    }
    vmall = vpall = vector("list", nsteps)
    pv = rep(NA, p)
    ci = miscov = matrix(NA, nrow = p, ncol = 2)
    for (kk in 1:length(which.steps)) {
        if (trace) 
            cat(c("step=", kk), fill = T)
        if (is.null(fixed.step)) 
            pp = sum(stepind <= which.steps[kk])
        if (!is.null(fixed.step)) 
            pp = nrow(A)
        mod = pred[1:which.steps[kk]]
        temp = (solve(t(x[, mod, drop = F]) \%*\% x[, mod, drop = F]) \%*\% 
            t(x[, mod, drop = F]))
        kkk = nrow(temp)
        if (!is.null(fixed.step)) 
            kkk = kk
        eta = as.vector(temp[kkk, ])
        bhat = sum(eta * y)
        if (one.sided) 
            eta = eta * sign(bhat)
        alp = as.vector(A \%*\% eta/sum(eta^2))
        alp[abs(alp) < SMALL] = 0
        vp = rep(Inf, pp)
        vm = rep(-Inf, pp)
        for (j in 1:pp) {
            if (alp[j] < 0) 
                vm[j] = (b[j] - (A \%*\% y)[j] + alp[j] * sum(eta * 
                  y))/alp[j]
            if (alp[j] > 0) 
                vp[j] = (b[j] - (A \%*\% y)[j] + alp[j] * sum(eta * 
                  y))/alp[j]
        }
        vmm = max(vm, na.rm = T)
        vpp = min(vp, na.rm = T)
        vmall[[kk]][jj] = vmm
        vpall[[kk]][jj] = vpp
        tt = sum(eta * y)
        sigma.eta = sigma * sqrt(sum(eta^2))
        u = 0
        pv[kk] = 1 - ptruncnorm(tt, vmm, vpp, u, sigma.eta)
        if (!one.sided) 
            pv[kk] = 2 * min(pv[kk], 1 - pv[kk])
        if (compute.ci) {
            vs = list(vm = vmm, vp = vpp)
            junk = selection.int(y, eta, sigma^2, vs, alpha)
            ci[kk, ] = junk$ci
            miscov[kk, ] = junk$miscov
        }
    }
    out = list(pv = pv, vm = vmall, vp = vpall, ci = ci, miscov = miscov, 
        which.steps = which.steps, pred = fsfit$pred, stepind = stepind, 
        alpha = alpha, one.sided = one.sided, A = A, b = b, call = this.call)
    class(out) = "forwardStepInf"
    return(out)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
