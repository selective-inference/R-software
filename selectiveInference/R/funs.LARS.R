 predict.lar=
function (object, newx, s, type = c("fit", "coefficients"), mode = c("step", 
    "fraction", "norm", "lambda"), ...) 
{
    mode <- match.arg(mode)
    type <- match.arg(type)

   
    if (missing(newx) & type == "fit") {
        warning("Type=fit with no newx argument; type switched to coefficients")
        type <- "coefficients"
    }
    betas <- object$beta
    
    if (object$type != "LASSO" && mode %in% c("fraction", "norm")) 
        betas = betabreaker(object)
    dimnames(betas) = list(NULL, dimnames(betas)[[2]])
    sbetas <- scale(betas, FALSE, 1/object$normx)
    kp <- dim(betas)
    k <- kp[1]
    p <- kp[2]
    steps <- seq(k)
    if (missing(s)) {
        s <- steps
        mode <- "step"
    }
    sbeta <- switch(mode, step = {
        if (any(s < 0) | any(s > k)) stop("Argument s out of range")
        steps
    }, fraction = {
        if (any(s > 1) | any(s < 0)) stop("Argument s out of range")
        nbeta <- drop(abs(sbetas) %*% rep(1, p))
        nbeta/nbeta[k]
    }, norm = {
        nbeta <- drop(abs(sbetas) %*% rep(1, p))
        if (any(s > nbeta[k]) | any(s < 0)) stop("Argument s out of range")
        nbeta
    }, lambda = {
        lambdas = object$lambda
        s[s > max(lambdas)] = max(lambdas)
        s[s < 0] = 0
        c(lambdas, 0)
    })
    sfrac <- (s - sbeta[1])/(sbeta[k] - sbeta[1])
    sbeta <- (sbeta - sbeta[1])/(sbeta[k] - sbeta[1])
    usbeta <- unique(sbeta)
    useq <- match(usbeta, sbeta)
    sbeta <- sbeta[useq]
    betas <- betas[useq, , drop = FALSE]
    coord <- approx(sbeta, seq(sbeta), sfrac)$y
    left <- floor(coord)
    right <- ceiling(coord)
    newbetas <- ((sbeta[right] - sfrac) * betas[left, , drop = FALSE] + 
        (sfrac - sbeta[left]) * betas[right, , drop = FALSE])/(sbeta[right] - 
        sbeta[left])
    newbetas[left == right, ] <- betas[left[left == right], ]
    robject <- switch(type, coefficients = list(s = s, fraction = sfrac, 
        mode = mode, coefficients = drop(newbetas)), fit = list(s = s, 
        fraction = sfrac, mode = mode, fit = drop(scale(newx, 
            object$meanx, FALSE) %*% t(newbetas)) + object$mu))
    robject
}



larInf=function(x,y,larfit,sigma,compute.si=TRUE,alpha=.10,one.sided=TRUE,nsteps = min(nrow(x), ncol(x))){
    this.call=match.call()
    
       checkargs(x=x,y=y,larfit=larfit,sigma=sigma,alpha=alpha,nsteps=nsteps)
SMALL=1e-6
p=ncol(x)
    n=nrow(x)
nk=larfit$nk
     if(is.null(sigma) & p>=n){cat("Warning: p ge n; sigma set to 1",fill=T)
                           sigma=1}
  if(is.null(sigma) & n>p){
      sigma=sqrt(sum(lsfit(x,y)$res^2)/(n-p))
  }
vmm=vpp=pv=sigma.eta=rep(NA,nsteps)
ci=miscov=matrix(NA,nsteps,2)
for(k in 1:nsteps){
    mod=larfit$act[1:k]
 temp=(solve(t(x[,mod,drop=F])%*%x[,mod,drop=F])%*%t(x[,mod,drop=F]))
    temp=temp[nrow(temp),]
    bhat=sum(temp*y)
  eta=as.vector(temp)
    if(one.sided)eta=eta*sign(bhat)
    flip=(one.sided & sign(bhat)==-1)
    tt=sum(eta*y)
    A= -larfit$Gam[1:nk[k],]
    pp=nrow(A)
    b=rep(0,pp)
  
    vs = compute.vmvp(y, eta, A, b, pp)
     vpp = vs$vp
     vmm = vs$vm
   
    sigma.eta=sigma*sqrt(sum(eta^2))
       u=0  #null
pv[k]=1-rob.ptruncnorm(tt, vmm, vpp, u, sigma.eta)
if(!one.sided)  pv[k]=2*min(pv[k],1-pv[k])
  if(compute.si)
      {
          vs=list(vm=vmm,vp=vpp)
           junk=selection.int(y,eta,sigma, vs, alpha,flip=flip)
          ci[k,]=junk$ci
          miscov[k,]=junk$miscov
      }
 
}

pv.spacing=spacing.pval.asymp.list(y,larfit,nsteps,sigma=sigma)
   junk=covtest(larfit,x,y,sigma,nsteps)
    pv.cov=1-pexp(junk,1)

     forwardStopHat=forwardStop(pv,alpha)
    out=list(pv=pv,ci=ci,tailarea=miscov,vm=vmm,vp=vpp,sigma=sigma,alpha=alpha,act=larfit$act,pv.spacing=pv.spacing, pv.cov=pv.cov, forwardStopHat= forwardStopHat)
    out$call=this.call
   class(out)="larInference"
return(out)
}


print.larInference=function(x,digits = max(3, getOption("digits") - 3),...){
      cat("\nCall: ", deparse(x$call), "\n\n")
      cat(c("alpha=",x$alpha),fill=T)
      cat("",fill=T)
      nsteps=length(x$pv)
tab=cbind(1:nsteps,x$act[1:nsteps],round(x$pv,3),round(x$ci,3),round(x$tailarea,3),round(x$pv.spacing,3),round(x$pv.cov,3))
      dimnames(tab)=list(NULL,c("step","pred","exactPv","lowerConfPt","upperConfPt","lowerArea","upperArea","spacingPv","covtestPv"))
      print(tab)
      cat("",fill=T)
cat(c("Estimated stopping point from forwardStop rule=", x$forwardStopHat),fill=T)
         cat("",fill=T)
  }
covtest=function(fitobj,x,y,sigma,nsteps = min(nrow(x), ncol(x))) {

    n = nrow(x)
    p = ncol(x)
    my = mean(y)
      betas <- fitobj$beta
    lambda.min.ratio = ifelse(nrow(x) < ncol(x), 0.1, 1e-04)
    jlist = unlist(fitobj$act)
        lamlist = c(fitobj$lambda, 0)
    nsteps.call = nsteps
    nsteps = length(jlist)
    nsteps = min(nsteps, which(lamlist == 0))
    nsteps = min(nsteps, nsteps.call)
    jlist = jlist[1:nsteps]
    cov0 = cov = sig = rep(NA, nsteps)
    yy = y - my
    for (j in 1:nsteps) {
            lambda = lamlist[j + 1]
                yhat = predict.lar(fitobj, x, s = lambda, type = "fit", mode = "lam")$fit   
            cov[j] = sum(yy * yhat)
            if (j == 1) {
                cov0[j] = 0
            }
            if (j > 1) {
                tt0 = which(betas[j, ] != 0)
             #     aa = update(fitobj, x = x[, tt0, drop = F])
          # aa=lar(y, x = x[, tt0, drop = F], normalize=fitobj$normalize)# use this once ryan adds the option to lar
                 aa=lar(x[, tt0, drop = F], y,normalize=fitobj$normalize,intercept=fitobj$intercept)
       #         aa$beta=cbind(aa$beta,lsfit(x[,tt0,drop=F],y)$coef[-1]) # NOTE these and next few  only needed until ryans fixes lar
    #           aa$meanx=colMeans(x[,tt0,drop=F]);aa$normx=rep(1,length(tt0))
    #            aa$mu=mean(y)
   #             aa$normalize=FALSE
                  yhat0 = predict.lar(aa, x[, tt0,drop=F], type = "fit", s = lambda, mode = "lam")$fit
                cov0[j] = sum(yy * yhat0)
            }}
        
   
  
    tt = ((cov - cov0)/sigma^2)
    
return(tt)
}

# Asymptotic spacing LAR p values list
spacing.pval.asymp.list= function(y,out,k,sigma=1) {
  if (length(out$lambda)==1) stop("The LAR path needs to be run for at least 2 steps.")
  pvals = numeric(k)
  for (i in 1:k) {
    v = out$Gamma[out$nk[i],]
    denom = sigma*sqrt(sum(v*v))
    if (i==1 && i<length(out$lambda)) {
      pvals[1] = (1-pnorm(out$lambda[i]/denom))/
        (1-pnorm(out$lambda[i+1]/denom))
    }
    else if (i<length(out$lambda)) {
    #  pvals[i] = (pnorm(out$lambda[i-1]/denom)-pnorm(out$lambda[i]/denom))/
     #   (pnorm(out$lambda[i-1]/denom)-pnorm(out$lambda[i+1]/denom))
        pvals[i]=1-ptruncnorm(out$lambda[i],out$lambda[i+1],out$lambda[i-1],0,denom)
    }
    else {
     # pvals[i] = (pnorm(out$lambda[i-1]/denom)-pnorm(out$lambda[i]/denom))/
      #  (pnorm(out$lambda[i-1]/denom)-pnorm(0))
         pvals[i]=1-ptruncnorm(out$lambda[i],out$lambda[i+1],0,0,denom)
    }
  }
  return(pvals)
}
 plot.lar=
function (x, xvar = c("norm", "df", "arc.length", "step"), breaks = TRUE, 
    plottype = c("coefficients", "Cp"), omit.zeros = TRUE, eps = 1e-10, 
    ...) 
{
    object <- x
    plottype <- match.arg(plottype)
    xvar <- match.arg(xvar)
    coef1 <- object$beta
    if (x$type != "LASSO" && xvar == "norm") 
        coef1 = betabreaker(x)
    stepid = trunc(as.numeric(dimnames(coef1)[[1]]))
    coef1 <- scale(coef1, FALSE, 1/object$normx)
    if (omit.zeros) {
        c1 <- drop(rep(1, nrow(coef1)) %*% abs(coef1))
        nonzeros <- c1 > eps
        cnums <- seq(nonzeros)[nonzeros]
        coef1 <- coef1[, nonzeros, drop = FALSE]
    }
    else cnums <- seq(ncol(coef1))
    s1 <- switch(xvar, norm = {
        s1 <- apply(abs(coef1), 1, sum)
        s1/max(s1)
    }, df = object$df, arc.length = cumsum(c(0, object$arc.length)), 
        step = seq(nrow(coef1)) - 1)
    xname <- switch(xvar, norm = "|beta|/max|beta|", df = "Df", 
        arc.length = "Arc Length", step = "Step")
    if (plottype == "Cp") {
        Cp <- object$Cp
        plot(s1, Cp, type = "b", xlab = xname, main = object$type, 
            ...)
    }
    else {
        matplot(s1, coef1, xlab = xname, ..., type = "b", pch = "*", 
            ylab = "Standardized Coefficients")
        title(object$type, line = 2.5)
        abline(h = 0, lty = 3)
        axis(4, at = coef1[nrow(coef1), ], labels = paste(cnums), 
            cex = 0.8, adj = 0)
        if (breaks) {
            axis(3, at = s1, labels = paste(stepid), cex = 0.8)
            abline(v = s1)
        }
    }
    invisible()
}
