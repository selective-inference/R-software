require(truncnorm)
require(MASS)
require(intervals)
require(glmnet)

mytruncnorm <- function(etay, vpos, vneg, sigma, etamu){
	# if range is too many sds away from mu, then there
	# will be numerical errors from using truncnorm
	if(max(vneg-etamu,etamu-vpos)/sigma < 7){
		return(ifelse(etay < etamu, 
		       ptruncnorm(etay, vneg, vpos, etamu, sigma),
		       1 - ptruncnorm(etay, vneg, vpos, etamu, sigma)))
	}
	else{
	    return(ifelse(etay < etamu,
	           1 - pexp(vpos-etay, etamu-vpos)/
	               pexp(vpos-vneg, etamu-vpos),
	           1 - pexp(etay-vneg, vneg-etamu)/
	               pexp(vpos-vneg, vneg-etamu)
	    ))	
    }
}

buildALasso <- function(x, y, beta, lam){
    n = length(y)
    tmp <- -t(x)%*%(x%*%beta - y)/lam
    E <- sapply(beta, all.equal, current=0)!=TRUE # instead of beta!=0
    P_E <- x[,E, drop=FALSE]%*%ginv(x[,E, drop=FALSE])
    xTxinv <- solve(crossprod(x[,E, drop=FALSE]))
    A0 <- rbind(t(x[,!E, drop=FALSE])%*%(diag(n) - P_E), 
                -t(x[,!E, drop=FALSE])%*%(diag(n) - P_E))/lam
    A1 <- -diag(tmp[E], nrow=sum(E))%*%xTxinv%*%t(x[,E, drop=FALSE])
    return(rbind(A0, A1))
}

buildbLasso <- function(x, y, beta, lam){
	tmp <- -t(x)%*%(x%*%beta - y)/lam
	E <- sapply(beta, all.equal, current=0)!=TRUE # instead of beta!=0
	mp_E <- ginv(t(x[,E, drop=FALSE]))
	xTxinv <- solve(crossprod(x[,E, drop=FALSE]))
	b0 <- c(rep(1, sum(!E)) - t(x[,!E, drop=FALSE])%*%mp_E%*%tmp[E], 
	        rep(1, sum(!E)) + t(x[,!E, drop=FALSE])%*%mp_E%*%tmp[E])
    b1 <- -lam*diag(tmp[E], nrow=sum(E))%*%xTxinv%*%tmp[E]
    return(c(b0, b1))
}


calcp <- function(A, b, eta, y, etamu=0){
	alpha = A%*%eta/sum(eta^2)
    Ay = A%*%y
    etay = sum(eta*y)
    vs = sapply(1:length(alpha), function(ii){
        (b[ii]-Ay[ii]+alpha[ii]*etay)/alpha[ii]})
    vneg = max(vs[alpha < 0])
    vpos = min(vs[alpha > 0])
    sigma = sqrt(sum(eta^2))
    mytruncnorm(etay, vpos, vneg, sigma, etamu)
}


# calculate the feasible intervals implied by selection event
calcintervals <- function(AmVnl, AmVdl, Amd, b, ct, trace=FALSE){
    allintervals = lapply(1:length(b), function(ii){
                                        #if(ii%%100 == 0){print(ii)}
        tmpf = functionfact(AmVnl[ii], AmVdl[ii], Amd[ii], b[ii])
        ta = AmVnl[ii]
        tb = Amd[ii] - b[ii]
        intervals = NULL
        if(sign(ta) != sign(tb) & abs(tb) > abs(ta))
            {root1 = tryCatch(uniroot(tmpf, interval= c(0, ta^2/(tb^2-ta^2)))$root, error = function(e){NULL})
             root2 = tryCatch(uniroot(tmpf, interval= c(ta^2/(tb^2-ta^2), ct+20))$root, error = function(e){NULL})
             if(!is.null(root1))    #convert a root into the interval where ct is allowed to fall
                 {s1 = sign(tmpf(0))
                  s2 = sign(tmpf(ta^2/(tb^2-ta^2)))
                  if(s1 > s2) intervals = rbind(intervals, c(root1, ifelse(is.null(root2),Inf,root2)))
                  if(s2 > s1) intervals = rbind(intervals, c(0, root1))
              }
             if(!is.null(root2))    #convert a root into the interval where ct is allowed to fall
                 {s1 = sign(tmpf(ta^2/(tb^2-ta^2)))
                  s2 = sign(tmpf(ct+20))
                  if(s1 > s2) intervals = rbind(intervals, c(root2, Inf))
                  if(s2 > s1) intervals = rbind(intervals, c(ifelse(is.null(root1),0,root1), root2))
              }
             return(intervals)
         }
        tryCatch({root1 = uniroot(tmpf, interval= c(0, ct+20))$root
                  intervals = Intervals(closed=c(F,F))
                  if(!is.null(root1))    #convert a root into the interval where ct is allowed to fall
                      {s1 = sign(tmpf(0))
                       s2 = sign(tmpf(ct + 20))
                       #cat(s1,s2)
                       if(s1 > s2) intervals = rbind(intervals, c(root1, Inf))
                       if(s2 > s1) intervals = rbind(intervals, c(0, root1))
                       return(intervals)
                   }
              }        
                ,error = function(e){NULL})
         })
    notnull = which(!sapply(allintervals, is.null))
    ints = lapply(allintervals[notnull], function(intmat){Intervals(intmat,closed =  c(TRUE,TRUE), type="R")})
    finalinterval <- do.call(interval_intersection,ints)
    if(trace == TRUE){print(finalinterval)}
    return(as.vector(finalinterval))
}

# build projection onto complement
Proj <- function(mat){diag(nrow(mat)) - mat%*%ginv(mat)}

# special equality test bc glmnet is weird
iseq <- function(x, vec){sapply(vec, all.equal, current=x)==TRUE}

# function to create functions whose roots imply the truncation
functionfact <- function(AmVnl, AmVdl, Amd, b){
  function(x){sqrt(x)*AmVnl + AmVdl + sqrt(1+x)*(Amd-b)}
}

# is Ay <= b?
checkOK <- function(x, y, beta, lam){
    A = buildALasso(x, y, beta, lam)
    b = buildbLasso(x, y, beta, lam)
    return(max(A%*%y - b) <= 0)
}


# here goes!  using 1/2RSS +lam*L1 formulation, so to convert to glmnet use lam/n in the glmnet call
# z = null tests using X_E against no model at all
# z = 1 tests against using ybar to predict
# z = z will test against using the columns in z (make sure to add an intercept if you want one!)
# does not work in the case where truncation is not a simple interval (theoretically possible, but I've never seen it).
prevallassop <- function(x, y, z, lam, trace=FALSE, debug=FALSE, fallbackthr = 1e-15){
    this.call=match.call()
    n <- length(y)
    tmpmod <- cv.glmnet(x, y, standardize=FALSE, intercept=FALSE)
    tmpmod <- cv.glmnet(x, y, standardize=FALSE, intercept=FALSE, lambda=sort(c(lam/n,tmpmod$lambda), decreasing=TRUE))
    beta <- as.vector(coef(tmpmod, s=lam/n)[-1])
    nz <- sum(!iseq(0, beta))
    if(nz >= n){stop("lam too small!")}
    if(nz == 0){stop("lam too big!")}
    A <- buildALasso(x, y, beta, lam)
    b <- buildbLasso(x, y, beta, lam)
    if(!checkOK(x,y,beta,lam)){
        if(trace){print(paste("Whoops!  Needed to refit, one second..."))}
        tmpmod <- cv.glmnet(x, y, standardize=FALSE, intercept=FALSE, lambda=tmpmod$lambda, thresh = fallbackthr)
        beta <- as.vector(coef(tmpmod, s=lam/n)[-1])
        A <- buildALasso(x, y, beta, lam)
        b <- buildbLasso(x, y, beta, lam)
        if(!checkOK(x,y,beta,lam)){stop("set an even smaller fallbackthr")}
    }
    nz <- sum(!iseq(0, beta))
    if(trace){print(paste("Fitting complete, number of nonzero coefs: ", nz, ".", sep=""))}
    if(is.null(z)){
        oPz <- diag(n)
    }
    else{
        if(isTRUE(all.equal(z,1))){
            z <-  rep(1,n)
        }
        oPz = Proj(z)
    }
    oPm = Proj(cbind(x[,!iseq(0, beta)],z))
    R1 = oPz%*%y
    R2 = oPm%*%y
    df1 = sum(svd(oPz)$d) - sum(svd(oPm)$d)
    df2 = sum(svd(oPm)$d)
    delta = y - R1
    l = sqrt(sum(R1^2))
    myc = df1/df2
    ct = (sum(R1^2)-sum(R2^2))/sum(R2^2)
    t = ct/myc
    Vn = (R1 - R2)/sqrt(sum((R1-R2)^2))
    Vd = R2/sqrt(sum((R2)^2))
    AmVnl = A %*% Vn * l
    Amd = A%*%delta
    AmVdl = A %*% Vd *l
    myint = as.vector(calcintervals(AmVnl, AmVdl, Amd, b, ct))
    if(debug == TRUE){browser()}
    if(length(myint)==0){
        Fmin <- 0;Fmax <- Inf}
    else{
        Fmin = min(myint)/myc
        Fmax = max(myint)/myc}
    #cat(Fmin, Fmax)
    if(trace){print(paste("Done calculating truncation.  (Fmin, F, Fmax) = (", Fmin,",", t,",", Fmax, ").", sep=""))}
    pval <- (pf(t,df1,df2, lower.tail=FALSE) - pf(Fmax, df1, df2, lower.tail = FALSE))/
               (pf(Fmin,df1,df2, lower.tail=FALSE) - pf(Fmax, df1, df2, lower.tail = FALSE))
    out <- list(pval = pval, nz = nz, beta = beta, truncation = c(Fmin, Fmax), F = t, call = this.call, lam = lam, df1 = df1, df2 = df2)
    class(out) <- "LassoFInf"
    return(out)
}

print.LassoFInf <-  function(x, digits = max(3, getOption("digits") - 3), ...){
    cat("\nCall: ", deparse(x$call), "\n")
    cat(c("lambda =", round(x$lam, digits)),fill=T)
    cat("",fill=T)
    cat(c("F =", round(x$F, digits),
                "and is constrained to the interval (",
                round(x$truncation[1], digits),",",
                round(x$truncation[2], digits),")."), fill=TRUE)
    cat(c("df1 =",round(x$df1),"and df2 =",round(x$df2),"degrees of freedom."),fill=T)
    cat("",fill=T)
    cat(c("P-value is", round(x$pval,digits)),fill=T)
}


#n = 50; p = 100; pz = 10; b = 0; bz  = sqrt(10);lam = 5
#set.seed(1)
#z <- matrix(rnorm(n*pz),n)
#z <- scale(z)
#x <- matrix(rnorm(n*(p)),n)
#x <- scale(x)
#eps = rnorm(n)
#y = bz*apply(z, 1, mean) + b*apply(x[,1:pz], 1, mean) + eps
#tmp <- prevallassop(x,y,z,lam)
#tmp                               
