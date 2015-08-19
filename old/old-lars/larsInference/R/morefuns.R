
# code from Jon Taylor
lars_inference = function(X, Y, GL, num_strong_signals=0, sigma=NULL, mu=NULL, debugging=FALSE, ls_step=4, num_resid_contrasts=80) {

    # mu should be an n-vector, i.e. length(mu) == nrow(Z)
    LF = lars_functionals(X, Y, GL, 
                          ls_step=ls_step,
                          mu=mu)

    theta = 1. / LF$weights
    theta_exact = 1 * theta
    theta_cumsum = sqrt(cumsum(theta[(num_strong_signals+1):length(theta)]^2))
    theta[(num_strong_signals+1):length(theta)] = theta_cumsum
    nstep = length(LF$lambda)

    pvalues = c()
    tails = c()
    Z = c()
    critical_values = c(Inf, LF$lambda, 0)
    
    if (is.null(mu)) {
       mu = rep(0, nrow(X))
    }

    Delta = LF$gamma %*% mu

    if (is.null(sigma)) {
        # if n > p, also use some residuals
	# this is rather inefficient here
	
	if (nrow(X) > ncol(X)) {
	    # svd method -- often seem to get LAPACK errors
            # maybe just use average of MAD and sigma^2
            svdX = svd(X)
            Presid = diag(1, nrow(X)) - svdX$u %*% t(svdX$u)
            ncomp = nrow(X) - ncol(X)
            Uresid = svd(Presid, nu=max(ncomp, num_resid_contrasts),
                         nv=max(ncomp, num_resid_contrasts))$u
            Zresid = t(Uresid) %*% Y
	}
        else {
            Zresid = c()
        } 

        Z = (LF$lambda - Delta) * theta_exact
        sigma = c()
	for (i in 1:10) {
	    rZ = (2*rbinom(length(Z),1,0.5)-1) * Z
            sigma = c(sigma,mad(c(Zresid,rZ)))
        }
        sigma = mean(sigma)
        Z = Z / sigma
    }
    else {
        Z = (LF$lambda - Delta) * theta_exact / sigma
    } 

    for (i in 2:nstep) {
        Vm = critical_values[i-1]
        L = critical_values[i]
        Vp = critical_values[i+1]
        cumsum_w = sigma / theta[i-1]
	delta = sum(LF$gamma[i-1,]*mu)

        exact_w = sigma / theta_exact[i-1]
	pm = pnorm((Vm-delta)/exact_w, lower.tail=FALSE, log.p=TRUE)
	pL = pnorm((L-delta)/exact_w, lower.tail=FALSE, log.p=TRUE)
	pp = pnorm((Vp-delta)/exact_w, lower.tail=FALSE, log.p=TRUE)
        pval_exact = (exp(pm-pL) - 1) / (exp(pm-pL) - exp(pp-pL))
	tails = c(tails, 2*exp(pL))

        pval_covtest_corrected = exp(-L*(L-Vp)/cumsum_w^2) 
	pval_covtest_uncorrected = exp(-L*(L-Vp)/exact_w^2)
        pvalues = rbind(pvalues, 
                        c(pval_exact, 
                          pval_covtest_corrected,
                          pval_covtest_uncorrected,
                          cumsum_w^2)) # corrected exponential scaling
    }		     

    bounds = LF$bounds
    bounds[,4] = bounds[,4] * sigma

    check = LF$check
    slack = LF$slack

    return(list(pvalues=pvalues,
                pivots=LF$pivots,
                bounds=bounds,
		check=check,
		weights=LF$weights,
		gamma=LF$gamma,
		tails=tails,
		Z=Z,
		sigma=sigma,
		actions=LF$actions,
		lambda=LF$lambda))
}



lars_functionals = function(X, Y, GL, ls_step=4,
                            debugging=FALSE,
                            mu=NULL) {
    # Parameters
    # ==========

    # X : design matrix used in type="lar" fit
    # Y : response used in type="lar" fit
    # output from genlasso
    # ls_step : compute linear functionals for least 
    #           squares coefficients of which step of type="lar" fit?
    # object returned by genlasso

    # Returns
    # =======
    # weights : standard deviation of linear functionals compute lambda values
    # gamma : matrix of functionals that should compute each lambda value
    # LS : matrix that forms the least squares solution based on variables 
    #      specified
    #      by ls_step
    # bounds : information used to form selection interval
    #      
    # pivots : pvalue for each least squares coefficient of model ls_step     

    weights = c()

    #GL = genlasso(as.numeric(Y),
    #              X,
    #              diag(1,ncol(X)),
    #              maxsteps=min(nrow(X),ncol(X),ls_step+1),
    #              approx=TRUE)
    FNAL = -diff(GL$gamma)
    weights = sqrt(diag(GL$gamma %*% t(GL$gamma)))

    active_set = GL$actions[,'variable'][1:ls_step]
    active_X = X[,active_set]
    LSeta = solve(t(active_X) %*% active_X) %*% t(active_X)

    # alternate computation for weights and functionals

    bounds = c()
    pivots = c()

    for (j in 1:nrow(LSeta)) {

        constraint_info = interval_constraints(FNAL,
                                               rep(0, nrow(FNAL)),
                                               Y,
                                               LSeta[j,],
                                               debugging=debugging,
                                               mu=mu)

        pivots = c(pivots, constraint_info$pval)
        bounds = rbind(bounds, 
                       c(constraint_info$lower,
                         constraint_info$center,
                         constraint_info$upper,
                         constraint_info$sigma))

        check = constraint_info$check
        slack = constraint_info$slack

        if (debugging==TRUE) {
            print('slack')
            print(slack)
        }
    }

    return(list(weights=weights,
		gamma=GL$gamma,
                bounds=bounds,
                check=check,
                pivots=pivots,
                lambda=GL$lambda,
		actions=active_set))
}



selint=function(larsfit,j,k,alpha,sigma,trace=F,del=alpha/5,maxit=5000,larsinf.obj=NULL){
# compute lower and upper selection intervals of coverage 1-alpha, at kth lars step for jth predictor chosen
# larsfit is the result of call to lars
#   if method can't get tail area below alpha/2, corresponding point is set to - or + Inf
# del is stepsize for sequential serach  
# returns selection interval and area1,area2- achieved tails areas on left and right
#  flag1 and flag2 are error flags for lower and upper point: 0=no error, 1=error in computing tail area
#         in latter case we would set the endpoint to +-infty

 act=unlist(larsfit$act)
 if(is.null(larsinf.obj)){larsinf.obj=lars_inference(x,y,sigma=sigma,ls_step=k)}
 vp=larsinf.obj$bounds[,1]
 betahat=larsinf.obj$bounds[,2]
 vm=larsinf.obj$bounds[,3]
 sigma.eta=larsinf.obj$bounds[,4]

junk2=selint0(vp[j],vm[j],betahat[j],sigma.eta[j],alpha,2,trace=trace,maxit=maxit)
junk1=selint0(vp[j],vm[j],betahat[j],sigma.eta[j],alpha,1,trace=trace,maxit=maxit)

cup=junk2$sol;
clow=junk1$sol;


junk11=areaf(clow,vp[j],vm[j],betahat[j],sigma.eta[j])
area1=junk11$area;flag1=junk11$flag
junk22=areaf(cup,vp[j],vm[j],betahat[j],sigma.eta[j])
area2=1-junk22$area;flag2=junk22$flag


clow[area1>(alpha/2)+del]=-Inf
cup[area2>(alpha/2)+del]=Inf
pv=g(0,vp[j],vm[j],betahat[j],sigma.eta[j])
pv=2*min(pv,1-pv)
return(list(clow=clow,cup=cup,area1=area1,area2=area2,flag1=flag1,flag2=flag2,act=unlist(larsfit$act),vp=vp,vm=vm,betahat=betahat,sigma.eta=sigma.eta,pv=pv))
}
                                


selint0=function(vp,vm,betahat,sigma.eta,alpha,dir,maxit=NULL,del=.01,trace=FALSE){
# sequential search method
eps=1e-6
s=sigma.eta
alp2=alpha/2
ii=0
 go=T
 a=betahat
a1=pnorm((vm-a)/s)
a2=pnorm((betahat-a)/s)
a3=pnorm((vp-a)/s)
if (dir==1) val=a1*(1-alp2)-a2+alp2*a3
if(dir==2) val=a1*alp2-a2+(1-alp2)*a3

while(go & ii<maxit){
   ii=ii+1
   go=F
   aold=a
   if(dir==1) a=a-del
   if(dir==2) a=a+del
   a1=pnorm((vm-a)/s)
   a2=pnorm((betahat-a)/s)
   a3=pnorm((vp-a)/s)
   valold=val
   if(dir==1) val=a1*(1-alp2)-a2+alp2*a3
   if(dir==2) val=a1*alp2-a2+(1-alp2)*a3

   if (dir==1 & val >0) go=T
   if (dir==2 & val <0) go=T
   if(trace) cat(c(dir,a,val),fill=T)
}

#linearly interpolate backwards
slope=(val-valold)/(a-aold)
if(abs(slope)>eps){
  a=aold-valold/slope
}
area=areaf(a,vp,vm,betahat,sigma.eta)
return(list(sol=a,area=area,niter=ii))
}


areaf=function(a,vp,vm,betaval,sigma.eta){
# Compute Prob_a(betahat*>betaval| betahat^* in [vp,vm])

eps=1e-7
s=sigma.eta
num=(pnorm((vm-a)/s)- pnorm((betaval-a)/s))
den=(pnorm((vm-a)/s)- pnorm((vp-a)/s))
val=num/den
flag=0
if(is.na(den)){flag==1}
if(!is.na(den)){
if(abs(den)<eps) flag=1
}
return( list(area=val,flag=flag))
}

 
# the rest of these functions are not currently used

foo2=function(vp,vm,betahat,sigma.eta,alpha,dir,maxit=NULL,del=.01,rang=2,trace=FALSE){
# binary seaech approach
ii=0
if(dir==1) {low=betahat-rang;up=betahat}
if(dir==2) {low=betahat;up=betahat+rang}
go=TRUE
ii=0
while(go & ii< maxit){
ii=ii+1
go=F
f1=hh(low,vp,vm,betahat,sigma.eta,alpha,dir)
f3=hh(up,vp,vm,betahat,sigma.eta,alpha,dir)
if(sign(f1)==sign(f3) & dir==1 ){rang=rang*1.5;low=betahat-rang;go=T}
if(sign(f1)==sign(f3) & dir==2 ){rang=rang*1.5;up=betahat+rang;go=T}
}
if(ii==maxit) cat("could not bracket sol",fill=T)

err=FALSE
  go=T
  while(go & ii<maxit){
 ii=ii+1
  go=F
a=(low+up)/2
f1=hh(low,vp,vm,betahat,sigma.eta,alpha,dir)
 f2=hh(a,vp,vm,betahat,sigma.eta,alpha,dir)
f3=hh(up,vp,vm,betahat,sigma.eta,alpha,dir)
if(sign(f1)==sign(f3)) err=TRUE
if(sign(f2)==sign(f1)){low=a}
if(sign(f2)==sign(f3)){up=a}
dif=abs(up-low)
if(abs(dif)>del) go=T
if(err) go=F
cat(c(low,a,up,f1,f2,f3),fill=T)
}
area=areaf(a,vp,vm,betahat,sigma.eta)
return(list(sol=a,area=area,niter=ii,err=err))
}

hh=function(a,vp,vm,betahat,sigma.eta,alpha,dir){
alp2=alpha/2
s=sigma.eta
a1=pnorm((vm-a)/s)
a2=pnorm((betahat-a)/s)
a3=pnorm((vp-a)/s)
if(dir==1) val=a1*(1-alp2)-a2+alp2*a3
if(dir==2) val=a1*alp2-a2+(1-alp2)*a3
return(val)
}

g2=function(a,vp,vm,betahat,sigma.eta,alpha,dir){
BIG=10e9
eps=1e-8
  s=sigma.eta
a1=(pnorm((vm-a)/s)- pnorm((betahat-a)/s))
a2= (pnorm((betahat-a)/s)-pnorm((vp-a)/s))
if (dir==1) val=abs(a1*(1-alpha)-a2*alpha) +BIG*(a1/(a1+a2) > alpha)+BIG*((abs(a1)<eps) + (abs(a2)<eps))
if (dir==2) val=abs(a2*(1-alpha)-a1*(alpha)) +BIG*(a2/(a1+a2) > alpha)+BIG*((abs(a1)<eps) + (abs(a2)<eps))
 return( val)
}

ff=function(a,vp,vm,betahat,sigma.eta,alpha,dir){
s=sigma.eta
BIG=10e8
a1=pnorm((vm-a)/s)
a2=pnorm((betahat-a)/s)
a3=pnorm((vp-a)/s)
if(dir==1) val=a1*(1-alp2)-a2+alp2*a3
if(dir==2) val=a1*(alp2)-a2+(1-alp2)*a3
return(val)
}

g=function(a,vp,vm,betahat,sigma.eta){
# computes pivot for beta
# this is Prob_a(betahat*>betahat)
eps=1e-9
        s=sigma.eta
        z1=betahat
 num=(pnorm((vm-a)/s)- pnorm((z1-a)/s))
 den=(pnorm((vm-a)/s)- pnorm((vp-a)/s))
 val=num/den
 if( (abs(num)<eps) & (abs(den)<eps)) val=0

 return( val)
}

simulate_from_constraints = function(A, b) {
    while (TRUE) {
    Y = rnorm(ncol(A))
    V = as.numeric(A%*%Y)+b
    if (all(V > 0)) {
       break
       }
    }
    return(Y)
}

interval_constraints = function(support_directions, 
                                support_offsets,
                                observed_data, 
                                direction_of_interest,
                                covariance=NULL,
                                tol = 1.e-4,
				mu = NULL,
                                debugging = FALSE,
				scale=1) {
    
#     Given an affine cone constraint $Ax+b \geq 0$ (elementwise)
#     specified with $A$ as `support_directions` and $b$ as
#     `support_offset`, a new direction of interest $w$, and
#     an observed Gaussian vector $X$ with some `covariance`, this
#     function returns $w^TX$ as well as an interval
#     bounding this value. 

#     The interval constructed is such that the endpoints are 
#     independent of $w^TX$, hence the $p$-value
#     of `Kac-Rice <http://arxiv.org/abs/1308.3020>`_
#     can be used to form an exact pivot.

#     `scale` is used as a multiplicative factor in the standard deviation

    # shorthand
    A = support_directions
    b = as.numeric(support_offsets)
    S = covariance
    X = observed_data
    w = as.numeric(direction_of_interest)

    U = as.numeric(A %*% X)
    U = U + b

    check=TRUE
    if (any(U <= -tol * max(abs(U)))) {
        check=FALSE
        warning('constraints not satisfied')
	if (debugging==TRUE) {
            print('constraints not satisfied')
	    print(U)
        }
    }

    if (!is.null(S)) {
        Sw = S %*% w
    }
    else {
        Sw = w
    }
    sigma = sqrt(sum(w*Sw))
    C = as.numeric(A %*% Sw / sigma^2)
    V = sum(w*X)
    if (is.null(mu)) {
    mu = rep(0, length(w))
    }
    delta = sum(w*mu)
    RHS = - as.numeric(A %*% X)
    RHS = (RHS + V * C - b) / C

    pos_coords = C > tol * max(abs(C))
    if (any(pos_coords)) {
        lower_bound = max(RHS[pos_coords])
    }
    else {
        lower_bound = -Inf
    }
    neg_coords = C < -tol * max(abs(C))
    if (any(neg_coords)) {
        upper_bound = min(RHS[neg_coords])
    }
    else {
        upper_bound = Inf
    }

    pvalue = compute_pvalue(lower_bound - delta, 
                            V - delta, 
                            upper_bound - delta, sigma*scale)


    return(list(lower=lower_bound-delta,center=V-delta,upper=upper_bound-delta,sigma=sigma,pvalue=pvalue,check=check,slack=U))
}

lars_fdr = function(lars_inference_obj,
                    sequential=FALSE,
		    q=0.05) {

    # shorthand
    LI = lars_inference_obj

    actions = as.numeric(LI$actions)
    qvalues = p.adjust(LI$tails, 'BH')

    if (sequential == FALSE) {
        # only use the tails declared significant by BH
        active = LI$actions[qvalues < q]
        inactive = LI$actions[qvalues >= q]

        active_steps = (1:length(qvalues))[qvalues < q]
        inactive_steps = (1:length(qvalues))[-active_steps]

        R = length(active_steps)
    }
    else {

        # take the last time BH called a tail significant
        active0 = LI$actions[qvalues < q]
        active_steps0 = (1:length(qvalues))[qvalues < q]
        R = max(active_steps0)
        
        active_steps = 1:R
        inactive_steps = (R+1):length(qvalues)

        active = LI$actions[active_steps]
        inactive = LI$actions[inactive_steps]

    }

    return(list(qvalues=qvalues,
                R=R,
                active=active,
		active_steps=active_steps,
		inactive=inactive,
		inactive_steps=inactive_steps))
}

compute_pvalue = function(lower_bound, Z, upper_bound, sigma) {
    pm = pnorm(upper_bound/sigma, lower.tail=FALSE, log.p=TRUE)
    pL = pnorm(Z/sigma, lower.tail=FALSE, log.p=TRUE)
    pp = pnorm(lower_bound/sigma, lower.tail=FALSE, log.p=TRUE)
    pvalue = (exp(pm-pL) - 1) / (exp(pm-pL) - exp(pp-pL))
    return(pvalue)
}
