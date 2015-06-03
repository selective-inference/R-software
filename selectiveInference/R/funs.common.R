
selection.int = function(y,eta,sigma,vs,alpha,del=1e-4,flip=F) {
  # compute selection intervals
    #Rob's version using grid search
    gridfac=50 
    etay=sum(eta*y)
     vm=vs$vm; vp=vs$vp
    fun = function(x,etay,vm,vp,sigma.eta)
       return(1-myptruncnorm(etay,a=vm,b=vp,mean=x,sd=sigma.eta))
  #  fun = function(x,etay,vm,vp,sigma.eta) return(1-ptruncnorm(etay,vm,vp,x,sigma.eta))
lo=-Inf
hi=Inf
covlo=covhi=0
     sigma.eta=sqrt(sum(eta^2))*sigma
#  cat(c(vm,etay,vp,sigma.eta),fill=T)
if( min(etay-vm,vp-etay)>0.001*sigma.eta){
  inc = sqrt(sum(eta^2)*sigma)*del
   
    xL=etay-gridfac*sigma.eta
    xR=etay+gridfac*sigma.eta
  lo = grid.search(fun,alpha/2,xL,xR,inc=inc,etay=etay,vm=vm,vp=vp,sigma.eta=sigma.eta)
  hi = grid.search(fun,1-alpha/2,xL,xR,inc=inc,etay=etay,vm=vm,vp=vp,sigma.eta=sigma.eta)
    covlo=fun(lo,etay,vm,vp,sigma.eta)
    covhi=1-fun(hi,etay,vm,vp,sigma.eta)
 
}
    if(flip){temp=hi;hi=-lo;lo=-temp;  temp2=covlo;covlo=covhi;covhi=temp2}
  return(list(ci=c(lo,hi), miscov=c(covlo,covhi)))
}

grid.search=function(fun, val, xL,xR, inc=0.01, tol=1e-2,etay=etay,vm=vm,vp=vp,sigma.eta=sigma.eta) {
#
G = 50000
  xg = seq(xL,xR,length=G)
 # vals = numeric(G)
#  for (i in 1:G) vals[i] = fun(xg[i],etay,vm,vp,sigma.eta)
vals = fun(xg,etay,vm,vp,sigma.eta)
o=vals<=val
v=abs(vals-val)
  return(xg[which.min(v)])

}


tnorm.cdf = function(x,mean,sd,a,b) {
  return((pnorm(x,mean,sd)-pnorm(a,mean,sd))/
         (pnorm(b,mean,sd)-pnorm(a,mean,sd)))
}

tnorm.surv = function(x,mean,sd,a,b) {
    #prob(X>x)
  return((pnorm(b,mean,sd)-pnorm(x,mean,sd))/
         (pnorm(b,mean,sd)-pnorm(a,mean,sd)))
}


#NOT CURRENTLY USED
bin.search = function(x, fun, val, inc=0.01, tol=1e-2) {  
  # Find left and right bounds
  xl = x
  while (fun(xl)>val) xl = xl-inc
  if (xl==x) {
    xr = x
    while (fun(xr)<val) xr = xr+inc
    xl = xr-inc
  }
  else xr = xl+inc
  
  # Grid search
  G = 100
  xg = seq(xl,xr,length=G)
  vals = numeric(G)
  for (i in 1:G) vals[i] = fun(xg[i])
  return(xg[which.min(abs(vals-val))])
 
}
rob.ptruncnorm=function(etay, vneg,vpos, etamu, sigma){
    	if(max(vneg-etamu,etamu-vpos)/sigma < 7){
		     return(ptruncnorm(etay, vneg, vpos, etamu, sigma))
                 }
		   
	else{
            if(etay > etamu) return(1)
            if(etay < etamu) return(0)
             if(etay== etamu) return(.5)
        }
    }


compute.vmvp=function(y,eta, A, b, pp,  SMALL = 1e-07){
alp = as.vector(A %*% eta/sum(eta^2))
    alp[abs(alp) < SMALL] = 0
    vp = rep(Inf, pp)
    vm = rep(-Inf, pp)
    temp = (b - (A %*% y) + alp * sum(eta * y))/alp
    vm[alp<0]=temp[alp<0]
    vp[alp>0]=temp[alp>0]
    vmm = max(vm, na.rm = T)
    vpp = min(vp, na.rm = T)
   return(list(vm = vmm, vp = vpp))

}

is.wholenumber <-
         function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol

checkargs=function(y,x=NULL,bhat=NULL,sigma=NULL,lambda=0,alpha=.1,tol=1e-5,nsteps=NULL
,larfit=NULL, fsfit=NULL, k=NULL,bh.q=NULL){
    # checks arguments for all user-callable functions.
  
 if(!is.null(x)){
    if(!is.matrix(x)) stop("x must be matrix")
    if(ncol(x)<2) stop("x must have at least 2 columns")
    if(length(y)!=nrow(x)) stop("Number of rows of x must equal length of y")
}
 
    if(!is.null(bhat)){
         if(length(bhat)!=ncol(x))  stop("bhat must have length equal to the number of cols of x")
     }
    if(!is.null(sigma)){
         if(sigma<=0) stop("sigma must be gt 0")
     }
    if(lambda<0) stop("lambda must be non-negative")
    if(alpha<=0 | alpha >=1) stop("alpha must be between 0 and 1")
      if(tol<=0) stop("tol must be gt 0")
    if(!is.null(nsteps)){
        if(nsteps<1 | nsteps> min(nrow(x),ncol(x))) stop("nsteps must be between 1 and min(nrow(x),ncol(x))")
        if(!is.wholenumber(nsteps)) stop("nsteps must be an integer")
    }
  
       if(!is.null(larfit)){
             if(class(larfit)!="lar") stop("larfit must be a lar object returned by lar (not lars)")
         }
    if(!is.null(fsfit)){
             if(class(fsfit)!="forwardStep") stop("fsfit must be an object returned by forwardStep")
         }
if(!is.null(bh.q)){
  if(bh.q<=0 | bh.q >=1) stop("bh.q must be between 0 and 1")
   }
   if(!is.null(k)){
        if(k<1 | k >length(y))  stop("k must be between 1 and length(y)")
    } 
}

estimateSigma=function(x,y){
    if(nrow(x)<10) stop("Number of observations must be at least 10 to run estimateSigma")
    cvfit=cv.glmnet(x,y,standardize=F)
    lamhat=cvfit$lambda.min
    fit=glmnet(x,y,standardize=F)
    yhat=predict(fit,x,s=lamhat)
    nz=sum(predict(fit,s=lamhat, type="coef")!=0)
    sigma=sqrt(sum((y-yhat)^2)/(length(y)-nz-1))
    return(list(sigmahat=sigma, df=nz))
}

myptruncnorm=function(z,a,b,mean=0,sd=1){
 # return Prob(Z<z | Z in [a,b])
 #"mean" can be a vector
#uses 
# A UNIFORM APPROXIMATION TO THE RIGHT NORMAL TAIL INTEGRAL
# W Bryc
#Applied Mathematics and Computation
#Volume 127, Issues 2–3, 15 April 2002, Pages 365–374
#https://math.uc.edu/~brycw/preprint/z-tail/z-tail.pdf
    
f=function(z){
(z^2+5.575192695*z+12.7743632)/(z^3*sqrt(2*pi)+ 14.38718147*z*z+31.53531977*z+2*12.77436324)
}
#first try standard  formula
a=(a-mean)/sd
b=(b-mean)/sd
z=(z-mean)/sd
term1=1-pnorm(a)
term2=1-pnorm(b)
term3=1-pnorm(z)
out=(term1-term3)/(term1-term2)
o=is.na(out)
#if any are NAs, apply modified method
if(sum(o)>0){
     zz=z[o]
     aa=a[o]
    bb=b[o]
   term1=exp(zz*zz)
    oo=aa>-Inf
   term1[oo]=f(aa[oo])*exp(-(aa[oo]^2-zz[oo]^2)/2)
   term2=0
   ooo=bb<Inf
  term2[ooo]=f(bb[ooo])*exp(-(bb[ooo]^2-zz[ooo]^2)/2)
   out[o]= (term1-f(zz))/ (term1-term2)
}
return(out)
}
