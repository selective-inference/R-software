 library(selectiveInference,lib.loc="mylib")
library(truncnorm)



mytruncnorm = function(etay, vneg,vpos, etamu, sigma){
    # From Sam Gross- uses exp approximation in extreme tails
	# if range is too many sds away from mu, then there
	# will be numerical errors from using truncnorm
	if(max(vneg-etamu,etamu-vpos)/sigma < 7){
		     return(ptruncnorm(etay, vneg, vpos, etamu, sigma))
                 }
		   
	else{
	   
	      return(1 - pexp(vpos-etay, etamu-vpos)/ pexp(vpos-vneg, etamu-vpos))
	        
          }
    }

  alpha=.1
sigma=1.199
sigma.eta=1.215
vm=.4454
vp=.5702
etay=.5066
del=1e-4
gridfac=50
    fun = function(x,etay,vm,vp,sigma.eta) return(1-ptruncnorm(etay,vm,vp,x,sigma.eta))
lo=-Inf
hi=Inf
covlo=covhi=0
if( min(etay-vm,vp-etay)>.05*sigma.eta){
    xL=etay-gridfac*sigma.eta
    xR=etay+gridfac*sigma.eta
  lo = grid.search(fun,alpha/2,xL,xR,etay=etay,vm=vm,vp=vp,sigma.eta=sigma.eta)
  hi = grid.search(fun,1-alpha/2,xL,xR,etay=etay,vm=vm,vp=vp,sigma.eta=sigma.eta)

    covlo=fun(lo,etay,vm,vp,sigma.eta)
    covhi=1-fun(hi,etay,vm,vp,sigma.eta)
}
 cat(c(lo,hi,covlo,covhi),fill=T)

x=w=seq(xL,xR,length=10000)

w=fun(x,etay,vm,vp,sigma.eta)

plot(x,w)
0.4454098 0.5066348 0.5702267 1.199096 1.215243 0.1 50
