library(selectiveInference)
set.seed(32103)

# parameters

n=20; p=12; nsim=50
sigma=1; alpha=.1; lam=0.04

# design

X=matrix(rnorm(n*p),n,p)
X=scale(X,TRUE,TRUE) / sqrt(n-1)

# truth

btrue=c(1,4,-4,rep(0,p-3)) 

sign=rep(NA,p)
ci=ci0=ci2=ci4=ci5=ci6=vlims0=array(NA,c(p,2,nsim))

betaall=matrix(NA,nsim,p)
btruepart=matrix(NA,p,nsim)
mc2=matrix(NA,nsim,p)
zalpha=-qnorm(alpha/2)


seeds=sample(1:99999, size=nsim)

for(ii in 1:nsim){
    set.seed(seeds[ii])
    mutrue = X%*%btrue
    y = mutrue+sigma*rnorm(n)
    y = y - mean(y)

    if (p>1) {   
        G = glmnet(X,y,standardize=FALSE)
        beta=as.numeric(coef(G, s=lam, exact=TRUE, x=X, y=y))[-1]
    } else if (p==1) {
        coef=lsfit(x,y)$coef[2]
        beta=sign(coef)*(abs(coef)-n*lam)*(abs(coef)>n*lam)
    }

    # any active
    if (sum(beta!=0) > 0){
        betaall[ii,]=beta
        a=lsfit(X[,beta!=0], y)
        aa=ls.diag(a)
        bhat=a$coef[-1]
        bhat0=a$coef[1]
        act=which(beta!=0)
        se=aa$std.err[-1]
        btruepart[,ii]=0
        btruepart[act,ii]=lsfit(X[,act,drop=F],mutrue)$coef[-1]

        #naive intervals
        ci0[beta!=0,1,ii]=bhat-zalpha*se
        ci0[beta!=0,2,ii]=bhat+zalpha*se

        #bonf-adj naive
        alpha4=alpha/p
        zalpha4=-qnorm(alpha4/2)
        ci4[beta!=0,1,ii]=bhat-zalpha4*se
        ci4[beta!=0,2,ii]=bhat+zalpha4*se

        #lee et al intervals
        lee = fixedLassoInf(X,
                            y,
                            beta,
                            lam*n,
			    alpha=alpha, 
			    family='gaussian',
                            type='partial', 
                            sigma=sigma)
        ci_lee = matrix(NA, p, 2)
	ci_lee[lee$vars,] = lee$ci    
        ci[,,ii] = ci_lee

        #randomized

        rand_lasso_soln = randomizedLasso(X, y, n*lam, noise_scale=0.5*sigma, ridge_term=0.5*sigma/sqrt(n))
        rand_E = rand_lasso_soln$active_set
        X_rand_E = X[,rand_E,drop=FALSE]
        rand_target = coef(glm(y~X_rand_E-1))
        cov_rand_E = solve(t(X_rand_E) %*% X_rand_E) * sigma^2
        cross_cov = matrix(0, p, length(rand_E))
        cross_cov[1:length(rand_E),1:length(rand_E)] = cov_rand_E
        targets = list(observed_target=rand_target,
                       cov_target=cov_rand_E, 
                       crosscov_target_internal=cross_cov)
        rand_inf = randomizedLassoInf(rand_lasso_soln, nsample=5000, burnin=1000, targets=targets)
        ci5[rand_lasso_soln$act,1,ii] = rand_inf$ci[,1]
        ci5[rand_lasso_soln$act,2,ii] = rand_inf$ci[,2]
    }




    mc0=mean(ci0[,1,1:ii]>btruepart[,1:ii] | ci0[,2,1:ii]<btruepart[,1:ii],na.rm=T)
    len0=mean(ci0[,2,1:ii]-ci0[,1,1:ii],na.rm=T)
    len0m=median(ci0[,2,1:ii]-ci0[,1,1:ii],na.rm=T)

    ninf=mean(abs(ci)==Inf,na.rm=T)
    ci[abs(ci)==Inf]=NA
    mc=mean(ci[,1,1:ii]>btruepart[,1:ii] | ci[,2,1:ii]<btruepart[,1:ii],na.rm=T)
    len=mean(ci[,2,1:ii]-ci[,1,1:ii],na.rm=T)
    lenm=median(ci[,2,1:ii]-ci[,1,1:ii],na.rm=T)


    mc4=mean(ci4[,1,1:ii]>btruepart[,1:ii] | ci4[,2,1:ii]<btruepart[,1:ii],na.rm=T)
    len4=mean(ci4[,2,1:ii]-ci4[,1,1:ii],na.rm=T)
    len4m=median(ci4[,2,1:ii]-ci4[,1,1:ii],na.rm=T)

    mc5=mean(ci5[,1,1:ii]>btruepart[,1:ii] | ci5[,2,1:ii]<btruepart[,1:ii],na.rm=T)
    len5=mean(ci5[,2,1:ii]-ci5[,1,1:ii],na.rm=T)
    len5m=median(ci5[,2,1:ii]-ci5[,1,1:ii],na.rm=T)



    cat(c("trial ", ii, " (cumulative)\n"))
    cat(c("ave# nz",mean(apply(betaall[1:ii,,drop=FALSE]!=0,1,sum,na.rm=T))),fill=T)

    miscov=c(mc0,mc,mc4,mc5)

    res=
    matrix(
        c(len0,len0m,0,
          len,lenm,ninf,
          len4,len4m,0,
            len5,len5m,0
          ),4,3,byrow=T)
    res=cbind(res,miscov)

    dimnames(res)=list(c("naive","Lee","bonf-naive","rand"),c("mean","median","propInf","miscov"))

    print(res)


}  
