
# The following code runs 4 tests on the sample_truncnorm sampler (sample_truncnorm.R). 
#
# For each of these runs, we compare the sample-mean to the approximated sample mean 
# based on the moments function (mtmvtnorm) of the tmvtnorm package. Scenarios 1 through 4 use d=5,
# whereas scenario 5 checks that the code does not break for scalars. 
#
# The scenarios are: 
# 1) Identity sigma, mu = 0, only lower bound
# 2) Colored sigma, mu = 0, only lower bound
# 3) Colored sigma, mu = 0, lower bound x_2,x_3,x_4 upper bound x_1,x_5
# 4) Colored sigma, mu >= 0, lower bound x_2,x_3,x_4 upper bound x_1,x_5
# 5) sigma = 2, mu = 1, lower bound = 2  

library(tmvtnorm)

#### Scenario 1:
# IID reference with mu = 0
SigmaIid = diag(5)
constr = thresh2constraints(5, lower = 2)

sample_1 = sample_truncnorm(ndraw = 1000,A= constr$A,b = constr$b, Sigma = SigmaIid, initial = rep(2.1,5), 
                            eta = rep(1,5),thinning = 500, how_often = 100)

exp_mean_1 = mtmvnorm( sigma = SigmaIid, mean = rep(0,5), lower = rep(2,5))$tmean
rms_err_1 = sqrt(mean((colMeans(sample_1) - exp_mean_1)^2))
cat(rms_err_1)

stopifnot(rms_err_1< 0.03 )

#### Scenario 2:
# Colored reference with mu = 0
SigmaCol = diag(5)
for (k in 1:4) {SigmaCol[k,k+1] = 0.3; SigmaCol[k+1,k] = 0.3}
SigmaCol[1,2] = SigmaCol[1,2] = SigmaCol[1,2] =SigmaCol[1,2]

sample_2 = sample_truncnorm(ndraw = 1000,A= constr$A,b = constr$b, Sigma = SigmaCol, initial = rep(2.1,5), 
                            eta = rep(1,5),thinning = 100, how_often = 100)

exp_mean_2 = mtmvnorm( sigma = SigmaCol, mean = rep(0,5), lower = rep(2,5))$tmean
rms_err_2 = sqrt(mean((colMeans(sample_2) - exp_mean_2)^2))
cat(rms_err_2)

stopifnot(rms_err_2< 0.03 )



#### Scenario 3:
# Colored reference with tied-down bump
constr_tied = thresh2constraints(5, lower =c(-Inf, 2,2,2, -Inf), upper = c(2, Inf,Inf,Inf,2) )
sample_3 = sample_truncnorm(ndraw = 1500,A= constr_tied$A,b = constr_tied$b, Sigma = SigmaCol, initial = 2+ c(-0.1,0.1,0.1,0.1,-0.1), 
                            eta = c(0,1,1,1,0),thinning = 100, how_often = 100)

exp_mean_3 = mtmvnorm( sigma = SigmaCol, mean = rep(0,5), lower =c(-Inf, 2,2,2, -Inf), upper = c(2, Inf,Inf,Inf,2) )$tmean
rms_err_3 = sqrt(mean((colMeans(sample_3) - exp_mean_3)^2))
cat(rms_err_3)

stopifnot(rms_err_3< 0.03 )


#### Scenario 4:
# Colored reference with tied-down bump and non-zero mean
constr_tied = thresh2constraints(5, lower =c(-Inf, 2,2,2, -Inf), upper = c(2, Inf,Inf,Inf,2) )
mu_vec = c(0,1,1,0.5,0)
sample_4 = sample_truncnorm(ndraw = 1500,A= constr_tied$A,b = constr_tied$b, Sigma = SigmaCol, mu = mu_vec, initial = 2+ c(-0.1,0.1,0.1,0.1,-0.1), 
                            eta = c(0,1,1,1,0),thinning = 100, how_often = 100)

exp_mean_4 = mtmvnorm( sigma = SigmaCol, mean = mu_vec, lower =c(-Inf, 2,2,2, -Inf), upper = c(2, Inf,Inf,Inf,2) )$tmean
rms_err_4 = sqrt(mean((colMeans(sample_4) - exp_mean_4)^2))
cat(rms_err_4)

stopifnot(rms_err_4< 0.03 )


#### Scenario 5:
sigma = 2 # The variance
mu = 1
lower_bound = 2
constr_1d = thresh2constraints(1, lower = lower_bound)
sample_5 = sample_truncnorm(ndraw = 2500,A= constr_1d$A,b = constr_1d$b, Sigma = sigma, mu = mu, initial = 2.1, 
                            eta = 1,thinning = 10, how_often = 100)

exp_mean_5 = mtmvnorm( sigma = sigma, mean = mu, lower = lower_bound)$tmean
rms_err_5 = abs(mean(sample_5) - exp_mean_5)
cat(rms_err_5)

stopifnot(rms_err_5< 0.03 )

