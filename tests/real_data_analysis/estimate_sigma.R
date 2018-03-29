library(glmnet)
library(selectiveInference)

setwd("/Users/Jelena/Dropbox/kevin/jelena/real_data")
data=readRDS("real_data.rds")
X=data$X
y=data$y

set.seed(1)

CV = cv.glmnet(X, y, standardize=FALSE, intercept=FALSE, family="gaussian")
sigma_est=selectiveInference:::estimate_sigma(X,y,coef(CV, s="lambda.min")[-1]) # sigma via Reid et al.
print(c("sigma est via Reid et al.", sigma_est))  "sigma est via Reid et al." "0.620021797352246" 

print(c("std error of y", sqrt(var(y))))  "std error of y"    "0.488729008112197"

n=nrow(X)
m=floor(n/2)
subsample = sample(1:n,m, replace=FALSE)
leftover = setdiff(1:n, subsample)
CV = cv.glmnet(X[subsample,], y[subsample], standardize=FALSE, intercept=FALSE, family="gaussian")
beta_hat = coef(CV, s="lambda.min")[-1]
selected = which(beta_hat!=0)
LM = lm(y[leftover]~X[leftover,][,selected])
sigma_est = sigma(LM)
print(c("sigma est from data splitting", sigma_est)) "0.481421791225578"  




