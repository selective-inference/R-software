
test_log_concave_sampler = function(){
  samples = log_concave_sampler(negative_log_density= function(x){x^2/2}, 
                                grad_negative_log_density=function(x){x},
                                constraints = t(as.matrix(c(2,3))),
                                observed = 2, nsamples=10000)
  mean(samples)
  hist(samples)
}


test_gaussian_sampler =function(){
  samples = gaussian_sampler(1, 1, 1, 0,10000)
  mean(samples)
  hist(samples)
}
