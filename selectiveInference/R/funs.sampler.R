# A no-rejection MCMC algorithm Jelena and Amir have been working on

log_concave_sampler = function(negative_log_density, 
                               grad_negative_log_density, 
                               constraints,
                               observed,
                               nsamples,
			       burnin){
  #print(constraints)
  constraints = as.matrix(constraints)
  dim = nrow(constraints)
  
  get_poisson_process = function(state){
    pos = as.matrix(state$pos)
    velocity = as.matrix(state$velocity)
    neg_velocity = velocity<0
    pos_velocity = velocity>0
    tau_min = 0
    tau_max = 10
    if (sum(neg_velocity)>0){
      R = (-constraints[neg_velocity,1]+pos[neg_velocity])/(-velocity[neg_velocity])
      tau_max = min(tau_max, min(R))
      L = (-constraints[neg_velocity,2]+pos[neg_velocity])/(-velocity[neg_velocity])
      tau_min = max(tau_min, max(L))
    }
    if (sum(pos_velocity)>0){
      R = (constraints[pos_velocity,2]-pos[pos_velocity])/velocity[pos_velocity]
      tau_max = min(tau_max, min(R))
      L = (constraints[pos_velocity,1]-pos[pos_velocity])/velocity[pos_velocity]
      tau_min = max(tau_min, max(L))
    }
    
    f=function(t){as.numeric(t(velocity) %*% grad_negative_log_density(pos+velocity*t))}
    tau_star = tau_max
    if (f(tau_min)*f(tau_max)<0){
      tau_star = uniroot(f, c(tau_min, tau_max))$root
    } else{
      if (negative_log_density(pos+velocity*tau_min)<negative_log_density(pos+velocity*tau_max)){
        tau_star = tau_min
      }
    }
    
    tau_min = max(tau_min, tau_star)
    
    RHS = negative_log_density(pos+velocity*tau_star)+rexp(1)
    g = function(t){negative_log_density(pos+velocity*t)-RHS}
    if (g(tau_min)*g(tau_max)<0){
      tau = uniroot(g, c(tau_min, tau_max))$root
    } else{
      tau = tau_max
    }
    return (tau)
  }
  
  update_velocity = function(){
    Z=rnorm(dim)
    return(Z/sqrt(sum(Z^2)))
  }
  
  compute_next = function(state){
    bounce_time = get_poisson_process(state)/2
    #print(paste("bounce time", bounce_time))
    next_pos = state$pos+state$velocity*bounce_time
    next_velocity=update_velocity()
    return(list(pos=next_pos, velocity=next_velocity))
  }
  
  state = list(pos=observed, velocity = update_velocity())
  samples = matrix(0, nrow = nsamples - burnin, ncol = dim)
  for (i in 1:nsamples){
    #print(paste("pos", toString(state$pos)))
    #print(paste("velocity", toString(state$velocity)))
    if (i > burnin) {
        samples[i - burnin,]=state$pos
    }
    state = compute_next(state)
  }
  return (samples)
}

gaussian_sampler = function(noise_scale, 
                            observed, 
                            linear_term, 
                            offset_term, 
                            constraints,
                            nsamples=10000,
			    burnin=2000){
  
  negative_log_density = function(x) {
    recon = linear_term %*% x+offset_term
    return(as.numeric(t(recon)%*%recon/(2*noise_scale^2)))
  }
  grad_negative_log_density=function(x){
    recon = linear_term %*% x+offset_term
    return(t(linear_term) %*% recon/(noise_scale^2))
  }
  
  return(log_concave_sampler(negative_log_density, 
                             grad_negative_log_density,
                             constraints,
                             observed,
                             nsamples,
			     burnin))
}
