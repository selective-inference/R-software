set.seed(43)

    n = 100                                                                                                                                                                              
    p = 200                                                                                                                                                                              
    lam = 10                                                                                                                                                                             
    X = matrix(rnorm(n*p), n, p)                                                                                                                                                         
    Y = rnorm(n)                                                                                                                                                                         
    library(selectiveInference)                                                                                                                                                          
    p = ncol(X)                                                                                                                                                                          
    soln_R = rep(0, p)                                                                                                                                                                   
    grad = -t(X) %*% Y                                                                                                                                                                   
    ever_active = c(1, rep(0, p-1))                                                                                                                                                      
    nactive = as.integer(1)                                                                                                                                                              
    kkt_tol = 1.e-12                                                                                                                                                                     
    objective_tol = 1.e-12                                                                                                                                                               
    maxiter = 500                                                                                                                                                                        
    Xtheta = rep(0, n)                                                                                                                                                                   
    soln_R = selectiveInference:::solve_QP(t(X) %*% X, lam, maxiter, soln_R, -t(X) %*% Y, grad, ever_active, nactive, kkt_tol, objective_tol, p)$soln                                    
    print(soln_R)                                                                                                                                                                        
    # test wide solver                                                                                                                                                                   
                                                                                                                                                                                       
    soln_R_wide = selectiveInference:::solve_QP_wide(X, lam, maxiter, soln_R, -t(X) %*% Y, grad, Xtheta, ever_active, nactive, kkt_tol, objective_tol, p)
    print(soln_R_wide)