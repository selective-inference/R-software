set.seed(43)

    n = 100                                                                                                                                                                              
    p = 200                                                                                                                                                                              
    lam = 0.2                                                                                                                                                                             
    X = matrix(rnorm(n*p), n, p)                                                                                                                                                         
    Y = rnorm(n)                                                                                                                                                                         
    library(selectiveInference)                                                                                                                                                          
    p = ncol(X)                                                                                                                                                                          
    soln_R = rep(0, p)                                                                                                                                                                   
    grad = -t(X) %*% Y / n                                                                                                                                                                   
    ever_active = as.integer(c(1, rep(0, p-1)))
    nactive = as.integer(1)                                                                                                                                                              
    kkt_tol = 1.e-12                                                                                                                                                                     
    objective_tol = 1.e-16                                                                                                                                                               
    maxiter = 500                                                                                                                                                                        
    soln_R = selectiveInference:::solve_QP(t(X) %*% X / n, lam, maxiter, soln_R, -t(X) %*% Y / n, grad, ever_active, nactive, kkt_tol, objective_tol, p)
    print('active')
    print(nactive)
    print(ever_active)
    print(soln_R$ever_active)
    soln_R = soln_R$soln
    soln_R_old = soln_R
    print(soln_R)                                                                                                                                                                        
    Xtheta = rep(0, n)                                                                                                                                                                   
    nactive = as.integer(1)
    ever_active = as.integer(c(1, rep(0, p-1)))
    soln_R = rep(0, p)
    grad = - t(X) %*% Y / n
    # test wide solver                                                                                                                                                                   
    soln_R_wide = selectiveInference:::solve_QP_wide(X, lam, maxiter, soln_R*1., -t(X) %*% Y / n, grad, Xtheta, ever_active, nactive, kkt_tol, objective_tol, p)
    print(nactive)
    print(soln_R_wide$ever_active)

    print('diff')
    print(soln_R_wide$soln - soln_R_old)
    print(soln_R_wide$gradient[soln_R_wide$ever_active])
    print(max(abs(soln_R_wide$gradient[-soln_R_wide$ever_active])))
    print(soln_R_wide$kkt_check)
    print(soln_R_wide$iter)
#    print(Xtheta - X %*% soln_R_wide$soln)
#     print(Xtheta)

     soln_R_wide = selectiveInference:::solve_QP_wide(X, 0.7 * lam, maxiter, soln_R_wide$soln, -t(X) %*% Y / n, grad, Xtheta, ever_active, nactive, kkt_tol, objective_tol, p)
#     print(Xtheta - X %*% soln_R_wide$soln)
#     print(soln_R_wide$soln)
     soln_R_wide = selectiveInference:::solve_QP_wide(X, 0.5 * lam, maxiter, soln_R_wide$soln, -t(X) %*% Y / n, grad, Xtheta, ever_active, nactive, kkt_tol, objective_tol, p)