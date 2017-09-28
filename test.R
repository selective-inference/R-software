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
    soln_R = rep(0, p)
    grad = -t(X) %*% Y / n
    ever_active = as.integer(c(1, rep(0, p-1)))
    nactive = as.integer(1)
    kkt_tol = 1.e-12
    objective_tol = 1.e-16
    maxiter = 500
    soln_R = selectiveInference:::solve_QP(t(X) %*% X / n, lam, maxiter, soln_R, -t(X) %*% Y / n, grad, ever_active, nactive, kkt_tol, objective_tol, p)$soln

    # test wide solver
    Xtheta = rep(0, n)
    nactive = as.integer(1)
    ever_active = as.integer(c(1, rep(0, p-1)))
    soln_R_wide = rep(0, p)
    grad = - t(X) %*% Y / n
    soln_R_wide = selectiveInference:::solve_QP_wide(X, lam, maxiter, soln_R_wide, -t(X) %*% Y / n, grad, Xtheta, ever_active, nactive, kkt_tol, objective_tol, p)$soln

#     soln_R = selectiveInference:::solve_QP(t(X) %*% X / n, lam, maxiter, soln_R, -t(X) %*% Y / n, grad, ever_active, nactive, kkt_tol, objective_tol, p)
#     print('active')
#     print(nactive)
#     print(ever_active)
#     print(soln_R$ever_active)
#     soln_R = soln_R$soln
#     soln_R_old = soln_R
#     print(soln_R)                                                                                                                                                                        
#     Xtheta = rep(0, n)                                                                                                                                                                   
#     nactive = as.integer(1)
#     ever_active = as.integer(c(1, rep(0, p-1)))
#     soln_R = rep(0, p)
#     grad = - t(X) %*% Y / n
#     # test wide solver                                                                                                                                                                   
#     soln_R_wide = selectiveInference:::solve_QP_wide(X, lam, maxiter, soln_R*1., -t(X) %*% Y / n, grad, Xtheta, ever_active, nactive, kkt_tol, objective_tol, p)
#     print(nactive)
#     print(soln_R_wide$ever_active)

     print('diff')
    print(soln_R_wide - soln_R)
#     print(soln_R_wide$gradient[soln_R_wide$ever_active])
#     print(max(abs(soln_R_wide$gradient[-soln_R_wide$ever_active])))
#     print(soln_R_wide$kkt_check)
#     print(soln_R_wide$iter)
#     print(max(abs(Xtheta - X %*% soln_R_wide$soln)))
# #     print(Xtheta)

#       print('R objective')
#       print(0.5 * sum(Xtheta^2)/n - sum(Xtheta*Y)/n + lam * 0.7 * sum(abs(soln_R_wide$soln)))
#       print(max(abs(soln_R_wide$gradient - t(X) %*% X %*% soln_R_wide$soln / n + t(X) %*% Y / n)))
#       print(lam)
#      print(max(abs(soln_R_wide$gradient[soln_R_wide$soln != 0])))
#      print(which(soln_R_wide$soln != 0))
#      soln_R_wide = selectiveInference:::solve_QP_wide(X, 0.7 * lam, maxiter, soln_R_wide$soln, -t(X) %*% Y / n, grad, Xtheta, ever_active, nactive, kkt_tol, objective_tol, p)
# #     print(Xtheta - X %*% soln_R_wide$soln)
# #     print(soln_R_wide$soln)
#      soln_R_wide = selectiveInference:::solve_QP_wide(X, 0.5 * lam, maxiter, soln_R_wide$soln, -t(X) %*% Y / n, grad, Xtheta, ever_active, nactive, kkt_tol, objective_tol, p)