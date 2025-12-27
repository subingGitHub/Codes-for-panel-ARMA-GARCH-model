library(Rsolnp)
library(VGAM)


beta<-3

phi<-0.3

psi<-0.3

tau<-0.2

nu<-0.5


n<-N<-100

T<-100

sample_size <- 500  

l_N<-matrix(rep(1,N))

l_T<-l<-matrix(rep(1,T))

J_N<-diag(N) - l_N%*%t(l_N)/N

J_T<-diag(T) - l_T%*%t(l_T)/T


lag<-4



our_method <- our_method_jacknife <- matrix(c(beta, phi, psi, tau, nu),nrow=1)


proflie_method <- proflie_method_jacknife <- matrix(c(beta, phi, psi, tau, nu),nrow=1)
 



#########################################
#########################################
#####
fun_A<-function(phi,T_used = T){
  
  A <- diag(T_used)
  diag(A[-1,-T_used]) <- -phi
  
  
  return(A)
}
#####

fun_A_dot1<-function(phi){
  
  A_dot1 <- diag(T)*0
  diag(A_dot1[-1,-T]) <- -1
  return(A_dot1)
}
######
######
fun_B<-function(psi,T_used = T){
  
  B <- diag(T_used)
  diag(B[-1,-T_used]) <- psi
  
  
  return(B)
}
########
fun_B_dot1<-function(psi){
  
  B_dot1 <- diag(T)*0
  diag(B_dot1[-1,-T]) <- 1
  return(B_dot1)
}
###############
fun_Sigma<-function(psi,T_used = T){
  B<-fun_B(psi,T_used)
  return(B%*%t(B))
}



#########################################
#########################################

omega<-matrix(runif(N,1,3))
mu<-rnorm(N)

while(dim(our_method)[1]<sample_size){ 
  
  

  U<-c()
  
  for(i in 1:N){
    eps0<-0
    h0<-omega[i]
    eps<-c() 
    
    for(t in 1:(2*T)){
      h0<- omega[i]*(1-nu-tau) + nu*h0 + tau*eps0^2  
      
      eps0<-sqrt(h0)*rnorm(1)
      #eps0<-sqrt(h0)*sqrt(0.6)*rt(1,5)
      #eps0<-sqrt(h0)*runif(1,-sqrt(3),sqrt(3))
      
      eps<-rbind(eps,eps0)
    }
    
    U<-rbind(U,eps)
    
  }
  
  
  
  X_all<-matrix(rnorm(N*T*2))
  
  
  
  
  A_phi<-fun_A(phi,2*T)
  B_psi<-fun_B(psi,2*T)
  Y_all<-(diag(N)%x%solve(A_phi))%*%( X_all*beta + mu%x%matrix(rep(1,2*T)) + (diag(N)%x%B_psi)%*%U )
  
  Y<-c()
  X<-c()
  for(i in 1:N){
    Y<-rbind(Y,matrix(Y_all[((i-1)*2*T+1):((i-1)*2*T+T)]))
    X<-rbind(X,matrix(X_all[((i-1)*2*T+1):((i-1)*2*T+T)]))
    
  }
  
  
  #########################################
  #########################################
  
  
  fun_OLS <- function(par){
    
    beta<-par[1]
    phi<-par[2]
    psi<-par[3]
    
    
    A<-fun_A(phi)
    
    Sigma<-fun_Sigma(psi)
    
    Sigma_<-solve(Sigma)
    
    
    l <- rep(1,T)
    
    C<-Sigma_%*%l%*%t(l)%*%Sigma_/c(t(l)%*%Sigma_%*%l)
    
    V<-(diag(N)%x%A)%*%Y - X*beta
    
    value<-0
    for(i in 1:N){
      
      y_i<-matrix(Y[((i-1)*T+1):(i*T)])
      
      x_i<-matrix(X[((i-1)*T+1):(i*T)])
      v_i<-A%*%y_i - x_i*beta
      
      value<- value + t(v_i)%*%(Sigma_ - C)%*%v_i
    }
    
    return(value)   
    
  }
  
  
  estml<-try(suppressWarnings(solnp(c(beta,phi,psi), fun = fun_OLS,
                                    # ineqfun = fun_ineq, 
                                    # ineqLB = c(-0.999), 
                                    # ineqUB = c(0.999),
                                    LB = c(-Inf,-0.999,-0.999), 
                                    UB = c(Inf,0.999,0.999),
                                    control = list(trace = 0)))
             ,silent = TRUE)
  
  if('try-error' %in% class(estml)){next} 
  
   
  est_gamma<-estml$pars
  
  est_A<-fun_A(est_gamma[2])
  
  est_Sigma<-fun_Sigma(est_gamma[3])
  
  est_Sigma_<-solve(est_Sigma)
  
  c_Sigma<-c(t(l)%*%est_Sigma_%*%l)
  
  est_V<-(diag(N)%x%est_A)%*%Y - X*est_gamma[1]
  
  est_mu<-(diag(N)%x%(1/c_Sigma*t(l_T)%*%est_Sigma_))%*%est_V
  
  ############################################
  ###########################################
  
  est_B<-fun_B(est_gamma[3])
  
  est_U <- (diag(N)%x%solve(est_B)) %*% (est_V - est_mu%x%l_T)
  #
  ##########################
  est_omega<-c()
 
  
  for(i in 1:N){
    est_omega<-c(est_omega,mean((est_U[((i-1)*T+1):(i*T)])^2))
 
    
  }
  
    
  ml_GARCH<-function(par){
    
    tau<-par[1]
    nu<-par[2]
    
    
    value<-0
    for(i in 1:N){
      h0<- est_omega[i]*(1-nu-tau) + nu*est_omega[i]
      value<-value + log(dnorm(est_U[(i-1)*T+1],sd=sqrt(h0)))
      
      for(t in 2:T){
        
        h0<- est_omega[i]*(1-nu-tau) + nu*h0 + tau*est_U[(i-1)*T+t-1]^2
        
        
        value<-value+log(dnorm(est_U[(i-1)*T+t],sd=sqrt(h0)))
      }
      
      
    }
    
    return(-value)
  }
  
  fun_ineq<-function(par){
    tau<-par[1]
    nu<-par[2]
    
    return(c(nu+tau))
  }
  
  estml_GARCH<-try(suppressWarnings(solnp(c(tau,nu), fun = ml_GARCH,
                                          ineqfun = fun_ineq,
                                          ineqLB = c(0.001),
                                          ineqUB = c(0.999),
                                          LB = c(0.001,0.001),
                                          UB = c(0.999,0.999),
                                          control = list(trace = 0))) 
                   ,silent = TRUE)
  if('try-error' %in% class(estml_GARCH)){next}
  
  est_zeta<-estml_GARCH$pars
  
  ###########################################
  ############################################


  # Jackknife
  
  fun_OLS_1 <- function(par){
    
    beta<-par[1]
    phi<-par[2]
    psi<-par[3]
    
    
    A<-fun_A(phi,T/2)
    
    Sigma<-fun_Sigma(psi,T/2)
    
    Sigma_<-solve(Sigma)
    
    l<-rep(1,T/2)
    
    C<-Sigma_%*%l%*%t(l)%*%Sigma_/c(t(l)%*%Sigma_%*%l)
    
    
    value<-0
    for(i in 1:N){
      
      y_i<-matrix(Y[((i-1)*T+1):((i-1)*T+T/2)])
      
      x_i<-matrix(X[((i-1)*T+1):((i-1)*T+T/2)])
      v_i<-A%*%y_i - x_i*beta
      
      value<- value + t(v_i)%*%(Sigma_ - C)%*%v_i
    }
    
    return(value)   
    
  }
  
  estml<-try(suppressWarnings(solnp(c(beta,phi,psi), fun = fun_OLS_1,
                                    # ineqfun = fun_ineq, 
                                    # ineqLB = c(-0.999), 
                                    # ineqUB = c(0.999),
                                    LB = c(-Inf,-0.999,-0.999), 
                                    UB = c(Inf,0.999,0.999),
                                    control = list(trace = 0))) 
             ,silent = TRUE)
  
  if('try-error' %in% class(estml)){next} 
  
  est_gamma_1<-estml$pars
  
  
  fun_OLS_2 <- function(par){
    
    beta<-par[1]
    phi<-par[2]
    psi<-par[3]
    
    
    A<-fun_A(phi,T/2)
    
    Sigma<-fun_Sigma(psi,T/2)
    
    Sigma_<-solve(Sigma)
    
    l<-rep(1,T/2)
    
    C<-Sigma_%*%l%*%t(l)%*%Sigma_/c(t(l)%*%Sigma_%*%l)
    
    
    value<-0
    for(i in 1:N){
      
      y_i<-matrix(Y[((i-1)*T+T/2+1):(i*T)])
      
      x_i<-matrix(X[((i-1)*T+T/2+1):(i*T)])
      v_i<-A%*%y_i - x_i*beta
      
      value<- value + t(v_i)%*%(Sigma_ - C)%*%v_i
    }
    
    return(value)   
    
  }
 
  estml<-try(suppressWarnings(solnp(c(beta,phi,psi), fun = fun_OLS_2,
                                    # ineqfun = fun_ineq, 
                                    # ineqLB = c(-0.999), 
                                    # ineqUB = c(0.999),
                                    LB = c(-Inf,-0.999,-0.999), 
                                    UB = c(Inf,0.999,0.999),
                                    control = list(trace = 0))) 
             ,silent = TRUE)
  
  if('try-error' %in% class(estml)){next} 
  
  est_gamma_2<-estml$pars
     
  ###########################
  est_gamma_unbias_J<- 2*est_gamma - (est_gamma_1+est_gamma_2)/2
  
  est_A<-fun_A(est_gamma_unbias_J[2])
  
  est_Sigma<-fun_Sigma(est_gamma_unbias_J[3])
  
  est_Sigma_<-solve(est_Sigma)
  
  c_Sigma<-c(t(l)%*%est_Sigma_%*%l)
  
  est_V<-(diag(N)%x%est_A)%*%Y - X*est_gamma_unbias_J[1]
  
  est_mu<-(diag(N)%x%(1/c_Sigma*t(l_T)%*%est_Sigma_))%*%est_V
  
  ############################################
  ###########################################
  
  est_B<-fun_B(est_gamma_unbias_J[3])
  
  est_U <- (diag(N)%x%solve(est_B)) %*% (est_V - est_mu%x%l_T)
  
  
  ##########################
  est_omega<-c()
  
  for(i in 1:N){
    est_omega<-c(est_omega,mean((est_U[((i-1)*T+1):(i*T)])^2))
  
  }
  
  
  ##################################
  
  estml_GARCH_0<-try(suppressWarnings(solnp(c(tau,nu), fun = ml_GARCH,
                                            ineqfun = fun_ineq,
                                            ineqLB = c(0.001),
                                            ineqUB = c(0.999),
                                            LB = c(0.001,0.001),
                                            UB = c(0.999,0.999),
                                            control = list(trace = 0))) 
                     ,silent = TRUE)
  if('try-error' %in% class(estml_GARCH_0)){next}
  
  est_zeta_0<-estml_GARCH_0$pars
  
  
  
  
  
  ##########################
  
  est_omega_1<-c()
  for(i in 1:N){
    est_omega_1<-c(est_omega_1,mean((est_U[((i-1)*T+1):((i-1)*T+T/2)])^2))
    
  }
  
  ml_GARCH_1<-function(par){
    
    tau<-par[1]
    nu<-par[2]
    
    value<-0
    for(i in 1:N){
      h0<- est_omega_1[i]*(1-nu-tau) + nu*est_omega_1[i]
      value<-value + log(dnorm(est_U[(i-1)*T+1],sd=sqrt(h0)))
      
      for(t in 2:(T/2)){
        
        h0<- est_omega_1[i]*(1-nu-tau) + nu*h0 + tau*est_U[(i-1)*T+t-1]^2
        
        value<-value+log(dnorm(est_U[(i-1)*T+t],sd=sqrt(h0)))
      }
      
      
    }
    
    return(-value)
  }
  
  
  estml_GARCH<-try(suppressWarnings(solnp(c(tau,nu), fun = ml_GARCH_1,
                                          ineqfun = fun_ineq,
                                          ineqLB = c(0.001),
                                          ineqUB = c(0.999),
                                          LB = c(0.001,0.001),
                                          UB = c(0.999,0.999),
                                          control = list(trace = 0))) 
                   ,silent = TRUE)
  if('try-error' %in% class(estml_GARCH)){next}
  
  est_zeta_1<-estml_GARCH$pars
  
  est_omega_2<-c()
  for(i in 1:N){
    est_omega_2<-c(est_omega_2,mean((est_U[((i-1)*T+T/2+1):((i-1)*T+T)])^2))
    
  }
  
  ml_GARCH_2<-function(par){
    
    tau<-par[1]
    nu<-par[2]
    
    
    value<-0
    for(i in 1:N){
      h0<- est_omega_2[i]*(1-nu-tau) + nu*est_omega_2[i]
      value<-value + log(dnorm(est_U[(i-1)*T+T/2+1],sd=sqrt(h0)))
      
      for(t in (T/2 + 2) : T){
        
        h0<- est_omega_2[i]*(1-nu-tau) + nu*h0 + tau*est_U[(i-1)*T+t-1]^2
        
        
        value<-value+log(dnorm(est_U[(i-1)*T+t],sd=sqrt(h0)))
      }
      
      
    }
    
    return(-value)
  }
  
  
  estml_GARCH<-try(suppressWarnings(solnp(c(tau,nu), fun = ml_GARCH_2,
                                          ineqfun = fun_ineq,
                                          ineqLB = c(0.001),
                                          ineqUB = c(0.999),
                                          LB = c(0.001,0.001),
                                          UB = c(0.999,0.999),
                                          control = list(trace = 0))) 
                   ,silent = TRUE)
  if('try-error' %in% class(estml_GARCH)){next}
  
  est_zeta_2<-estml_GARCH$pars
 
  
  est_zeta_unbias_J <- 2*est_zeta_0 - (est_zeta_1+est_zeta_2)/2


 
  ###########################################
  ############################################
  # ---- deps ----
  suppressPackageStartupMessages({
    library(Rsolnp)
    library(parallel)
  })
    
  # ---- 核心估计函数 ----
  estimate_ARMA_GARCH_FE <- function(y_mat, x_mat, N, T,
                                    fun_A, fun_B,
                                    max_outer_iter = 200,
                                    beta, phi, psi, tau, nu,
                                    tol_outer = 1e-3,
                                    max_inner_iter = 200,
                                    verbose = TRUE) {



 
    
    diff         <- 100
    outer_it     <- 0

    n_cores <- max(1, 25)
    cl <- makeCluster(n_cores)
    on.exit(try(stopCluster(cl), silent = TRUE), add = TRUE)

    # 导出不会频繁变化的对象（动态变化的在循环里再 export）
    clusterExport(cl, varlist = c("y_mat", "x_mat", "N", "T" ), envir = environment())
    clusterEvalQ(cl, {suppressPackageStartupMessages(library(Rsolnp)); NULL})
 
    coef_est <- c(beta, phi, psi, tau, nu)

   
    pre_coef_est <- coef_est + 100
 
    while (diff > tol_outer && outer_it < max_outer_iter) {
      outer_it <- outer_it + 1
      # if (verbose) message(sprintf("[Outer %d] coef = %s", outer_it, paste(round(coef_est, 6), collapse = ", ")))

      beta <- coef_est[1]; phi <- coef_est[2]; psi <- coef_est[3]; tau <- coef_est[4]; nu <- coef_est[5]

      # 4) 在本轮外层参数下，预先计算 A, B（供内层用）
      A_now <- fun_A(phi, T)
      B_now <- fun_B(psi, T)

      # 内层目标函数工厂：给定某个 i 的 y_i, x_i，返回一个仅关于 par = (mu_i, varpi_i) 的函数
      make_inner_obj <- function(y_i, x_i, beta, tau, nu, A, B, T) {
        function(par) {
          mu_i    <- par[1]
          varpi_i <- par[2]

          # u_i = solve(B, A %*% y_i - x_i * beta - mu_i)
          rhs <- A %*% y_i - x_i * beta - mu_i
          u_i <- solve(B, rhs)  # 避免显式求逆

          # GARCH(1,1) 条件方差递推 + 正态对数似然
          h  <- numeric(T)
          h[1] <- varpi_i
          ll  <- dnorm(u_i[1], mean = 0, sd = sqrt(h[1]), log = TRUE)
          if (T > 1) {
            for (t in 2:T) {
              h[t] <- varpi_i + nu * h[t-1] + tau * u_i[t-1]^2
              ll   <- ll + dnorm(u_i[t], mean = 0, sd = sqrt(h[t]), log = TRUE)
            }
          }
          return(-ll) # solnp 做最小化
        }
      }

      # 5) 并行求解 N 个个体的 (mu_i, varpi_i)
      clusterExport(cl,
                    varlist = c("A_now", "B_now", "beta", "tau", "nu", "make_inner_obj"),
                    envir = environment())

      inner_res <- parLapply(cl, 1:N, function(i) {
        y_i <- matrix(y_mat[i, ], ncol = 1)
        x_i <- matrix(x_mat[i, ], ncol = 1)

        target_fun <- make_inner_obj(y_i, x_i, beta, tau, nu, A_now, B_now, T)

        # 初值取 (0, 1)，方差下界设为 1e-3
        ans <- try(suppressWarnings(
          solnp(pars = c(0, 1),
                fun  = target_fun,
                LB   = c(-Inf, 1e-3),
                UB   = c( Inf, Inf),
                control = list(trace = 0, tol = 1e-6, maxiter = 200))
        ), silent = TRUE)

        if (inherits(ans, "try-error")) return(c(NA_real_, NA_real_))
        c(ans$pars[1], ans$pars[2])
      })

      inner_mat <- do.call(rbind, inner_res)
      ok <- stats::complete.cases(inner_mat)
      if (!any(ok)) stop("Inner optimization failed for all i; check data/initialization.")
      mu_est    <- inner_mat[ok, 1]
      varpi_est <- inner_mat[ok, 2]
      # 同步 y/x 的有效样本（极少发生；仅当部分 i 失败时）
      y_eff <- y_mat[ok, , drop = FALSE]
      x_eff <- x_mat[ok, , drop = FALSE]
      N_eff <- nrow(y_eff)

      # 6) 外层目标函数（关于 beta,phi,psi,tau,nu）
      outer_obj <- function(par) {
        beta <- par[1]; phi <- par[2]; psi <- par[3]; tau <- par[4]; nu <- par[5]
        A <- fun_A(phi, T)
        B <- fun_B(psi, T)
        
        total_ll <- 0.0
        for (i in 1:N_eff) {
          mu_i    <- mu_est[i]
          varpi_i <- varpi_est[i]
          y_i <- matrix(y_eff[i, ], ncol = 1)
          x_i <- matrix(x_eff[i, ], ncol = 1)

          rhs <- A %*% y_i - x_i * beta - mu_i
          u_i <- solve(B, rhs)

          h  <- numeric(T)
          h[1] <- varpi_i
          ll  <- dnorm(u_i[1], mean = 0, sd = sqrt(h[1]), log = TRUE)
          if (T > 1) {
            for (t in 2:T) {
              h[t] <- varpi_i + nu * h[t-1] + tau * u_i[t-1]^2
              ll   <- ll + dnorm(u_i[t], mean = 0, sd = sqrt(h[t]), log = TRUE)
            }
          }
          total_ll <- total_ll + ll
        }
        return(-total_ll)
      }

      # 不等式约束：nu + tau ∈ [0.001, 0.999]
      ineq_fun <- function(par) par[4] + par[5]

      # 7) 外层优化
      result_coef <- try(suppressWarnings(
        solnp(pars   = coef_est,
              fun    = outer_obj,
              ineqfun = ineq_fun,
              ineqLB = 0.001, ineqUB = 0.999,
              LB = c(-Inf, -0.999, -0.999, 0.001, 0.001),
              UB = c( Inf,  0.999,  0.999, 0.999, 0.999),
              control = list(trace = 0, tol = 1e-6, maxiter = 1000))
      ), silent = TRUE)

      if (inherits(result_coef, "try-error")) {
        if (verbose) message("Outer solnp failed in this iteration; keep previous coef.")
        break
      }

      coef_est <- as.numeric(result_coef$pars)
      diff <- max(abs(coef_est - pre_coef_est))
      pre_coef_est <- coef_est

      # if (verbose) message(sprintf("  -> diff = %.6g", diff))
    }

    list(
      coef = setNames(coef_est, c("beta","phi","psi","tau","nu")),
      iterations = outer_it,
      converged = (diff <= tol_outer),
      last_diff = diff
    )
  }
  
  
  mat_Y <- matrix(Y, nrow = N, ncol = T, byrow = TRUE)
  mat_X <- matrix(X, nrow = N, ncol = T, byrow = TRUE)
 
  res <- estimate_ARMA_GARCH_FE(
    mat_Y, mat_X, N, T, 
    fun_A = fun_A, fun_B = fun_B,
    beta = beta, phi = phi, psi = psi, tau = tau, nu = nu,
    max_outer_iter = 200, tol_outer = 1e-3, verbose = TRUE
  )

 

  T_half <- T / 2
  
  mat_Y_half1 <- mat_Y[, 1:T_half]
  mat_X_half1 <- mat_X[, 1:T_half]
  
  mat_Y_half2 <- mat_Y[, (T_half+1):T]
  mat_X_half2 <- mat_X[, (T_half+1):T]

  res_1 <- estimate_ARMA_GARCH_FE(
    mat_Y_half1, mat_X_half1, N, T_half, 
    fun_A = fun_A, fun_B = fun_B,
    beta = beta, phi = phi, psi = psi, tau = tau, nu = nu,
    max_outer_iter = 200, tol_outer = 1e-3, verbose = TRUE
  )

  res_2 <- estimate_ARMA_GARCH_FE(
    mat_Y_half2, mat_X_half2, N, T_half, 
    fun_A = fun_A, fun_B = fun_B,
    beta = beta, phi = phi, psi = psi, tau = tau, nu = nu,
    max_outer_iter = 200, tol_outer = 1e-3, verbose = TRUE
  )
  
  est_profile_unbias_J <- c(2 * res$coef - (res_1$coef + res_2$coef) / 2)

  ###########################################
  
  our_method <- rbind(our_method, c(est_gamma, est_zeta))
  our_method_jacknife <- rbind(our_method_jacknife, c(est_gamma_unbias_J, est_zeta_unbias_J))
  
  proflie_method <- rbind(proflie_method, res$coef)
  proflie_method_jacknife <- rbind(proflie_method_jacknife, est_profile_unbias_J)
     
  print(dim(our_method)[1])


 

}
 



our_method[1,]


round(apply(our_method,2,mean)-our_method[1,],3)
round(apply(our_method_jacknife,2,mean)-our_method_jacknife[1,],3)

round(apply(proflie_method,2,mean)-proflie_method[1,],3)
round(apply(proflie_method_jacknife,2,mean)-proflie_method_jacknife[1,],3)

      

round(sqrt(apply(our_method,2,var)),3)
round(sqrt(apply(our_method_jacknife,2,var)),3)

round(sqrt(apply(proflie_method,2,var)),3) 
round(sqrt(apply(proflie_method_jacknife,2,var)),3)









