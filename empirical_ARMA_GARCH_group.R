



fit_group <- final_used$fit$leaf_groups


mean_error_1 <- mean_error_3 <- mean_error_6 <-c()
var_error_1 <- var_error_3 <- var_error_6 <-c()


logLik <- 0
for(i in 1:length(fit_group)){
  
  fit$leaf_groups[[i]] <- sort(fit$leaf_groups[[i]])
  
  # for(left in 1: (T - T_used - 3)){
  # 
  # ymat_used_i <- ymat[fit$leaf_groups[[j]], (left+1):(left + T_used), drop = FALSE]
  # ymat_pred_i <- ymat[fit$leaf_groups[[j]],(left + T_used + 1): T, drop = FALSE]

  
  ymat_used_i <- ymat_used[fit$leaf_groups[[i]],, drop = FALSE]
  
  y_pred_i <- ymat_pred[fit$leaf_groups[[i]],, drop = FALSE] 
 
  fun_OLS_ARMA <- function(par){
    
    phi<-par[1]
    psi<-par[2]
    
    A<-fun_A(phi,T_used)
    
    Sigma<-fun_Sigma(psi,T_used)
    
    Sigma_<-solve(Sigma)
    
    C<-Sigma_%*%l%*%t(l)%*%Sigma_/c(t(l)%*%Sigma_%*%l)
    
    
    
    value<-0
    
    
    for(i in 1: dim(ymat_used_i)[1]){
      
      y_i<-ymat_used_i[i,]
      
      
      v_i<-A%*%y_i 
      
      value<- value + t(v_i)%*%(Sigma_ - C)%*%v_i
    }
    
    return(c(value))   
    
  }
  
  
  # fun_ineq<-function(par){
  #   beta<-par[1]
  #   phi<-par[2]
  #   psi<-par[3]
  #   
  #   return(c(phi+psi))
  # }
  
  estml<-solnp(c(0,0), fun = fun_OLS_ARMA,
               # ineqfun = fun_ineq, 
               # ineqLB = c(-0.999), 
               # ineqUB = c(0.999),
               LB = c(-0.999, -0.999), 
               UB = c(0.999,0.999)) 
  
  #round(2 * (1 - pnorm(abs(estml$pars/sqrt(diag(solve(estml$hessian)))))),3)   
  
  est_gamma <- estml$pars
  if (abs(est_gamma[1]) > 0.9){
    est_gamma <- c(0, 0)
  }
  
  est_A<-fun_A(est_gamma[1],T_used)
  
  est_Sigma<-fun_Sigma(est_gamma[2],T_used)
  
  est_Sigma_<-solve(est_Sigma)
  
  c_Sigma<-c(t(l)%*%est_Sigma_%*%l)
  
  est_V<-(diag(dim(ymat_used_i)[1])%x%est_A)%*%as.vector(t(ymat_used_i)) 
  
  
  est_mu<-(diag(dim(ymat_used_i)[1])%x%(1/c_Sigma*t(l)%*%est_Sigma_))%*%est_V
  
  est_V_pred<- (diag(dim(ymat_used_i)[1])%x%fun_A(est_gamma[1],T-T_used))%*%as.vector(t(y_pred_i))
  
  
  
  est_pred_Sigma<-fun_Sigma(est_gamma[2],T-T_used)
  
  l_T_pred<-matrix(rep(1,T-T_used))
  
  c_Sigma<-c(t(l_T_pred)%*%solve(est_pred_Sigma)%*%l_T_pred)
  
  
  est_pred_mu<- (diag(dim(ymat_used_i)[1])%x% (1/c_Sigma*t(l_T_pred)%*%solve(est_pred_Sigma)))%*%est_V_pred
  
  est_pred_B<-fun_B(est_gamma[2],T-T_used)
  
  est_pred_U <- (diag(dim(ymat_used_i)[1])%x%solve(est_pred_B))%*% (est_V_pred - est_pred_mu%x%l_T_pred)
 
  mean_mse_fun<-function(par, pred_h){
    phi<-par[1]
    psi<-par[2]
    
    beta1<-0
    beta2<-0
    # beta3<-par[4]
    # beta4<-par[5]
    # beta5<-par[6]
    # beta6<-par[7]
    
    beta<-rbind(beta1,beta2)
    
    value<-c()
    
    if(T-T_used-pred_h <=6){
      return(NA)
    }
    for(i in 1:dim(y_pred_i)[1]){
      # mean<- X_pred[(i-1)*(T-T_used)+1,]%*%beta + est_mu[i]
      # 
      # value<-c(value, Y_pred[(i-1)*(T-T_used)+1] - mean) 
      
      for(t in 6:(T-T_used-pred_h)){
        mean<- (phi + psi) * y_pred_i[i, t-1]  + 1 / (1+ psi) * est_mu[i] - (psi + phi) * psi * y_pred_i[i, t-2] + 
          (psi + phi) * psi^2 * y_pred_i[i,t-3] - 
          (psi + phi) * psi^3 * y_pred_i[i,t-4] + 
          (psi + phi) * psi^4 * y_pred_i[i,t-5] 
        
        
        # + psi * last_pred_U[i] 
        t1<-t
        if(pred_h >= 1){
          for(t1 in (t+1):(t+pred_h)){
            mean<- phi*mean + est_mu[i]
          }
        }
        
        value<-c(value, y_pred_i[i, t1] - mean )
      }
      
    }
    
    # return(mean((value - mean(value))^2))
    return(value)
  }
  
  
  # mean_mse_fun<-function(par, pred_h){
  #   phi<-par[1]
  #   psi<-par[2]
  #   
  #   beta1<-0
  #   beta2<-0
  #   # beta3<-par[4]
  #   # beta4<-par[5]
  #   # beta5<-par[6]
  #   # beta6<-par[7]
  #   
  #   beta<-rbind(beta1,beta2)
  #   
  #   value<-c()
  #   
  #   if(T-T_used-pred_h <=6){
  #     return(NA)
  #   }
  #   for(i in 1:dim(y_pred_i)[1]){
  #     # mean<- X_pred[(i-1)*(T-T_used)+1,]%*%beta + est_mu[i]
  #     # 
  #     # value<-c(value, Y_pred[(i-1)*(T-T_used)+1] - mean) 
  #     
  #     for(t in 6:(pred_h)){
  #       mean<- (phi + psi) * y_pred_i[i, t-1]  + 1 / (1+ psi) * est_mu[i] - (psi + phi) * psi * y_pred_i[i, t-2] + 
  #         (psi + phi) * psi^2 * y_pred_i[i,t-3] - 
  #         (psi + phi) * psi^3 * y_pred_i[i,t-4] + 
  #         (psi + phi) * psi^4 * y_pred_i[i,t-5] 
  #       
  #       
  #       # + psi * last_pred_U[i] 
  #       t1<-t
  #       if(pred_h >= 1){
  #         for(t1 in (t+1):(t+pred_h)){
  #           mean<- phi*mean + est_mu[i]
  #         }
  #       }
  #       
  #       value<-c(value, y_pred_i[i, t1] - mean )
  #     }
  #     
  #   }
  #   
  #   # return(mean((value - mean(value))^2))
  #   return(value)
  # }
  # 
  
  
  
  # mean_error_1 <- mean_error_3 <- mean_error_6 <-c()
  mean_error_1 <- c(mean_error_1, mean_mse_fun(est_gamma, 0))
  mean_error_3 <- c(mean_error_3, mean_mse_fun(est_gamma, 2))
  mean_error_6 <- c(mean_error_6, mean_mse_fun(est_gamma, 5)) 
  

  
  est_B<-fun_B(est_gamma[2],T_used)
  
  est_U <- (diag(dim(ymat_used_i)[1])%x%solve(est_B)) %*% (est_V - est_mu%x%l) 
  
  ##########################
  est_omega<-c()
  est_mu4<-c()
  est_mu3<-c()
  for(i in 1:dim(ymat_used_i)[1]){
    est_omega<-c(est_omega,mean((est_U[((i-1)*T_used+1):(i*T_used)])^2))
    est_mu4<-c(est_mu4,mean((est_U[((i-1)*T_used+1):(i*T_used)])^4))
    est_mu3<-c(est_mu3,mean((est_U[((i-1)*T_used+1):(i*T_used)])^3))
  }
  
  
  
  
  ml_GARCH <- function(par){
    
    tau<-par[1]
    nu<-par[2]
    
    
    value<-0 
    for(i in 1:dim(ymat_used_i)[1]){
      h0<- est_omega[i] 
      value<-value + log(dnorm(est_U[(i-1)*T_used+1],sd=sqrt(h0)))
      
      for(t in 2:T_used){
        
        h0<- est_omega[i]*(1-nu-tau) + nu*h0 + tau*est_U[(i-1)*T_used+t-1]^2  
        
        
        value<-value+log(dnorm(est_U[(i-1)*T_used+t],sd=sqrt(h0)))
      }
      
      
    }
    
    return(-value)
  }
  
  
  fun_ineq<-function(par){
    tau<-par[1]
    nu<-par[2]
    
    return(c(nu+tau))
  }
  
  estml_GARCH<-solnp(c(0.1,0.6), fun = ml_GARCH,
                     ineqfun = fun_ineq,
                     ineqLB = c(0.001),
                     ineqUB = c(0.999),
                     LB = c(0.001,0.001),
                     UB = c(0.999,0.999))
  
  est_zeta<-estml_GARCH$pars
  
  
  
  fun_h<-function(par,y,T_used = T_used){
    omega<-par[1]
    tau<-par[2]
    nu<-par[3]
    
    
    h0<-omega 
    h<-c(h0)
    for(t in 2:T_used){
      
      h0<- omega*(1-nu-tau) + nu*h0 + tau*y[t-1]^2
      
      h<-c(h,h0)
    }
    
    return(h)
    
  }
  
  
  
  # pred_h<-0
  # true_U <- as.vector(t(ymat_pred - apply(cbind(ymat_used,ymat_pred), 1, mean)))
  
  # var_error_fun<-function(par, pred_h){
  #   tau<-par[1]
  #   nu<-par[2]
  #   value<-c()
  #   
  #   for(i in 1:dim(ymat_used_i)[1]){
  #     # h0<- est_omega[i]*(1-nu-tau) + nu*est_omega[i]
  #     # for(t in 2:(dim(ymat_used_i)[2]+1)){
  #     #   h0<- est_omega[i]*(1-nu-tau) + nu*h0 + tau*est_U[(i-1)*(dim(ymat_used_i)[2])+t-1]^2
  #     # } 
  #     
  #     h0<- est_omega[i]*(1-nu-tau) + nu*est_omega[i]
  #     # value<-rbind(value, (est_pred_U[(i-1)*(T-T_used)+1])^2 -  h0) 
  #     for(t in 2:(T-T_used-pred_h)){
  #       h0<- est_omega[i]*(1-nu-tau) + nu*h0 + tau*est_pred_U[(i-1)*(T-T_used)+t-1]^2
  #       
  #       h01<- h0
  #       t1<-t
  #       if(pred_h >= 1){
  #         for(t1 in (t+1):(t+pred_h)){
  #           h01<- est_omega[i]*(1-nu-tau) + (tau+nu)*h01
  #         }
  #       }
  #       
  #       value<-rbind(value, c((est_pred_U[(i-1)*(T-T_used)+t1])^2 - h01))
  #       #value<-c(value, abs(est_pred_U[(i-1)*(T-T_used)+t]) - sqrt(h0))
  #       
  #     }}
  #   
  #   return(value)
  #   
  # }
  # 
  # var_error_fun<-function(par, pred_h){
  #   tau<-par[1]
  #   nu<-par[2]
  #   value<-c()
  #   
  #   for(i in 1:N){
  #     h0<- est_omega[i]*(1-nu-tau) + nu*est_omega[i]
  #     # value<-rbind(value, (est_pred_U[(i-1)*(T-T_used)+1])^2 -  h0) 
  #     for(t in 2:(T-T_used-pred_h)){
  #       h0<- est_omega[i]*(1-nu-tau) + nu*h0 + tau*est_pred_U[(i-1)*(T-T_used)+t-1]^2
  #       
  #       h01<- h0
  #       t1<-t
  #       if(pred_h >= 1){
  #         for(t1 in 1:pred_h){
  #           h01<- est_omega[i]*(1-nu-tau) + (tau+nu)*h01
  #         }
  #       }
  #       
  #       value<-rbind(value, c(log(h01) + (est_pred_U[(i-1)*(T-T_used)+t1])^2 / h01))
  #       #value<-c(value, abs(est_pred_U[(i-1)*(T-T_used)+t]) - sqrt(h0))
  #       
  #     }
  #     
  #     
  #   }
  #   
  #   
  #   # return(mean((value - mean(value))^2))
  #   return((value))
  # } 
  
  var_error_fun<-function(par, pred_h){
    tau<-par[1]
    nu<-par[2]
    value<-c()
    
    for(i in 1:N){
      h0<- est_omega[i]*(1-nu-tau) + nu*est_omega[i]
      # value<-rbind(value, (est_pred_U[(i-1)*(T-T_used)+1])^2 -  h0) 
      for(t in 2:(T-T_used-pred_h)){
        h0<- est_omega[i]*(1-nu-tau) + nu*h0 + tau*est_pred_U[(i-1)*(T-T_used)+t-1]^2
        
        h01<- h0
        t1<-t
        if(pred_h >= 1){
          for(t1 in (t+1):(t+pred_h)){
            h01<- est_omega[i]*(1-nu-tau) + (tau+nu)*h01
          }
        }
        
        value<-rbind(value, c(log(h01) + (est_pred_U[(i-1)*(T-T_used)+t1])^2 / h01))
        #value<-c(value, abs(est_pred_U[(i-1)*(T-T_used)+t]) - sqrt(h0))
        
      }
      
      
    }
    
    
    # return(mean((value - mean(value))^2))
    return((value))
  } 
  
  var_error_1 <- c(var_error_1, var_error_fun(est_zeta, pred_h = 0))
  var_error_3 <- c(var_error_3, var_error_fun(est_zeta, pred_h = 2))
  var_error_6 <- c(var_error_6, var_error_fun(est_zeta, pred_h = 5))
  
  
  
  ml_ARMA_GARCH_value_i <- function(par, mu, omega){
    phi<-par[1]
    psi<-par[2]
    tau <- par[3]
    nu <- par[4]
    
    
    A <- fun_A(phi, T_used)
    B <- fun_B(psi, T_used)
    value <- 0
    for(i in 1:dim(ymat_used_i)[1]){
      y_i <- ymat_used_i[i,]
      mu_i <- mu[i]
      omega_i <- omega[i]
      
      u_i <- solve(B) %*% (A %*% y_i   - mu_i)
      
      h0 <- omega_i
      value <- value + log(dnorm(u_i[1], sd = sqrt(h0)))
      
      for (t in 2:T_used) {
        h0 <- omega_i * (1- tau -nu) + nu * h0 + tau * u_i[t - 1]^2
        value <- value + log(dnorm(u_i[t], sd = sqrt(h0)))
      }
    }
    
    return(-value)
  }
  
  
  
  
  
  logLik = logLik + ml_ARMA_GARCH_value_i(c(est_gamma,est_zeta), est_mu, est_omega)
}



round(c(sqrt(mean((mean_error_1)^2)), 
        sqrt(mean((mean_error_3)^2)), 
        sqrt(mean((mean_error_6)^2))) , 3) 

# 
# round(c(sqrt(mean((var_error_1)^2)), 
#         sqrt(mean((var_error_3)^2)), 
#         sqrt(mean((var_error_6)^2))) , 3)


round(c((mean(var_error_1, na.rm = TRUE)), 
        (mean(var_error_3, na.rm = TRUE)), 
        (mean(var_error_6, na.rm = TRUE))) , 3)



2 * logLik + log(T_used * N) * length(fit_group) * 4 


