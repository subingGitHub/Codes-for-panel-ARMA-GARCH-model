library(Rsolnp)
library(VGAM)

beta<-3

phi<-0.3

psi<-0.3

tau<-0.2

nu<-0.5


n<-N<-100

T<-100



l_N<-matrix(rep(1,N))

l_T<-l<-matrix(rep(1,T))

J_N<-diag(N) - l_N%*%t(l_N)/N

J_T<-diag(T) - l_T%*%t(l_T)/T

sample_size<-500 

a<-b<-c<-matrix(c(beta,phi,psi,tau,nu),nrow=1)

 
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
mu<-rnorm(N)
omega<-matrix(rep(1, N))
  #


mu_all <- matrix(mu,nrow=1)

omega_all <- matrix(omega,nrow=1)


se_mu_all <- c()

se_omega_all <- c()

while(dim(mu_all)[1]<sample_size){ 
  
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
  l<-l_T
  
  fun_OLS <- function(par){
    
    beta<-par[1]
    phi<-par[2]
    psi<-par[3]
    
    
    A<-fun_A(phi)
    
    Sigma<-fun_Sigma(psi)
    
    Sigma_<-solve(Sigma)
    
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
  
  
  # fun_ineq<-function(par){
  #   beta<-par[1]
  #   phi<-par[2]
  #   psi<-par[3]
  #   
  #   return(c(phi+psi))
  # }
  
  estml<-try(suppressWarnings(solnp(c(beta,phi,psi), fun = fun_OLS,
                                    # ineqfun = fun_ineq, 
                                    # ineqLB = c(-0.999), 
                                    # ineqUB = c(0.999),
                                    LB = c(-Inf,-0.999,-0.999), 
                                    UB = c(Inf,0.999,0.999)))
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
 
  est_omega<-c()
  est_mu4<-c()
  est_mu3<-c()

  for(i in 1:N){
    est_omega<-c(est_omega,mean((est_U[((i-1)*T+1):(i*T)])^2))
    est_mu4<-c(est_mu4,mean((est_U[((i-1)*T+1):(i*T)])^4))
    est_mu3<-c(est_mu3,mean((est_U[((i-1)*T+1):(i*T)])^3))

  }


  est_cov1<-c()
  for(i in 1:N){
    eu<-est_U[((i-1)*T+1):(i*T)]
    eu2<-eu[2:T]^2
    eu2_lag<-eu[1:(T-1)]^2

    est_cov1<-c(est_cov1,mean(eu2*eu2_lag))
  }

  est_cov2<-c()
  for(i in 1:N){
    eu<-est_U[((i-1)*T+1):(i*T)]
    eu2<-eu[3:T]^2
    eu2_lag<-eu[1:(T-2)]^2

    est_cov2<-c(est_cov2,mean(eu2*eu2_lag))
  }

  #c(omega)


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
                                          UB = c(0.999,0.999)))
                   ,silent = TRUE)
  if('try-error' %in% class(estml_GARCH)){next}

  est_zeta<-estml_GARCH$pars
  #
  #
  ################################################################
  ################################################################
  
   
  
  fir_deri_function<-function(par){
    
    beta<-par[1]
    phi<-par[2]
    psi<-par[3]
    
    
    A<-fun_A(phi)
    
    B<-fun_B(psi)
    
    A_dot1<-fun_A_dot1(phi)
    B_dot1<-fun_B_dot1(psi)
    
    Sigma<-fun_Sigma(psi)
    
    Sigma_<-solve(Sigma)
    
    Sigma_dot1<-B%*%t(B_dot1)+B_dot1%*%t(B)
    
    C<-Sigma_%*%l%*%t(l)%*%Sigma_/c(t(l)%*%Sigma_%*%l)
    
    
    exp_beta <- 0
    
    exp_phi <- - 2* sum(diag( (diag(est_omega)%x%(t(B)%*%C%*%A_dot1%*%solve(A)%*%B))))
    
    exp_psi <- 2* sum(diag( (diag(est_omega)%x%(t(B)%*%Sigma_%*%Sigma_dot1%*%C%*%B)))) - sum(diag(Sigma_dot1%*%C)) * sum(diag( (diag(est_omega)%x% (t(B)%*%C%*%B))))
    
    
    return(matrix(c(exp_beta,exp_phi,exp_psi)))
    
  }
  #
  
  sec_deri_function<-function(par){
    
    beta<-par[1]
    phi<-par[2]
    psi<-par[3]
    
    
    A<-fun_A(phi)
    
    B<-fun_B(psi)
    
    A_dot1<-fun_A_dot1(phi)
    B_dot1<-fun_B_dot1(psi)
    
    Sigma<-fun_Sigma(psi)
    
    Sigma_<-solve(Sigma)
    
    Sigma_dot1<-B%*%t(B_dot1)+B_dot1%*%t(B)
    
    B11<-B12<-B22<-diag(T)*0
    Sigma_ddot11 <- B_dot1%*%t(B_dot1) + B%*%t(B11)+B11%*%t(B) + B_dot1%*%t(B_dot1)
    
    C<-Sigma_%*%l%*%t(l)%*%Sigma_/c(t(l)%*%Sigma_%*%l)
    
    
    D <- 2*Sigma_%*%Sigma_dot1%*%C - Sigma_%*%Sigma_dot1%*%Sigma_ - sum(diag(Sigma_dot1%*%C))*C
    
    
    M<- t(solve(A))%*%t(A_dot1)%*%(Sigma_ - C)%*%A_dot1%*%solve(A)
    
    
    E<- Sigma_%*%Sigma_ddot11%*%Sigma_ - 2*Sigma_%*%Sigma_dot1%*%Sigma_%*%Sigma_dot1%*%Sigma_ + sum(diag(Sigma_ddot11%*%C))*C +
      2*(sum(diag(Sigma_dot1%*%C))*sum(diag(Sigma_dot1%*%C))*diag(T) - sum(diag(Sigma_dot1%*%Sigma_%*%Sigma_dot1%*%C))*diag(T)  - sum(diag(Sigma_dot1%*%C))*Sigma_%*%Sigma_dot1 + Sigma_%*%Sigma_dot1%*%Sigma_%*%Sigma_dot1 -
           Sigma_%*%Sigma_ddot11 + 2* Sigma_%*%Sigma_dot1%*%Sigma_%*%Sigma_dot1 -  sum(diag(Sigma_dot1%*%C))*Sigma_%*%Sigma_dot1)%*%C
    
    
    
    
    exp_beta_beta <- 2*t(X)%*%( diag(N)%x% (Sigma_ - C) )%*%X
    
    exp_beta_phi <- -2* t( X*beta + est_mu%x%l_T)%*% ( diag(N)%x% (t(solve(A))%*%t(A_dot1)%*%(Sigma_ - C)) )%*%X
    
    exp_beta_psi <- -2* ( (t(est_mu)%*%diag(N))%x% (t(l_T)%*%D))%*%X
    
    
    
    exp_phi_phi <- 2* t(X*beta + est_mu%x%l_T) %*% (diag(N) %x% M)%*%(X*beta + est_mu%x%l_T) + 2*sum(diag((diag(est_omega)%x%(t(B)%*%M%*%B))))
    
    
    exp_phi_psi <- 2* ( (t(est_mu)%*%diag(N))%x% (t(l_T)%*%D%*%A_dot1%*%solve(A)))%*%(X*beta + est_mu%x%l_T) + 2*sum(diag( (diag(est_omega)%x%(t(B)%*%D%*%A_dot1%*%solve(A)%*%B))))
    
    
    exp_psi_psi <- - (t(est_mu)%*%est_mu)*(t(l_T)%*%E%*%l_T) - sum(est_omega)*sum(diag(t(B)%*%E%*%B))
    
    # exp_phi_phi <- 2* t(X*beta + est_mu%x%l_T) %*% (diag(N) %x% M)%*%(X*beta + est_mu%x%l_T) + 2*sum(diag((diag(N)%x%(t(B)%*%M%*%B))%*%Sigma_U))
    #
    #
    # exp_phi_psi <- 2* ( (t(est_mu)%*%diag(N))%x% (t(l_T)%*%D%*%A_dot1%*%solve(A)))%*%(X*beta + est_mu%x%l_T) + 2*sum(diag( (diag(N)%x%(t(B)%*%D%*%A_dot1%*%solve(A)%*%B))%*%Sigma_U))
    #
    #
    # exp_psi_psi <- (t(est_mu)%*%diag(N)%*%mu)*(t(l_T)%*%E%*%l_T) + sum(diag( (diag(N)%x%(t(B)%*%E%*%B))%*%Sigma_U))
    #
    
    
    mat_sec<-rbind(cbind(exp_beta_beta, exp_beta_phi, exp_beta_psi),
                   cbind(0, exp_phi_phi, exp_phi_psi),
                   cbind(0, 0 , exp_psi_psi))
    
    
    colnames(mat_sec) <- NULL
    
    result<-mat_sec +t(mat_sec) -diag(diag(mat_sec))
    
    return(result)
    
  }
  
  
  # the se
  fun_mix_a<-function(e2,M){
    # par1 tau + nu
    
    value<-0
    for(t in 1:T){
      for(t_star in 1:T){
        if(t_star != t){
          value = value + (M[t,t]*M[t_star,t_star] + M[t,t_star]^2 + M[t,t_star]*M[t_star,t])* e2[abs(t-t_star)]
        }
      }
      
    }
    
    return(value)
  }
  
  fun_mix_a1<-function(e2,M,M1){
    # par1 tau + nu
    
    value<-0
    for(t in 1:T){
      for(t_star in 1:T){
        if(t_star != t){
          value = value + (M[t,t]*M1[t_star,t_star] + M[t,t_star]*M1[t,t_star] + M[t,t_star]*M1[t_star,t])* e2[abs(t-t_star)]
        }
      }
      
    }
    
    return(value)
  }
  
  
  fun_cov<-function(par,e2,M,b){
    
    # cmopute the covariance matrix for V'MV +b'V
    sigma2_i<-par[1]
    mu4_i<-par[2]
    mu3_i<-0 #par[3]
    
    value<-sigma2_i*t(b)%*%b + (mu4_i - 3*sigma2_i^2)*sum(diag(M)^2) + fun_mix_a(e2,M) + sigma2_i^2*(sum(diag(M%*%t(M)))+sum(diag(M%*%M)))   + 2*mu3_i*sum(diag(M)*b)
    
    
    return(c(value))
  }
  
  e2_fun<-function(data,lag){
    
    epsu<-data[(lag+1):length(data)]
    
    epsu_lag<-data[1:(length(data)-lag)]
    
    return(cov(epsu^2,epsu_lag^2))
    
  }
  
  # fun_cov<-function(par,par1,M,b){
  #
  #   # cmopute the covariance matrix for V'MV +b'V
  #   sigma2_i<-par[1]
  #   mu4_i<-par[2]
  #   mu3_i<-0 #par[3]
  #
  #   value<-sigma2_i*t(b)%*%b + (mu4_i - 3*sigma2_i^2)*sum(diag(M)^2) +fun_mix_a(par1,M)*(mu4_i - sigma2_i^2) + sigma2_i^2*(sum(diag(M%*%t(M)))+sum(diag(M%*%M)))   + 2*mu3_i*sum(diag(M)*b)
  #
  #
  #   return(c(value))
  # }
  
  
  se_lambda<-function(par){
    
    beta<-par[1]
    phi<-par[2]
    psi<-par[3]
    tau<-par[4]
    nu<-par[5]
    
    A<-fun_A(phi)
    B<-fun_B(psi)
    
    A_dot1<-fun_A_dot1(phi)
    B_dot1<-fun_B_dot1(psi)
    
    Sigma_dot1<-B%*%t(B_dot1)+B_dot1%*%t(B)
    
    Sigma<-fun_Sigma(psi)
    
    Sigma_<-solve(Sigma)
    
    C<-Sigma_%*%l%*%t(l)%*%Sigma_/c(t(l)%*%Sigma_%*%l)
    
    D <- 2*Sigma_%*%Sigma_dot1%*%C - Sigma_%*%Sigma_dot1%*%Sigma_ - sum(diag(Sigma_dot1%*%C))*C
    
    
    
    V<-(diag(N)%x%A)%*%Y - X*beta
    
    
    se_beta_beta<- 4* t(X)%*%(diag(est_omega)%x% (Sigma_ - C))%*%X
    
    se_phi_phi<-0
    
    se_psi_psi<-0
    
    se_beta_phi<-0
    
    se_beta_psi<-0
    
    se_phi_psi<-0
    
    for(i in 1:N){
      est_ui<-est_U[((i-1)*T+1):(i*T)]
      
      e2<-c()
      for(t in 1:T){
        if(t<lag){
          e2<-c(e2, e2_fun(est_ui, t))
        }
        else{
          e2<-c(e2,0)
        }
        
      }
      
      
      x<-X[((i-1)*T+1):(i*T)]
      M_phi<- t(B)%*%(Sigma_ - C)%*%A_dot1%*%solve(A)%*%B
      b_phi<- t(B)%*%(Sigma_ - C)%*%A_dot1%*%solve(A)%*%(est_mu[i]*l_T + x*beta)
      
      M_psi<- t(B)%*%D%*%B
      b_psi<- t(B)%*%(D+t(D))%*%l_T*est_mu[i]
      
      
      se_phi_phi<-se_phi_phi + 4 * fun_cov(c(est_omega[i],est_mu4[i]), e2, M_phi, b_phi)
      
      mu_mu_ldl<-0#est_mu[i]^2*c(t(l_T)%*%D%*%l_T)
      se_psi_psi<-se_psi_psi + fun_cov(c(est_omega[i],est_mu4[i]), e2, M_psi, b_psi) #+ # mu_mu_ldl^2 + 2*est_omega[i]*sum(diag(M_psi))*mu_mu_ldl
      
      ########################
      ########################
      
      b_beta_phi_1<-t(B)%*%(Sigma_  - C)%*%x
      b_beta_phi_2<-t(B)%*%(Sigma_  - C)%*%A_dot1%*%solve(A)%*%(est_mu[i]*l_T + x*beta)
      
      se_beta_phi<-se_beta_phi - 4*sum(b_beta_phi_1*b_beta_phi_2)*est_omega[i]
      
      
      
      b_beta_psi_1<- t(B)%*%(Sigma_  - C)%*%x
      b_beta_psi_2<- t(B)%*%(D+t(D))%*%l_T*est_mu[i]
      
      se_beta_psi<-se_beta_psi - 2*sum(b_beta_psi_1*b_beta_psi_2)*est_omega[i]
      
      
      ########################
      ########################
      
      b_phi_psi_1<- t(B)%*%(Sigma_ - C)%*%A_dot1%*%solve(A)%*%(est_mu[i]*l_T + x*beta)
      b_phi_psi_2<-  t(B)%*%(D+t(D))%*%l_T*est_mu[i]
      
      M_phi_psi_1<-t(B)%*%(Sigma_)%*%A_dot1%*%solve(A)%*%B
      
      M_phi_psi_2<-t(B)%*%D%*%B
      UMUUMU<-(est_mu4[i] - 3*est_omega[i]^2)*sum(diag(M_phi_psi_1)*diag(M_phi_psi_2)) + fun_mix_a1(e2, M_phi_psi_1, M_phi_psi_2)  + est_omega[i]^2*(sum(diag(M_phi_psi_1%*%t(M_phi_psi_2)))+sum(diag(M_phi_psi_1%*%M_phi_psi_2)) )
      
      se_phi_psi<- se_phi_psi  + 2*sum(b_phi_psi_1*b_phi_psi_2)*est_omega[i] + 2* UMUUMU #+ 2*sum(diag(M_phi_psi_1))*est_omega[i]*est_mu[i]^2*c(t(l_T)%*%D%*%l_T)
      
    }
    
    
    mat_sec<-rbind(cbind(se_beta_beta, se_beta_phi, se_beta_psi),
                   cbind(0, se_phi_phi, se_phi_psi),
                   cbind(0, 0 , se_psi_psi))
    
    colnames(mat_sec) <- NULL
    
    result<-mat_sec +t(mat_sec) -diag(diag(mat_sec))
    
    
    return(result)
  }
  
  
  ################################################
  ##################################################
  #GARCH-part
  
  fun_h<-function(par,y){
    omega<-par[1]
    tau<-par[2]
    nu<-par[3]
    
    
    h0<-omega*(1-nu-tau) + nu*omega
    h<-c(h0)
    for(t in 2:T){
      
      h0<- omega*(1-nu-tau) + nu*h0 + tau*y[t-1]^2
      
      h<-c(h,h0)
    }
    
    return(h)
    
  }
  
  
  eps_est<-c()
  
  for(r in 1:n){
    
    u_est1_<-est_U[((r-1)*T+1):(r*T)]
    
    
    eps_est<-c(eps_est,diag(1/sqrt(fun_h(c(est_omega[r],est_zeta),u_est1_)))%*%u_est1_)
    
  }
  
  mu4<-mean(eps_est^4)
  
  
  
  fun_h_nu<-function(par,y){
    omega<-par[1]
    tau<-par[2]
    nu<-par[3]
    
    
    h0<-0
    h<-c(h0)
    for(t in 2:T){
      
      h0<- -omega + nu*h0 + y[t-1]
      
      h<-c(h,h0)
    }
    
    return(h)
    
  }
  
  fun_h_tau<-function(par,y){
    omega<-par[1]
    tau<-par[2]
    nu<-par[3]
    
    
    h0<--omega
    h<-c(h0)
    for(t in 2:T){
      
      h0<- -omega + nu*h0 + y[t-1]^2
      
      h<-c(h,h0)
    }
    
    return(h)
    
  }
  
  
  
  Pi_function<-function(par){
    
    # the second derivate to zeta and lambda
    tau<-par[1]
    nu<-par[2]
    
    zeta<-c(tau,nu)
    
    P_tau_beta<-0
    P_tau_phi<-0
    P_tau_psi<-0
    
    P_nu_beta<-0
    P_nu_phi<-0
    P_nu_psi<-0
    
    for(r in 1:n){
      
      u_est1_<-est_U[((r-1)*T+1):(r*T)]
      
      x<-X[((r-1)*T+1):(r*T)]
      
      y<-Y[((r-1)*T+1):(r*T)]
      
      ##########################
      x_star<-x - mean(x)
      
      y_star<- y - mean(y)
      
      u_est1_star<-u_est1_ - mean(u_est1_)
      
      ##############################
      
      h_est<-fun_h(c(est_omega[r]-mean(u_est1_)^2,zeta),u_est1_star)
      
      h_tau<-fun_h_tau(c(est_omega[r]-mean(u_est1_)^2,zeta),u_est1_star)
      h_nu<-fun_h_nu(c(est_omega[r]-mean(u_est1_)^2,zeta),h_est)
      
      ###############################################
      
      h0_beta<- -(1-tau-nu)* mean(2*(u_est1_star)*x_star)
      h0_phi<- -(1-tau-nu)* mean(2*(u_est1_star[-1])*(y_star[-T]))
      h0_psi<- -(1-tau-nu)* mean(2*(u_est1_star[-1])*(u_est1_star[-T]))
      
      h_beta<-c(h0_beta)
      h_phi<-c(h0_phi)
      h_psi<-c(h0_psi)
      
      h0_beta <- - (1-tau-nu)* mean(2*(u_est1_star)*x_star) + nu*h0_beta - tau*2*u_est1_star[1]*x_star[1]
      
      h0_phi <- -(1-tau-nu)* mean(2*(u_est1_star[-1])*(y_star[-T])) + nu*h0_phi
      
      h0_psi <- -(1-tau-nu)* mean(2*(u_est1_star[-1])*(u_est1_star[-T])) + nu*h0_psi
      
      h_beta<-c(h_beta,h0_beta)
      h_phi<-c(h_phi,h0_phi)
      h_psi<-c(h_psi,h0_psi)
      
      
      for(t in 3:T){
        
        h0_beta <- - (1-tau-nu)* mean(2*(u_est1_star)*x_star) + nu*h0_beta - tau*2*u_est1_star[t-1]*x_star[t-1]
        
        
        h0_phi <- -(1-tau-nu)* mean(2*(u_est1_star[-1])*(y_star[-T])) + nu*h0_phi - tau*2*u_est1_star[t-1]*y_star[t-2]
        
        h0_psi <- -(1-tau-nu)* mean(2*(u_est1_star[-1])*(u_est1_star[-T])) + nu*h0_psi - tau*2*u_est1_star[t-1]*u_est1_star[t-2]
        
        h_beta<-c(h_beta,h0_beta)
        h_phi<-c(h_phi,h0_phi)
        h_psi<-c(h_psi,h0_psi)
        
        
      }
      
      ########################################################
      
      
      P_tau_beta<-P_tau_beta + sum(1/h_est^2*h_tau*h_beta)
      
      P_tau_phi<-P_tau_phi + sum(1/h_est^2*h_tau*h_phi)
      
      P_tau_psi<-P_tau_psi + sum(1/h_est^2*h_tau*h_psi)
      
      
      P_nu_beta<-P_nu_beta + sum(1/h_est^2*h_nu*h_beta)
      
      P_nu_phi<-P_nu_phi + sum(1/h_est^2*h_nu*h_phi)
      
      P_nu_psi<-P_nu_psi + sum(1/h_est^2*h_nu*h_psi)
      
      
    }
    
    
    P_est<-1/2 * rbind(c(P_tau_beta,P_tau_phi,P_tau_psi),
                       c(P_nu_beta,P_nu_phi,P_nu_psi))
    
    
    return(P_est)
    
  }
  
  
  
  Omega_function<-function(par){
    tau<-par[1]
    nu<-par[2]
    
    zeta<-c(tau,nu)
    O_11<-0
    O_12<-0
    O_22<-0
    for(r in 1:n){
      
      u_est1_<-est_U[((r-1)*T+1):(r*T)]
      
      h_est<-fun_h(c(est_omega[r],zeta),u_est1_)
      
      h_tau<-fun_h_tau(c(est_omega[r],zeta),u_est1_)
      h_nu<-fun_h_nu(c(est_omega[r],zeta),h_est)
      
      Pi_1<-1/h_est*h_tau - mean(1/h_est^2*h_tau)*h_est
      Pi_2<-1/h_est*h_nu - mean(1/h_est^2*h_nu)*h_est
      
      
      
      O_11<-O_11 + sum(Pi_1*Pi_1)
      O_12<-O_12 + sum(Pi_1*Pi_2)
      
      O_22<-O_22 + sum(Pi_2*Pi_2)
      
      # J_1_11<-J_1_11+ sum(  ( 2*u_est1_[-1]^2/(h_est[-1])^3 - 1/(h_est[-1])^2 ) *   (u_est1_[-T]^2  )^2)
      # J_1_12<-J_1_12+ sum(( 2*u_est1_[-1]^2/(h_est[-1])^3 - 1/(h_est[-1])^2 ) * (u_est1_[-T]^2 )*(h_est[-T] ) )
      #
      
    }
    
    O_est<-1/4*rbind(c(O_11,O_12),
                     c(O_12,O_22))
    
    # solve(J_1_est)
    return(O_est)
    
    
  }
  
  
  
  sec_deri_function_zeta<-function(par){
    
    tau<-par[1]
    nu<-par[2]
    zeta<-c(tau,nu)
    
    
    J_1_11<-0
    J_1_12<-0
    J_1_22<-0
    for(r in 1:n){
      
      u_est1_<-est_U[((r-1)*T+1):(r*T)]
      
      h_est<-fun_h(c(est_omega[r],zeta),u_est1_)
      
      h_tau<-fun_h_tau(c(est_omega[r],zeta),u_est1_)
      h_nu<-fun_h_nu(c(est_omega[r],zeta),h_est)
      #h_nu_omega<-fun_h_nu_omega(c(est_omega[r],zeta),h_omega,h_nu)
      
      
      J_1_11<-J_1_11+ sum(h_tau^2/(h_est)^2)
      J_1_12<-J_1_12+ sum(h_tau*h_nu/(h_est)^2)
      
      J_1_22<-J_1_22 + sum(h_nu^2/(h_est)^2)
      
    }
    
    J_1_est<--1/2*rbind(c(J_1_11,J_1_12),
                        c(J_1_12,J_1_22))
    
    # solve(J_1_est)
    return(J_1_est)
  }
  
  fir_deri_function_zeta<-function(par){
    
    tau<-par[1]
    nu<-par[2]
    zeta<-c(tau,nu)
    J_1<-0
    J_2<-0
    
    for(r in 1:n){
      
      u_est1_<-est_U[((r-1)*T+1):(r*T)]
      
      h_est<-fun_h(c(est_omega[r],zeta),u_est1_)
      
      h_tau_est<-fun_h_tau(c(est_omega[r],zeta),u_est1_)
      h_nu_est<-fun_h_nu(c(est_omega[r],zeta),h_est)
      
      exp_uh<- mu4-1
      
      #mean((u_est1_^2 - h_est)^2/h_est^2
      
      
      omega_est_omega<--1/(1-tau-nu)*((1-nu)*mean(u_est1_^2-h_est)  - 1/T*(tau+nu)*u_est1_[T]^2)
      
      sum_1<-sum( (h_tau_est) / (h_est)^2 )
      
      L_2<- 1/(1-tau-nu)*(tau+nu)/(1-nu)*mean((u_est1_^2-h_est)/h_est^2)*u_est1_[T]^2
      
      #
      
      c_1<- mean(h_tau_est / h_est^2)
      
      M_25<- -tau*mean(u_est1_^2)*( mean(u_est1_[-T]^2/(h_est[-1])^2)
                                    + nu^2*mean(u_est1_[-c(T,T-1)]^2/(h_est[-c(1,2)])^2)
                                    + nu^4*mean(u_est1_[-c(T,T-1,T-2)]^2/(h_est[-c(1,2,3)])^2)
                                    + nu^6*mean(u_est1_[-c(T,T-1,T-2,T-3)]^2/(h_est[-c(1,2,3,4)])^2))
      
      M_34<- mean(h_tau_est / (h_est)^2 )*mean(u_est1_^2)
      #mean( 1/tau* (h_est - est_omega[r]) / (h_est)^2)*mean(u_est1_^2)
      
      
      J_1<-J_1 -  1/(1-tau-nu) * exp_uh + L_2 + (1-tau-nu)/(1-nu)*omega_est_omega*sum_1 - c_1*sum(u_est1_^2-h_est) + M_25 +  M_34
      
      
      L_2_2 <- 1/(1-tau-nu)*tau/(1-nu)^2*(tau+nu)*mean((u_est1_^2-h_est)/h_est^2)*u_est1_[T]^2
      
      sum_2<-sum(h_nu_est / (h_est)^2 )
      
      
      c_2<- mean(h_nu_est / (h_est)^2 )
      
      M_25_2 <- -tau^2*mean(u_est1_^2)*(nu*mean(u_est1_[-c(T,T-1)]^2/(h_est[-c(1,2)])^2)
                                        + 2*nu^2*mean(u_est1_[-c(T,T-1,T-2)]^2/(h_est[-c(1,2,3)])^2)
                                        + 3*nu^3*mean(u_est1_[-c(T,T-1,T-2,T-3)]^2/(h_est[-c(1,2,3,4)])^2)
                                        + 4*nu^4*mean(u_est1_[-c(T,T-1,T-2,T-3,T-4)]^2/(h_est[-c(1,2,3,4,5)])^2)
                                        + 5*nu^5*mean(u_est1_[-c(T,T-1,T-2,T-3,T-4,T-5)]^2/(h_est[-c(1,2,3,4,5,6)])^2))
      
      
      M_34_2 <- (1-tau-nu)/(1-nu)*mean(h_nu_est / (h_est)^2)*mean(u_est1_^2)
      J_2<-J_2 - 1/(1-tau-nu)*tau/(1-nu) * exp_uh + L_2_2 + (1-tau-nu)/(1-nu)*omega_est_omega*sum_2 - c_2*sum(u_est1_^2-h_est) + M_25_2 +  M_34_2  #(1-tau-nu)/(1-nu)* ((1-tau-nu)/(1-nu) -1 -nu/(1-nu)+nu*(1-tau-nu)/(1-nu)^2)  * var(u_est1_^2)* mean(1/h_est^2)
      
    }
    
    J_est<- 1/2*rbind(J_1,J_2)
    
    return(J_est)
    
  }
  
  Lambda_function<-function(par){
    beta<-par[1]
    phi<-par[2]
    psi<-par[3]
    tau<-par[4]
    nu<-par[5]
    
    zeta<-c(tau,nu)
    
    # give the derivative of tilde Q / lambda
    A<-fun_A(phi)
    
    B<-fun_B(psi)
    
    A_dot1<-fun_A_dot1(phi)
    B_dot1<-fun_B_dot1(psi)
    
    Sigma<-fun_Sigma(psi)
    
    Sigma_<-solve(Sigma)
    
    Sigma_dot1<-B%*%t(B_dot1)+B_dot1%*%t(B)
    
    C<-Sigma_%*%l%*%t(l)%*%Sigma_/c(t(l)%*%Sigma_%*%l)
    
    
    D <- 2*Sigma_%*%Sigma_dot1%*%C - Sigma_%*%Sigma_dot1%*%Sigma_ - sum(diag(Sigma_dot1%*%C))*C
    
    M_beta_tau<-0
    M_beta_nu<-0
    
    M_phi_tau<-0
    M_phi_nu<-0
    
    M_psi_tau<-0
    M_psi_nu<-0
    ######################
    for(i in 1:n){
      
      x<-X[((i-1)*T+1):(i*T)]
      y<-Y[((i-1)*T+1):(i*T)]
      
      
      v<-A%*%y - x*beta
      
      u <- solve(B)%*%(v - est_mu[i]) 
      
      
      tildeQ_beta<- -2*diag(c(u))%*%(Sigma_ - C)%*%x
      
      tildeQ_phi<- 2*diag(c(u))%*%lower.tri(t(B)%*%(Sigma_ - C)%*%A_dot1%*%solve(A)%*%B)%*%u + 2*diag(c(u))%*%t(B)%*%(Sigma_ - C)%*%A_dot1%*%solve(A)%*%B%*%(est_mu[i] + x*beta)
      
      tildeQ_psi<- diag(c(u))%*%lower.tri(t(B)%*%D%*%B)%*%u + 2 * diag(c(u))%*%t(B)%*%D%*%l* est_mu[i]
      
      # diag(c(v))%*%D%*%v - 2 * est_omega[i]*(sum(diag((t(B)%*%Sigma_%*%Sigma_dot1%*%C%*%B))) - sum(diag(Sigma_dot1%*%C)) * sum(diag((t(B)%*%C%*%B))))
      
      u_est1_<-est_U[((i-1)*T+1):(i*T)]
      ##########################
      
      
      # x_star<-x - t(l)%*%Sigma_%*%x/c(t(l)%*%Sigma_%*%l)
      #
      # y_star<- y -  t(l)%*%Sigma_%*%y/c(t(l)%*%Sigma_%*%l)
      #
      # u_est1_star<-u_est1_ -  t(l)%*%Sigma_%*%u_est1_/c(t(l)%*%Sigma_%*%l)
      
      ##############################
      
      h_est<-fun_h(c(est_omega[i],zeta),u_est1_)
      
      h_tau<-fun_h_tau(c(est_omega[i],zeta),u_est1_)
      h_nu<-fun_h_nu(c(est_omega[i],zeta),h_est)
      
      
      
      
      M_beta_tau<- M_beta_tau + sum( tildeQ_beta * 1/2*(u_est1_^2 -h_est) * (1/h_est^2 *h_tau - mean(1/h_est^2 *h_tau)))
      M_beta_nu<-  M_beta_nu + sum( tildeQ_beta * 1/2*(u_est1_^2 -h_est)* (1/h_est^2 *h_nu - mean(1/h_est^2 *h_nu)))  
      
      M_phi_tau<-  M_phi_tau + sum( tildeQ_phi * 1/2*(u_est1_^2 -h_est)* (1/h_est^2 *h_tau - mean(1/h_est^2 *h_tau)))
      M_phi_nu<-  M_phi_nu + sum( tildeQ_phi * 1/2*(u_est1_^2 -h_est)* (1/h_est^2 *h_nu - mean(1/h_est^2 *h_nu)))  
      
      M_psi_tau<-  M_psi_tau + sum( tildeQ_psi * 1/2*(u_est1_^2 -h_est)* (1/h_est^2 *h_tau - mean(1/h_est^2 *h_tau)))
      M_psi_nu<-  M_psi_nu + sum( tildeQ_psi * 1/2*(u_est1_^2 -h_est)* (1/h_est^2 *h_nu - mean(1/h_est^2 *h_nu)))  
      
    }
    
    M_est<-rbind(c(M_beta_tau,M_beta_nu),
                 c(M_phi_tau,M_phi_nu),
                 c(M_psi_tau,M_psi_nu))
    
    
    
    
    return(M_est/N/T)
  }
  
   
  ##################################
  ##################################
  # analytical bias correction
  est_gamma_unbias<-est_gamma+solve(sec_deri_function(est_gamma))%*%fir_deri_function(est_gamma)
  
  est_A<-fun_A(est_gamma_unbias[2])
  
  est_Sigma<-fun_Sigma(est_gamma_unbias[3])
  
  est_Sigma_<-solve(est_Sigma)
  
  c_Sigma<-c(t(l)%*%est_Sigma_%*%l)
  
  est_V<-(diag(N)%x%est_A)%*%Y - X*est_gamma_unbias[1]
  
  est_mu<-(diag(N)%x%(1/c_Sigma*t(l_T)%*%est_Sigma_))%*%est_V
  
  
  
  
  
  
  
  ############################################
  ###########################################
  
  est_B<-fun_B(est_gamma_unbias[3])
  
  est_U <- (diag(N)%x%solve(est_B)) %*% (est_V - est_mu%x%l_T)
  
  
  ##########################
  est_omega<-c()
  est_mu4<-c()
  est_mu3<-c()
  for(i in 1:N){
    est_omega<-c(est_omega,mean((est_U[((i-1)*T+1):(i*T)])^2))
    est_mu4<-c(est_mu4,mean((est_U[((i-1)*T+1):(i*T)])^4))
    est_mu3<-c(est_mu3,mean((est_U[((i-1)*T+1):(i*T)])^3))
    
  }
  

  
  ##################################
  
  estml_GARCH_0<-try(suppressWarnings(solnp(c(tau,nu), fun = ml_GARCH,
                                            ineqfun = fun_ineq,
                                            ineqLB = c(0.001),
                                            ineqUB = c(0.999),
                                            LB = c(0.001,0.001),
                                            UB = c(0.999,0.999)))
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
                                          UB = c(0.999,0.999)))
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
                                          UB = c(0.999,0.999)))
                   ,silent = TRUE)
  if('try-error' %in% class(estml_GARCH)){next}
  
  est_zeta_2<-estml_GARCH$pars
  
  
  
  
  
  eps_est<-c()
  
  moment2_h_est <- c()
  moment2_uh_est <- c()
  
  for(r in 1:n){
    
    u_est1_<-est_U[((r-1)*T+1):(r*T)]
   
      
    eps_est<-c(eps_est,diag(1/sqrt(fun_h(c(est_omega[r],est_zeta_1),u_est1_)))%*%u_est1_) 
    h_est <- fun_h(c(est_omega[r],est_zeta_1),u_est1_)
    
    moment2_h_est <- rbind(moment2_h_est, mean(h_est^2))
    moment2_uh_est <- rbind(moment2_uh_est, mean((u_est1_^2 - h_est)^2))
    
  }
  
  mu4<-mean(eps_est^4)
  
  est_zeta_unbias_J <- 2*est_zeta_0 - (est_zeta_1+est_zeta_2)/2
  
  if(sum(est_zeta_unbias_J)>1){
 
    next
  }
 
  
  tau_est <- est_zeta_unbias_J[1]
  nu_est <- est_zeta_unbias_J[2]
  
  se_omega <- (1 - nu_est)/ (1 - tau_est - nu_est) * sqrt( moment2_uh_est)

  
  
  mu_all <- rbind(mu_all, c(est_mu))
  
  se_mu <- sqrt(T / c_Sigma * est_omega)
  
  se_mu_all <- rbind(se_mu_all, se_mu)
  
  
  omega_all <- rbind(omega_all, c(est_omega))
  
  
  se_omega_all <- rbind(se_omega_all, se_omega)
  
 
  # 
  
  
  print(dim(mu_all)[1])
  
 
}


a[1,]

#apply(a,2,mean)

 

mu_all_bias = apply(mu_all,2,mean)-mu_all[1,]

hist(mu_all_bias, freq=FALSE, main="mu_bias", xlab="value", border="black")
# 计算密度估计
dens <- density(mu_all_bias)
# 在同一图上添加密度线
lines(dens, col="red", lwd=2)


omega_all_bias = apply(omega_all,2,mean)-omega_all[1,]

hist(omega_all_bias, freq=FALSE, main="omega", xlab="value", border="black")
# 计算密度估计
dens <- density(omega_all_bias)
# 在同一图上添加密度线
lines(dens, col="red", lwd=2)





mu_SD_AD <- round(sqrt(apply(mu_all,2,var)),3)  / (round(apply(se_mu_all,2,mean),3) / sqrt(T))


hist(mu_SD_AD, freq=FALSE, main="mu_SD_AD", xlab="value", border="black")
# 计算密度估计
dens <- density(mu_SD_AD)
# 在同一图上添加密度线
lines(dens, col="red", lwd=2)


# 
omega_SD_AD <- round(sqrt(apply(omega_all,2,var)),3)  / (round(apply(se_omega_all,2,mean),3) / sqrt(T))


hist(omega_SD_AD, freq=FALSE, main="omega_SD_AD", xlab="value", border="black")
# 计算密度估计
dens <- density(omega_SD_AD)
# 在同一图上添加密度线
lines(dens, col="red", lwd=2)

 

 
