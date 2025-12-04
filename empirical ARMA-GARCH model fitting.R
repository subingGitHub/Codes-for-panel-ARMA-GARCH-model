library(Rsolnp)


T<-dim(ymat)[2]

n<-N<-dim(ymat)[1]

l_N<-matrix(rep(1,N))

l_T<-l<-matrix(rep(1,T))

J_N<-diag(N) - l_N%*%t(l_N)/N

J_T<-diag(T) - l_T%*%t(l_T)/T




Y<-as.vector(t(ymat))

X<-cbind(as.vector(t(xmat1)), as.vector(t(xmat4)) )
# include zmat5
# remove xmat2  xmat3 zmat4 zmat7


# phi<-psi<-0.1
# 
# beta1<-beta2<-beta3<-beta4<-beta5<-0.2


fun_A<-function(phi){
  
  A <- diag(T)
  diag(A[-1,-T]) <- -phi
  
  return(A)
}
############

fun_A_dot1<-function(phi){
  
  A_dot1 <- diag(T)*0
  diag(A_dot1[-1,-T]) <- -1
  return(A_dot1)
}
######
######
fun_B<-function(psi){
  
  B <- diag(T)
  diag(B[-1,-T]) <- psi
  
  return(B)
}
########
fun_B_dot1<-function(psi){
  
  B_dot1 <- diag(T)*0
  diag(B_dot1[-1,-T]) <- 1
  return(B_dot1)
}
###############
fun_Sigma<-function(psi){
  B<-fun_B(psi)
  return(B%*%t(B))
}


fun_OLS <- function(par){
  
  phi<-par[1]
  psi<-par[2]
  
  
  beta1<-par[3]
  beta2<-par[4]
  
  #beta4<-par[6]
  # beta5<-par[6]
  # beta6<-par[7]
  
  
  beta<-rbind(beta1,beta2)
  
  A<-fun_A(phi)
  
  Sigma<-fun_Sigma(psi)
  
  Sigma_<-solve(Sigma)
  
  C<-Sigma_%*%l%*%t(l)%*%Sigma_/c(t(l)%*%Sigma_%*%l)
  
  value<-0
  for(i in 1:N){
    
    y_i<-ymat[i,]
    
    
    x_i <-X[((i-1)*T+1):(i*T),]
    
    
    v_i<-A%*%y_i - x_i%*%beta
    
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

estml<-solnp(c(0,0, 0,0   ), fun = fun_OLS,
             # ineqfun = fun_ineq, 
             # ineqLB = c(-0.999), 
             # ineqUB = c(0.999),
             LB = c(-0.999,-0.999,-Inf,-Inf ), 
             UB = c(0.999,0.999,Inf,Inf  )) 

#round(2 * (1 - pnorm(abs(estml$pars/sqrt(diag(solve(estml$hessian)))))),3)   

est_gamma<-estml$pars

est_A<-fun_A(est_gamma[1])

est_Sigma<-fun_Sigma(est_gamma[2])

est_Sigma_<-solve(est_Sigma)

c_Sigma<-c(t(l)%*%est_Sigma_%*%l)

est_V<-(diag(N)%x%est_A)%*%Y - X%*%matrix(est_gamma[-c(1,2)])

est_mu<-(diag(N)%x%(1/c_Sigma*t(l_T)%*%est_Sigma_))%*%est_V

dim(est_V)


# mean_mse1<-0
# for(i in 1:N){
#   est_V_<-est_V[((i-1)*T+170):(i*T)]
#   mean_mse1<- mean_mse1 + mean(est_V_^2)/N
# }
# 
# mean_mse1

############################################
###########################################

est_B<-fun_B(est_gamma[2])

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

#c(omega)


ml_GARCH<-function(par){
  
  tau<-par[1]
  nu<-par[2]
  
  
  value<-0 
  for(i in 1:N){
    h0<- est_omega[i]
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

estml_GARCH<-solnp(c(0.1,0.6), fun = ml_GARCH,
                   ineqfun = fun_ineq,
                   ineqLB = c(0.001),
                   ineqUB = c(0.999),
                   LB = c(0.001,0.001), 
                   UB = c(0.999,0.999))

est_zeta<-estml_GARCH$pars

# var_mse<-function(par){
#   
#   tau<-par[1]
#   nu<-par[2]
#   
#   
#   mse<-0
#   for(i in 1:N){
#     
#     value<-c()
#     h0<- est_omega[i]*(1-nu-tau) + nu*est_omega[i]
#     value<-c(value,est_U[(i-1)*T+1]^2 - h0)
#     
#     for(t in 2:T){
#       
#       h0<- est_omega[i]*(1-nu-tau) + nu*h0 + tau*est_U[(i-1)*T+t-1]^2  
#       
#       value<-c(value, est_U[(i-1)*T+t]^2 - h0) 
#     }
#     
#     mse<- mse +  mean(value[170:182]^2)/N
#   }
#   
#   return(mse)
# }
# 
# sqrt(var_mse(est_zeta_unbias))



################################################################
################################################################
# mean part




fir_deri_function<-function(par){
  
  phi<-par[1]
  psi<-par[2]
  
  beta1<-par[3]
  beta2<-par[4]
  
  #beta4<-par[6]
  # beta5<-par[6]
  # beta6<-par[7]
  # 
  
  
  beta<-rbind(beta1,beta2)
  
  A<-fun_A(phi)
  
  B<-fun_B(psi)
  
  A_dot1<-fun_A_dot1(phi)
  B_dot1<-fun_B_dot1(psi)
  
  Sigma<-fun_Sigma(psi)
  
  Sigma_<-solve(Sigma)
  
  Sigma_dot1<-B%*%t(B_dot1)+B_dot1%*%t(B)
  
  C<-Sigma_%*%l%*%t(l)%*%Sigma_/c(t(l)%*%Sigma_%*%l)
  
  
  exp_beta <- rep(0,dim(X)[2])
  
  exp_phi <- - 2* sum(diag( (diag(est_omega)%x%(t(B)%*%C%*%A_dot1%*%solve(A)%*%B))))
  
  exp_psi <- 2* sum(diag( (diag(est_omega)%x%(t(B)%*%Sigma_%*%Sigma_dot1%*%C%*%B)))) - sum(diag(Sigma_dot1%*%C)) * sum(diag( (diag(est_omega)%x% (t(B)%*%C%*%B))))
  
  
  return(matrix(c(exp_phi,exp_psi ,exp_beta)))
  
}


sec_deri_function<-function(par){
  
  
  phi<-par[1]
  psi<-par[2]
  
  beta1<-par[3]
  beta2<-par[4]
  
  #beta4<-par[6]
  # beta5<-par[6]
  # beta6<-par[7]
  
  
  beta<-rbind(beta1,beta2)
  
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
  
  exp_beta_phi <- -2* t( X%*%beta + est_mu%x%l_T)%*% ( diag(N)%x% (t(solve(A))%*%t(A_dot1)%*%(Sigma_ - C)) )%*%X
  
  exp_beta_psi <- -2* ( (t(est_mu)%*%diag(N))%x% (t(l_T)%*%D))%*%X
  
  
  
  exp_phi_phi <- 2* t(X%*%beta + est_mu%x%l_T) %*% (diag(N) %x% M)%*%(X%*%beta + est_mu%x%l_T) + 2*sum(diag((diag(est_omega)%x%(t(B)%*%M%*%B))))
  
  
  exp_phi_psi <- 2* ( (t(est_mu)%*%diag(N))%x% (t(l_T)%*%D%*%A_dot1%*%solve(A)))%*%(X%*%beta + est_mu%x%l_T) + 2*sum(diag( (diag(est_omega)%x%(t(B)%*%D%*%A_dot1%*%solve(A)%*%B))))
  
  
  exp_psi_psi <- - (t(est_mu)%*%est_mu)*(t(l_T)%*%E%*%l_T) - sum(est_omega)*sum(diag(t(B)%*%E%*%B)) 
  
  # exp_phi_phi <- 2* t(X%*%beta + est_mu%x%l_T) %*% (diag(N) %x% M)%*%(X%*%beta + est_mu%x%l_T) + 2*sum(diag((diag(N)%x%(t(B)%*%M%*%B))%*%Sigma_U))
  # 
  # 
  # exp_phi_psi <- 2* ( (t(est_mu)%*%diag(N))%x% (t(l_T)%*%D%*%A_dot1%*%solve(A)))%*%(X%*%beta + est_mu%x%l_T) + 2*sum(diag( (diag(N)%x%(t(B)%*%D%*%A_dot1%*%solve(A)%*%B))%*%Sigma_U))
  # 
  # 
  # exp_psi_psi <- (t(est_mu)%*%diag(N)%*%mu)*(t(l_T)%*%E%*%l_T) + sum(diag( (diag(N)%x%(t(B)%*%E%*%B))%*%Sigma_U))
  # 
  
  result<-rbind(cbind(exp_phi_phi,exp_phi_psi,exp_beta_phi), 
                cbind(exp_phi_psi,exp_psi_psi,exp_beta_psi), 
                cbind(t(exp_beta_phi),t(exp_beta_psi),exp_beta_beta))
  
  
  
  
  
  return(result)
  
}


# the se

fun_cov<-function(par,M,b){
  
  # cmopute the covariance matrix for V'MV +b'V
  sigma2_i<-par[1]
  mu4_i<-par[2]
  mu3_i<-0 #par[3]
  
  value<-sigma2_i*t(b)%*%b + (mu4_i - 3*sigma2_i^2)*sum(diag(M)^2)  + sigma2_i^2*(sum(diag(M%*%t(M)))+sum(diag(M%*%M)) + sum(diag(M))^2)   + 2*mu3_i*sum(diag(M)*b)
  
  
  return(c(value))
}



se_lambda<-function(par){
  
  phi<-par[1]
  psi<-par[2]
  
  beta1<-par[3]
  beta2<-par[4]
  
  #beta4<-par[6]
  # beta5<-par[6]
  # beta6<-par[7]
  
  # phi<-psi<-0.1
  # 
  # beta<-rbind(0.1,0.1,0.1,0.1)
  beta<-rbind(beta1,beta2)
  
  
  
  A<-fun_A(phi)
  B<-fun_B(psi)
  
  A_dot1<-fun_A_dot1(phi)
  B_dot1<-fun_B_dot1(psi)
  
  Sigma_dot1<-B%*%t(B_dot1)+B_dot1%*%t(B)
  
  Sigma<-fun_Sigma(psi)
  
  Sigma_<-solve(Sigma)
  
  C<-Sigma_%*%l%*%t(l)%*%Sigma_/c(t(l)%*%Sigma_%*%l)
  
  
  D <- 2*Sigma_%*%Sigma_dot1%*%C - Sigma_%*%Sigma_dot1%*%Sigma_ - sum(diag(Sigma_dot1%*%C))*C
  
  
  V<-(diag(N)%x%A)%*%Y - X%*%beta
  
  
  se_beta_beta<- 4* t(X)%*%(diag(est_omega)%x% (Sigma_ - C))%*%X
  
  se_phi_phi<-0
  
  se_psi_psi<-0
  
  se_beta_phi<-0
  
  se_beta_psi<-0
  
  
  se_phi_psi<-0
  
  
  for(i in 1:N){
    
    x<-X[((i-1)*T+1):(i*T),]
    M_phi<-  t(B)%*%(Sigma_ - C)%*%A_dot1%*%solve(A)%*%B
    b_phi<- t(B)%*%(Sigma_ - C)%*%A_dot1%*%solve(A)%*%(est_mu[i]*l_T + x%*%beta)
    
    M_psi<- t(B)%*%D%*%B
    b_psi<- t(B)%*%(D+t(D))%*%l_T*est_mu[i] 
    
    
    se_phi_phi<-se_phi_phi + 4 * fun_cov(c(est_omega[i],est_mu4[i]), M_phi, b_phi)
    
    mu_mu_ldl<-0#est_mu[i]^2*c(t(l_T)%*%D%*%l_T)
    se_psi_psi<-se_psi_psi + fun_cov(c(est_omega[i],est_mu4[i]), M_psi, b_psi) #+ # mu_mu_ldl^2 + 2*est_omega[i]*sum(diag(M_psi))*mu_mu_ldl
    
    ########################
    ########################
    
    b_beta_phi_1<-t(B)%*%(Sigma_ - C)%*%x
    b_beta_phi_2<-t(B)%*%(Sigma_ - C)%*%A_dot1%*%solve(A)%*%(est_mu[i]*l_T + x%*%beta)
    
    
    se_beta_phi<-se_beta_phi - 4*apply(diag(c(b_beta_phi_2))%*%b_beta_phi_1,2,sum)*est_omega[i]
    
    
    
    
    
    b_beta_psi_1<- t(B)%*%(Sigma_ - C)%*%x
    b_beta_psi_2<- t(B)%*%(D+t(D))%*%l_T*est_mu[i]
    
    se_beta_psi<-se_beta_psi - 2*apply(diag(c(b_beta_psi_2))%*%b_beta_psi_1,2,sum)*est_omega[i]
    
    
    ########################
    ########################
    
    b_phi_psi_1<- t(B)%*%(Sigma_ - C)%*%A_dot1%*%solve(A)%*%(est_mu[i]*l_T + x%*%beta)
    b_phi_psi_2<-  t(B)%*%(D+t(D))%*%l_T*est_mu[i]
    
    M_phi_psi_1<-t(B)%*%(Sigma_ - C)%*%A_dot1%*%solve(A)%*%B
    
    M_phi_psi_2<-t(B)%*%D%*%B
    UMUUMU<-(est_mu4[i] - 3*est_omega[i]^2)*sum(diag(M_phi_psi_1)*diag(M_phi_psi_2))  + est_omega[i]^2*(sum(diag(M_phi_psi_1%*%t(M_phi_psi_2)))+sum(diag(M_phi_psi_1%*%M_phi_psi_2)) + sum(diag(M_phi_psi_1))*sum(diag(M_phi_psi_2)))
    
    se_phi_psi<- se_phi_psi  + 2*sum(b_phi_psi_1*b_phi_psi_2)*est_omega[i] + 2* UMUUMU #+ 2*sum(diag(M_phi_psi_1))*est_omega[i]*est_mu[i]^2*c(t(l_T)%*%D%*%l_T)
    
  }
  
  
  # result<-rbind(cbind(se_phi_phi, t(se_beta_phi)), 
  #               cbind(matrix(se_beta_phi),se_beta_beta))
  
  result<-rbind(cbind(se_phi_phi,se_phi_psi,t(se_beta_phi)), 
                cbind(se_phi_psi,se_psi_psi,t(se_beta_psi)), 
                cbind(matrix(se_beta_phi),matrix(se_beta_psi),se_beta_beta))
  
  
  return(result)
}





##########################
# SE-mean
est_gamma_unbias_J<-c(0.311516489 ,0.124514136 ,0.006905309 ,0.011582521)

Gamma_<--solve(sec_deri_function(est_gamma))
Lambda<-se_lambda(est_gamma)

tvalue_gamma<- est_gamma/sqrt(diag(Gamma_%*%Lambda%*%Gamma_))




est_gamma_unbias<-est_gamma+solve(sec_deri_function(est_gamma))%*%fir_deri_function(est_gamma)


Gamma_unbias_<--solve(sec_deri_function(est_gamma_unbias))
Lambda_unbias<-se_lambda(est_gamma_unbias) 


tvalue_gamma_unbias<- est_gamma_unbias/sqrt(diag(Gamma_unbias_%*%Lambda_unbias%*% Gamma_unbias_))


round(2 * (1 - pnorm(abs(tvalue_gamma))),3)  
c(round(2 * (1 - pnorm(abs(tvalue_gamma_unbias))),3))

round(est_gamma,4)
c(round(est_gamma_unbias,4))

estml_GARCH$pars



#

##################################
# bias correction

est_gamma_unbias<-est_gamma + solve(sec_deri_function(est_gamma))%*%fir_deri_function(est_gamma)


est_A<-fun_A(est_gamma_unbias[1])

est_Sigma<-fun_Sigma(est_gamma_unbias[2])

est_Sigma_<-solve(est_Sigma)

c_Sigma<-c(t(l)%*%est_Sigma_%*%l)

est_V<-(diag(N)%x%est_A)%*%Y - X%*%est_gamma_unbias[-c(1,2)]

est_mu<-(diag(N)%x%(1/c_Sigma*t(l_T)%*%est_Sigma_))%*%est_V

############################################
###########################################

est_B<-fun_B(est_gamma_unbias[2])

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

estml_GARCH_1<- solnp(est_zeta, fun = ml_GARCH,
                      ineqfun = fun_ineq,
                      ineqLB = c(0.001),
                      ineqUB = c(0.999),
                      LB = c(0.001,0.001), 
                      UB = c(0.999,0.999)) 


est_zeta<-estml_GARCH_1$pars


################################################
##################################################
# GARCH-part

fun_h<-function(par,y){
  omega<-par[1]
  tau<-par[2]
  nu<-par[3]
  
  
  h0<-omega 
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

#mu4<-mean(eps_est^4)
mu4<-mean(c(eps_est^4)[which(eps_est^4<=500)])


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
  
  
  h0<-0
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
    
    x<-X[((r-1)*T+1):(r*T),]
    
    y<-Y[((r-1)*T+1):(r*T)]
    
    ##########################
    x_star<-cbind(x[,1]-mean(x[,1]),x[,2]-mean(x[,2]) )
    
    
    y_star<- y - mean(y)
    
    u_est1_star<-u_est1_ - mean(u_est1_)
    
    ##############################
    
    h_est<-fun_h(c(est_omega[r]-mean(u_est1_)^2,zeta),u_est1_star)
    
    h_tau<-fun_h_tau(c(est_omega[r]-mean(u_est1_)^2,zeta),u_est1_star)
    h_nu<-fun_h_nu(c(est_omega[r]-mean(u_est1_)^2,zeta),h_est)
    
    ###############################################
    
    h0_beta<- -(1-tau-nu)* apply(diag(2*c(u_est1_star))%*%x_star,2,mean)
    
    h0_phi<- -(1-tau-nu)* mean(2*(u_est1_star[-1])*(y_star[-T]))
    h0_psi<- -(1-tau-nu)* mean(2*(u_est1_star[-1])*(u_est1_star[-T]))
    
    h_beta<-c(h0_beta)
    h_phi<-c(h0_phi)
    h_psi<-c(h0_psi)
    
    h0_beta <- - (1-tau-nu)*  apply(diag(2*c(u_est1_star))%*%x_star,2,mean) + nu*h0_beta - tau*2*u_est1_star[1]*x_star[1,]
    
    h0_phi <- -(1-tau-nu)* mean(2*(u_est1_star[-1])*(y_star[-T])) + nu*h0_phi  
    
    h0_psi <- -(1-tau-nu)* mean(2*(u_est1_star[-1])*(u_est1_star[-T])) + nu*h0_psi  
    
    h_beta<-rbind(h_beta,h0_beta)
    h_phi<-c(h_phi,h0_phi)
    h_psi<-c(h_psi,h0_psi)
    
    
    for(t in 3:T){
      
      h0_beta <- - (1-tau-nu)* apply(diag(2*c(u_est1_star))%*%x_star,2,mean) + nu*h0_beta - tau*2*u_est1_star[t-1]*x_star[t-1,]
      
      
      h0_phi <- -(1-tau-nu)* mean(2*(u_est1_star[-1])*(y_star[-T])) + nu*h0_phi - tau*2*u_est1_star[t-1]*y_star[t-2]
      
      h0_psi <- -(1-tau-nu)* mean(2*(u_est1_star[-1])*(u_est1_star[-T])) + nu*h0_psi - tau*2*u_est1_star[t-1]*u_est1_star[t-2]
      
      h_beta<-rbind(h_beta,h0_beta)
      h_phi<-c(h_phi,h0_phi)
      h_psi<-c(h_psi,h0_psi)
      
      
    }
    
    ########################################################
    
    
    P_tau_beta<-P_tau_beta +  apply(diag(c(1/h_est^2*h_tau))%*%h_beta,2,sum)
    
    P_tau_phi<-P_tau_phi + sum(1/h_est^2*h_tau*h_phi)
    
    P_tau_psi<-P_tau_psi + sum(1/h_est^2*h_tau*h_psi)
    
    
    P_nu_beta<-P_nu_beta + apply(diag(c(1/h_est^2*h_nu))%*%h_beta,2,sum)
    
    P_nu_phi<-P_nu_phi + sum(1/h_est^2*h_nu*h_phi)
    
    P_nu_psi<-P_nu_psi + sum(1/h_est^2*h_nu*h_psi)
    
    
  }
  
  
  
  P_est<-1/2 * rbind(c(P_tau_phi,P_tau_psi,P_tau_beta),
                     c(P_nu_phi,P_nu_psi,P_nu_beta))
  
  
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
    
    u_est1_tau<-c(0)
    u_est1_tau_<-0
    for(t1 in 2:T){
      u_est1_tau_<- nu * u_est1_tau_ + u_est1_[t1] - (tau+ nu) * u_est1_[t1 - 1]
      u_est1_tau<-c(u_est1_tau,u_est1_tau_ )
    }
    
    u_est1_nu<-u_est1_tau
    #   c(0,0)
    # u_est1_nu_lag<-u_est1_nu_lagg<-0
    # for(t1 in 3:T){
    #   u_est1_nu_<- 2  * nu * u_est1_nu_lag - nu^2 * u_est1_nu_lagg + tau * u_est1_[t1 - 2]
    #   u_est1_nu<- c(u_est1_nu,u_est1_nu_)
    #   u_est1_nu_lagg<- u_est1_nu_lag
    #   u_est1_nu_lag <- u_est1_nu_
    # }
    
    
    Pi_1<- (u_est1_^2 -h_est) * (1/h_est^2 *h_tau - mean(1/h_est^2 *h_tau))  - 2 *  mean(1/h_est^2 *h_tau * u_est1_tau) * u_est1_
    Pi_2<- (u_est1_^2 -h_est) * (1/h_est^2 *h_nu - mean(1/h_est^2 *h_nu)) - 2 * mean(1/h_est^2 *h_tau* u_est1_nu ) * u_est1_
    
    
    
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


# phi<-psi<-beta1<-beta2<-tau<-nu<-0.2

Lambda_function<-function(par){
  
  phi<-par[1]
  psi<-par[2]
  
  
  beta1<-par[3]
  beta2<-par[4]
  # beta3<-par[4]
  # beta4<-par[5]
  # beta5<-par[6]
  # beta6<-par[7]
  
  
  tau<-par[5]
  nu<-par[6]
  
  beta<-rbind(beta1,beta2)
  
  
  # give the derivative of tilde Q / lambda
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
    
    x<-X[((i-1)*T+1):(i*T),]
    y<-Y[((i-1)*T+1):(i*T)]
    
    
    v<-A%*%y - x%*%beta
    
    u <- solve(B)%*%(v - est_mu[i]) 
    
    
    tildeQ_beta<- -2*diag(c(u))%*%(Sigma_ - C)%*%x
    
    tildeQ_phi<- 2*diag(c(u))%*%lower.tri(t(B)%*%(Sigma_ - C)%*%A_dot1%*%solve(A)%*%B)%*%u + 2*diag(c(u))%*%t(B)%*%(Sigma_ - C)%*%A_dot1%*%solve(A)%*%B%*%(est_mu[i] + x%*%beta)
    
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
    
    
    u_est1_tau<-c(0)
    u_est1_tau_<-0
    for(t1 in 2:T){
      u_est1_tau_<- nu * u_est1_tau_ + u_est1_[t1] - (tau+ nu) * u_est1_[t1 - 1]
      u_est1_tau<-c(u_est1_tau,u_est1_tau_ )
    }
    
    u_est1_nu<-u_est1_tau
    
    # u_est1_nu<-c(0,0)
    # u_est1_nu_lag<-u_est1_nu_lagg<-0
    # for(t1 in 3:T){
    #   u_est1_nu_<- 2  * nu * u_est1_nu_lag - nu^2 * u_est1_nu_lagg + tau * u_est1_[t1 - 2]
    #   u_est1_nu<- c(u_est1_nu,u_est1_nu_)
    #   u_est1_nu_lagg<- u_est1_nu_lag
    #   u_est1_nu_lag <- u_est1_nu_
    # }
    
    
    Pi_1<- (u_est1_^2 -h_est) * (1/h_est^2 *h_tau - mean(1/h_est^2 *h_tau))  - 2 *  mean(1/h_est^2 *h_tau * u_est1_tau) * u_est1_
    Pi_2<- (u_est1_^2 -h_est) * (1/h_est^2 *h_nu - mean(1/h_est^2 *h_nu)) - 2 * mean(1/h_est^2 *h_tau* u_est1_nu ) * u_est1_
    
    
    
    
    M_beta_tau<- M_beta_tau + apply( tildeQ_beta * 1/2*Pi_1,2,sum)
    M_beta_nu<-  M_beta_nu + apply(tildeQ_beta * 1/2*Pi_2,2,sum)
    
    M_phi_tau<-  M_phi_tau +apply( tildeQ_phi * 1/2*Pi_1,2,sum)
    M_phi_nu<-  M_phi_nu +apply( tildeQ_phi * 1/2*Pi_2,2,sum) 
    
    M_psi_tau<-  M_psi_tau + apply( tildeQ_psi * 1/2*Pi_1,2,sum)
    M_psi_nu<-  M_psi_nu + apply( tildeQ_psi * 1/2*Pi_2,2,sum)  
    
  }
  
  M_est<-rbind(c(M_phi_tau,M_phi_nu) ,
               c(M_psi_tau,M_psi_nu) ,
               cbind(M_beta_tau,M_beta_nu))
  
  
  
  
  return(M_est/N/T)
}


#mu4<-mean(eps_est^4)
# mu4<-mean(c(eps_est^4)[which(eps_est^4<=500)])


# est_zeta_unbias<- c(0.2147652 , 0.4004814)


#est_zeta + solve(sec_deri_function_zeta(est_zeta))%*%fir_deri_function_zeta(est_zeta)

#######################################
#######################################
# Gamma_unbias_<--solve(sec_deri_function(est_gamma_unbias))
# Lambda_unbias<-se_lambda(est_gamma_unbias) 
# 
# 
# tvalue_gamma_unbias<- est_gamma_unbias/sqrt(diag(Gamma_unbias_%*%Lambda_unbias%*% Gamma_unbias_))

# Gamma_zeta_unbias_<--solve(sec_deri_function_zeta(est_zeta_unbias))
# 
# Omega_zeta_unbias<-Omega_function(est_zeta_unbias)
# 
# Pi_zeta_unbias <- Pi_function(est_zeta_unbias)
# 
# Lambda_zeta_unbias<- Lambda_function(c(est_gamma_unbias,est_zeta_unbias))
# 
# 
# # Xi_unbias<-(mu4-1)*Omega_zeta_unbias + Pi_zeta_unbias%*%Gamma_unbias_%*%Lambda_unbias%*%Gamma_unbias_%*%t(Pi_zeta_unbias) + Pi_zeta_unbias%*%Gamma_unbias_%*%Lambda_zeta_unbias + t(Pi_zeta_unbias%*%Gamma_unbias_%*%Lambda_zeta_unbias)
# 
# Xi_unbias <- Omega_zeta_unbias + Pi_zeta_unbias%*%Gamma_unbias_%*%Lambda_unbias%*%Gamma_unbias_%*%t(Pi_zeta_unbias) + Pi_zeta_unbias%*%Gamma_unbias_%*%Lambda_zeta_unbias + t(Pi_zeta_unbias%*%Gamma_unbias_%*%Lambda_zeta_unbias)
# 
# 
# tvalue_zeta_unbias<- est_zeta_unbias/sqrt(diag(Gamma_zeta_unbias_%*%Xi_unbias%*%Gamma_zeta_unbias_))




round(c(est_gamma,est_zeta),4)
# round(c(est_gamma_unbias, est_zeta_unbias),6)


# round(2 * (1 - pnorm(abs(tvalue_gamma_unbias))),3)   
# round(2 * (1 - pnorm(abs(tvalue_zeta_unbias))),3)  
# c(round(sqrt(diag(Gamma_unbias_%*%Lambda_unbias%*% Gamma_unbias_)),5) )

 



