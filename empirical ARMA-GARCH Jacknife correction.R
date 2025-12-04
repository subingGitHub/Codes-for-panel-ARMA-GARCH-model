
fun_A<-function(phi,T_used){
  
  A <- diag(T_used)
  diag(A[-1,-T_used]) <- -phi
  
  return(A)
}
############

fun_A_dot1<-function(phi,T_used){
  
  A_dot1 <- diag(T_used)*0
  diag(A_dot1[-1,-T_used]) <- -1
  return(A_dot1)
}
######
######
fun_B<-function(psi,T_used){
  
  B <- diag(T_used)
  diag(B[-1,-T_used]) <- psi
  
  return(B)
}
########
fun_B_dot1<-function(psi,T_used){
  
  B_dot1 <- diag(T_used)*0
  diag(B_dot1[-1,-T_used]) <- 1
  return(B_dot1)
}
###############
fun_Sigma<-function(psi,T_used){
  B<-fun_B(psi,T_used)
  return(B%*%t(B))
}




fun_OLS_1 <- function(par){
  
  
  phi<-par[1]
  psi<-par[2]
  
  
  beta1<-par[3]
  beta2<-par[4]
  #beta3<-par[4]
  # beta4<-par[5]
  # beta5<-par[6]
  # beta6<-par[7]
  
  
  beta<-rbind(beta1,beta2 )
  
  
  A<-fun_A(phi,floor(T_used/2))
  
  Sigma<-fun_Sigma(psi,floor(T_used/2))
  
  Sigma_<-solve(Sigma)
  
  l<-rep(1,floor(T_used/2))
  
  C<-Sigma_%*%l%*%t(l)%*%Sigma_/c(t(l)%*%Sigma_%*%l)
  
  
  value<-0
  for(i in 1:N){
    
    y_i<-matrix(Y[((i-1)*T_used+1):((i-1)*T_used+floor(T_used/2))])
    
  
    x_i<-X[((i-1)*T_used+1):((i-1)*T_used+floor(T_used/2)),]
    v_i<-A%*%y_i - x_i%*%beta
    
    value<- value + t(v_i)%*%(Sigma_ - C)%*%v_i
  }
  
  return(value)   
  
}

estml<-solnp(est_gamma, fun = fun_OLS_1,
                                  # ineqfun = fun_ineq, 
                                  # ineqLB = c(-0.999), 
                                  # ineqUB = c(0.999),
                                  LB = c(-0.999,-0.999,-Inf,-Inf  ), 
                                  UB = c(0.999,0.999,Inf,Inf  )) 

 

est_gamma_1<-estml$pars


fun_OLS_2 <- function(par){

  phi<-par[1]
  psi<-par[2]
  
  
  beta1<-par[3]
  beta2<-par[4]
  #beta3<-par[4]
  # beta4<-par[5]
  # beta5<-par[6]
  # beta6<-par[7]
  
  
  beta<-rbind(beta1,beta2 )
  
  
  A<-fun_A(phi,floor(T_used/2))
  
  Sigma<-fun_Sigma(psi,floor(T_used/2))
  
  Sigma_<-solve(Sigma)
  
  l<-rep(1,floor(T_used/2))
  
  C<-Sigma_%*%l%*%t(l)%*%Sigma_/c(t(l)%*%Sigma_%*%l)
  
  
  value<-0
  for(i in 1:N){
    
    y_i<-matrix(Y[((i-1)*T_used+floor(T_used/2)+1):(i*T_used)])
    
    x_i<-X[((i-1)*T_used+floor(T_used/2)+1):(i*T_used),]
    v_i<-A%*%y_i - x_i%*%beta
    
    value<- value + t(v_i)%*%(Sigma_ - C)%*%v_i
  }
  
  return(value)   
  
}

estml<- solnp(est_gamma, fun = fun_OLS_2,
                                  # ineqfun = fun_ineq, 
                                  # ineqLB = c(-0.999), 
                                  # ineqUB = c(0.999),
              LB = c(-0.999,-0.999,-Inf,-Inf  ), 
              UB = c(0.999,0.999,Inf,Inf  )) 
            
 

est_gamma_2<-estml$pars

###########################
est_gamma_unbias_J <- 2 * est_gamma - (est_gamma_1+est_gamma_2)/2

# est_gamma_unbias <-c(0.250893,0.151333,0.014291,0.010970)
# Gamma_unbias_<--solve(sec_deri_function(est_gamma_unbias_J))
# Lambda_unbias<-se_lambda(est_gamma_unbias_J) 
#  
#  
# tvalue_gamma_unbias<- est_gamma_unbias_J/sqrt(diag(Gamma_unbias_%*%Lambda_unbias%*% Gamma_unbias_))




est_omega_1<-c()
for(i in 1:N){
  est_omega_1<-c(est_omega_1,mean((est_U[((i-1)*T_used+1):((i-1)*T_used+floor(T_used/2))])^2))
  
}

ml_GARCH_1<-function(par){
  
  tau<-par[1]
  nu<-par[2]
  
  value<-0
  for(i in 1:N){
    h0<- est_omega_1[i] 
    value<-value + log(dnorm(est_U[(i-1)*T_used+1],sd=sqrt(h0)))
    
    for(t in 2:floor(T_used/2)){
      
      h0<- est_omega_1[i]*(1-nu-tau) + nu*h0 + tau*est_U[(i-1)*T_used+t-1]^2
      
      
      value<-value+log(dnorm(est_U[(i-1)*T_used+t],sd=sqrt(h0)))
    }
    
    
  }
  
  return(-value)
}


estml_GARCH<-try(suppressWarnings(solnp(c(0.1,0.8), fun = ml_GARCH_1,
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
  est_omega_2<-c(est_omega_2,mean((est_U[((i-1)*T_used+floor(T_used/2)+1):((i-1)*T_used+T_used)])^2))
  
}

ml_GARCH_2<-function(par){
  
  tau<-par[1]
  nu<-par[2]
  
  value<-0
  for(i in 1:N){
    h0<- est_omega_2[i] 
    value<-value + log(dnorm(est_U[(i-1)*T_used+floor(T_used/2)+1],sd=sqrt(h0)))
    
    for(t in (floor(T_used/2) + 2) : T_used){
      
      h0<- est_omega_2[i]*(1-nu-tau) + nu*h0 + tau*est_U[(i-1)*T_used+t-1]^2
      
      
      value<-value+log(dnorm(est_U[(i-1)*T_used+t],sd=sqrt(h0)))
    }
    
    
  }
  return(-value)
  
}


estml_GARCH<-try(suppressWarnings(solnp(c(0.1,0.8), fun = ml_GARCH_2,
                                        ineqfun = fun_ineq,
                                        ineqLB = c(0.001),
                                        ineqUB = c(0.999),
                                        LB = c(0.001,0.001),
                                        UB = c(0.999,0.999)))
                 ,silent = TRUE)
if('try-error' %in% class(estml_GARCH)){next}

est_zeta_2<-estml_GARCH$pars




est_zeta_unbias_J <- 2*est_zeta - (est_zeta_1+est_zeta_2)/2    

