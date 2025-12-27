library(Rsolnp)
library(VGAM)

beta<-3

phi<-0.3

psi<-0.3

tau<-0.2

nu<-0.5


n<-N<-100

T<-100

sample_size<-500

l_N<-matrix(rep(1,N))

l_T<-l<-matrix(rep(1,T))

J_N<-diag(N) - l_N%*%t(l_N)/N

J_T<-diag(T) - l_T%*%t(l_T)/T

 


lag<-4

 
 
initial_0_gamma<- matrix(c(beta,phi,psi),nrow=1)

initial_Ey_gamma<- matrix(c(beta,phi,psi),nrow=1)

initial_fixedeffect_gamma<- matrix(c(beta,phi,psi),nrow=1)

initial_0_gamma_unbias_J <- matrix(c(beta,phi,psi),nrow=1)

initial_Ey_gamma_unbias_J <- matrix(c(beta,phi,psi),nrow=1)


initial_fixedeffect_gamma_unbias_J <- matrix(c(beta,phi,psi),nrow=1)


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


omega<-matrix(runif(N,1,3))
mu<-rnorm(N)

#########################################
#########################################
while(dim(initial_0_gamma)[1]<sample_size){ 
  
  
  
  
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
                                    UB = c(Inf,0.999,0.999)))
             ,silent = TRUE)
  
  if('try-error' %in% class(estml)){next} 
  
  
  initial_0_gamma_est <-estml$pars
  
  
  #######################################
  #######################################
  
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
                                    UB = c(Inf,0.999,0.999)))
             ,silent = TRUE)
  
  if('try-error' %in% class(estml)){next} 
  
  initial_0_gamma_est_1 <-estml$pars
  
  
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
                                    UB = c(Inf,0.999,0.999)))
             ,silent = TRUE)
  
  if('try-error' %in% class(estml)){next} 
  
  initial_0_gamma_est_2 <-estml$pars

  initial_0_gamma_unbias_J_est <- 2*initial_0_gamma_est - (initial_0_gamma_est_1 + initial_0_gamma_est_2 )/2
  
  
  #########################################################
  #########################################################
  
  
  
  fun_OLS_inital_Ey <- function(par){
    
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
       
      v_i<-A%*%y_i - x_i*beta - c(phi + psi, rep(0, T-1)) * mean(y_i)      
      
      value<- value + t(v_i)%*%(Sigma_ - C)%*%v_i   
    }
    
    return(value)   
    
  }
  
  
  
 
  
  
  
  estml<-try(suppressWarnings(solnp(c(beta,phi,psi), fun = fun_OLS_inital_Ey,
                                    # ineqfun = fun_ineq, 
                                    # ineqLB = c(-0.999), 
                                    # ineqUB = c(0.999),
                                    LB = c(-Inf,-0.999,-0.999), 
                                    UB = c(Inf,0.999,0.999)))
             ,silent = TRUE)
  
  if('try-error' %in% class(estml)){next} 
  
  
  initial_Ey_gamma_est <- estml$pars
  
  
 
  
  # Jackknife
  
  fun_OLS_inital_Ey_1 <- function(par){
    
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
      v_i<-A%*%y_i - x_i*beta  - c(phi + psi, rep(0,  dim(y_i)[1] -1)) * mean(y_i)      
      
      value<- value + t(v_i)%*%(Sigma_ - C)%*%v_i
    }
    
    return(value)   
    
  }
  
  estml<-try(suppressWarnings(solnp(c(beta,phi,psi), fun = fun_OLS_inital_Ey_1,
                                    # ineqfun = fun_ineq, 
                                    # ineqLB = c(-0.999), 
                                    # ineqUB = c(0.999),
                                    LB = c(-Inf,-0.999,-0.999), 
                                    UB = c(Inf,0.999,0.999)))
             ,silent = TRUE)
  
  if('try-error' %in% class(estml)){next} 
  
  initial_Ey_gamma_est_1 <-estml$pars
  
  
  fun_OLS_inital_Ey_2 <- function(par){
    
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
      v_i<-A%*%y_i - x_i*beta  - c(phi + psi, rep(0,  dim(y_i)[1] -1)) * mean(y_i)      
      
      value<- value + t(v_i)%*%(Sigma_ - C)%*%v_i
    }
    
    return(value)   
    
  }
  
  estml<-try(suppressWarnings(solnp(c(beta,phi,psi), fun = fun_OLS_inital_Ey_2,
                                    # ineqfun = fun_ineq, 
                                    # ineqLB = c(-0.999), 
                                    # ineqUB = c(0.999),
                                    LB = c(-Inf,-0.999,-0.999), 
                                    UB = c(Inf,0.999,0.999)))
             ,silent = TRUE)
  
  if('try-error' %in% class(estml)){next} 
  
  initial_Ey_gamma_est_2 <-estml$pars
  
  initial_Ey_gamma_unbias_J_est <- 2 * initial_Ey_gamma_est - (initial_Ey_gamma_est_1 + initial_Ey_gamma_est_2)/2
  
  
  
  ##################################################
   
  
  
  
  #########################################################
  #########################################################
  
  
  
  fun_OLS_inital_fixedeffect <- function(par){
    
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
      
      Sigma_star <- Sigma_ - C
      
      l_star <- c(phi + psi, rep(0,  dim(y_i)[1] -1)) 
      
      C_star <- c(solve(t(l_star)%*%Sigma_star%*%l_star)) * Sigma_star%*%l_star%*%t(l_star)%*%Sigma_star
      
      value<- value + t(v_i)%*%(Sigma_star - C_star)%*%v_i   
    }
    
    return(value)   
    
  }
  
  
  
  
  
  
  
  estml<-try(suppressWarnings(solnp(c(beta,phi,psi), fun = fun_OLS_inital_fixedeffect,
                                    # ineqfun = fun_ineq, 
                                    # ineqLB = c(-0.999), 
                                    # ineqUB = c(0.999),
                                    LB = c(-Inf,-0.999,-0.999), 
                                    UB = c(Inf,0.999,0.999)))
             ,silent = TRUE)
  
  if('try-error' %in% class(estml)){next} 
  
  
  initial_fixedeffect_gamma_est <- estml$pars
   
  
  
  ##############################################
  # Jackknife
  
  fun_OLS_inital_fixedeffect_1 <- function(par){
    
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
      
      Sigma_star <- Sigma_ - C
      
      l_star <- c(phi + psi, rep(0,  dim(y_i)[1] -1)) 
      
      C_star <- c(solve(t(l_star)%*%Sigma_star%*%l_star)) * Sigma_star%*%l_star%*%t(l_star)%*%Sigma_star
      
      value<- value + t(v_i)%*%(Sigma_star - C_star)%*%v_i    
 
    }
    
    return(value)   
    
  }
  
  estml<-try(suppressWarnings(solnp(c(beta,phi,psi), fun = fun_OLS_inital_fixedeffect_1,
                                    # ineqfun = fun_ineq, 
                                    # ineqLB = c(-0.999), 
                                    # ineqUB = c(0.999),
                                    LB = c(-Inf,-0.999,-0.999), 
                                    UB = c(Inf,0.999,0.999)))
             ,silent = TRUE)
  
  if('try-error' %in% class(estml)){next} 
  
  initial_fixedeffect_gamma_est_1 <-estml$pars
  
  
  fun_OLS_inital_fixedeffect_2 <- function(par){
    
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
      v_i<-A%*%y_i - x_i*beta  - c(phi + psi, rep(0,  dim(y_i)[1] -1)) * mean(y_i)      
      
      Sigma_star <- Sigma_ - C
      
      l_star <- c(phi + psi, rep(0,  dim(y_i)[1] -1)) 
      
      C_star <- c(solve(t(l_star)%*%Sigma_star%*%l_star)) * Sigma_star%*%l_star%*%t(l_star)%*%Sigma_star
      
      value<- value + t(v_i)%*%(Sigma_star - C_star)%*%v_i   
    }
    
    return(value)   
    
  }
  
  estml<-try(suppressWarnings(solnp(c(beta,phi,psi), fun = fun_OLS_inital_fixedeffect_2,
                                    # ineqfun = fun_ineq, 
                                    # ineqLB = c(-0.999), 
                                    # ineqUB = c(0.999),
                                    LB = c(-Inf,-0.999,-0.999), 
                                    UB = c(Inf,0.999,0.999)))
             ,silent = TRUE)
  
  if('try-error' %in% class(estml)){next} 
  
  initial_fixedeffect_gamma_est_2 <-estml$pars
  
  initial_fixedeffect_gamma_unbias_J_est <- 2 * initial_fixedeffect_gamma_est - (initial_fixedeffect_gamma_est_1 + initial_fixedeffect_gamma_est_2)/2
  
  
  
  initial_0_gamma <- rbind(initial_0_gamma, initial_0_gamma_est)
  initial_0_gamma_unbias_J <- rbind(initial_0_gamma_unbias_J, initial_0_gamma_unbias_J_est)
  
  
  initial_Ey_gamma <- rbind(initial_Ey_gamma, initial_Ey_gamma_est)
  initial_Ey_gamma_unbias_J <- rbind(initial_Ey_gamma_unbias_J, initial_Ey_gamma_unbias_J_est)  
  
  
  
  initial_fixedeffect_gamma <- rbind(initial_fixedeffect_gamma, initial_fixedeffect_gamma_est)
  initial_fixedeffect_gamma_unbias_J <- rbind(initial_fixedeffect_gamma_unbias_J, initial_fixedeffect_gamma_unbias_J_est)    
  
  
  
  print(dim(initial_0_gamma)[1])
  
  
  
}


initial_0_gamma[1,]
 

round(apply(initial_0_gamma,2,mean)-initial_0_gamma[1,],3)
round(apply(initial_Ey_gamma,2,mean)-initial_Ey_gamma[1,],3)
round(apply(initial_fixedeffect_gamma,2,mean)-initial_fixedeffect_gamma[1,],3)


round(apply(initial_0_gamma_unbias_J,2,mean)-initial_0_gamma[1,],3)
round(apply(initial_Ey_gamma_unbias_J,2,mean)-initial_Ey_gamma[1,],3)
round(apply(initial_fixedeffect_gamma_unbias_J,2,mean)-initial_fixedeffect_gamma[1,],3)



round(sqrt(apply(initial_0_gamma,2,var)),3)
round(sqrt(apply(initial_Ey_gamma,2,var)),3)
round(sqrt(apply(initial_fixedeffect_gamma,2,var)),3)



round(sqrt(apply(initial_0_gamma_unbias_J,2,var)),3)
round(sqrt(apply(initial_Ey_gamma_unbias_J,2,var)),3)
round(sqrt(apply(initial_fixedeffect_gamma_unbias_J,2,var)),3)
