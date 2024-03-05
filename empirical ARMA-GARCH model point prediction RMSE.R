k <- 12

T_used<- 96

# only pred_h = 0
pred_h<-0

library(forecast)
library(readr)


setwd("C:/Users/subing/OneDrive - The University of Hong Kong - Connect/桌面/monthly data all")


#Short_term_interest_rates <- read_csv("Short-term interest rates.csv")
# 
#Long_term_interest_rates <- read_csv("Long-term interest rates.csv")
#CPI_total <- read_csv("CPI_total.csv")
#total_production <- read_csv("total production.csv")
# Construction <- read_csv("Construction.csv")
# PPI_manuf <- read_csv("PPI_manuf.csv")
# PPI_total <- read_csv("PPI_total.csv")
# CPI_food <- read_csv("CPI_food.csv")
# CPI_other <- read_csv("CPI_other.csv")

PPI_domestic <- read_csv("PPI_domestic.csv")

CPI_energy <- read_csv("CPI_energy.csv")
Manufacturing <- read_csv("Manufacturing.csv")


PPI_domestic <- na.omit(PPI_domestic)


colname<-Reduce(intersect, list(CPI_energy$Location, PPI_domestic$Location, #PPI_total$Location,
                                Manufacturing$Location))


# colname<-setdiff(colname, c("Estonia"))

library('dse')

fun_data<-function(data,difflog = TRUE){
  data<-as.data.frame(data)
  
  data<-data[data[,1] %in% colname,]
  
  # remove location
  data<- data[,-1]
  
  data<- as.matrix(data)
  data<-apply(data,2, as.numeric)
  
  for(i in 1:dim(data)[1]){
    
    data_i <-as.numeric(data[i,])
    if(difflog){
      data_i_diff<-diff(log(data_i))
    }else{
      data_i_diff<-diff((data_i))
    }
    
    if(is.na(mean(data_i_diff,na.rm=T))){
      data_i_diff[is.na(data_i_diff)] <- mean(data,na.rm=T)
    }else{
      data_i_diff[is.na(data_i_diff)] <- mean(data_i_diff,na.rm=T)
    }
    
    # remove 0
    data[i,-1] <- data_i_diff
    
  }
  
  data<- data[,-1]
  
  # 10 years data
  return(data[,((12 * (k -1) + 1) + 5 ):( (12 * (k -1 + 10)) + 6)])
  #return(data[,((12 * (k -1) + 1) - 5 ):( (12 * (k -1 + 10))  )])
  
}

# total has 24 years data 

# Share_prices_<-fun_data(Share_prices)

# CPI_total<-fun_data(CPI_total)

# Short_term_interest_rates<-fun_data(Short_term_interest_rates, FALSE)
# 
# Long_term_interest_rates<-fun_data(Long_term_interest_rates, FALSE)
# total_production<- fun_data(total_production)

# PPI_manuf <- fun_data(PPI_manuf)

# CPI_total <-  fun_data(CPI_total)

# Construction <- fun_data(Construction)
# CPI_food <-fun_data(CPI_food)
#CPI_other<-fun_data(CPI_other)

# PPI_total <-  fun_data(PPI_total)

PPI_domestic <- fun_data(PPI_domestic)
CPI_energy <-fun_data(CPI_energy)
Manufacturing<-fun_data(Manufacturing)



ymat<- PPI_domestic[,-1]

T_len <- dim(ymat)[2] + 1

xmat1<- CPI_energy[,- T_len] - apply(CPI_energy[,- T_len], 1, mean)

xmat4<- Manufacturing[,- T_len]  - apply(Manufacturing[,- T_len], 1, mean)



# xmat1<- CPI_food[,- T_len] - apply(CPI_food[,- T_len], 1, mean)
# xmat2<-  CPI_other[,- T_len]  - apply(CPI_other[,- T_len], 1, mean)
#xmat3<-  Short_term_interest_rates[,- T_len]  -  apply(Short_term_interest_rates[,- T_len], 1, mean)
#Long_term_interest_rates[,- T_len]  -  apply(Long_term_interest_rates[,- T_len], 1, mean)


# xmat2 <- xmat4 <- xmat3


library(Rsolnp)


T<-dim(ymat)[2]

n<-N<-dim(ymat)[1]

Y<-as.vector(t(ymat))

X<-cbind(as.vector(t(xmat1)), as.vector(t(xmat4)))




T_pred<-T - T_used 

dim(ymat)

 
l_N<-matrix(rep(1,N))

l_T_used<-l<-matrix(rep(1,T_used))

J_N<-diag(N) - l_N%*%t(l_N)/N

J_T_used<-diag(T_used) - l_T_used%*%t(l_T_used)/T_used
# include zmat5
# remove xmat2  xmat3 zmat4 zmat7


# phi<-psi<-par[2].1
# 
# beta1<-beta2<-beta3<-beta4<-beta5<-0.2

#est_value<-c(est_gamma ,est_zeta )

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

est_U_pred_2<-c()
error_pred_2<-var_error_2<-c()
 

mean_mse_one<-c()

 
for(i in 1:N){
  
  
  y_i<-as.matrix(ymat[i,1:T_used])
  
  y_pred_i<-as.matrix(ymat[i,(T_used+1):T])
  
  x1_i<-xmat1[i,1:T_used]
  # x2_i<-xmat2[i,1:T_used]
  # x3_i<-xmat3[i,1:T_used]
  x4_i<-xmat4[i,1:T_used]
  # 
  
  x1_pred_i<-xmat1[i,(T_used+1):T]
  # x2_pred_i<-xmat2[i,(T_used+1):T]
  # x3_pred_i<-xmat3[i,(T_used+1):T]
  x4_pred_i<-xmat4[i,(T_used+1):T]
  # 
  
  
  
  x_i <- cbind(x1_i,x4_i )
  
  x_pred_i <- cbind(x1_pred_i,x4_pred_i )
  
  fun_OLS <- function(par){
    
    phi<-par[1]
    psi<-par[2]
    
    
    beta1<-par[3]
    beta2<-par[4]
    # beta3<-par[4]
    # beta4<-par[5]
    # beta5<-par[6]
    # beta6<-par[7]
    
    
    beta<-rbind(beta1,beta2)
    
    A<-fun_A(phi,T_used)
    
    Sigma<-fun_Sigma(psi,T_used)
    
    Sigma_<-solve(Sigma)
    
    C<-Sigma_%*%l%*%t(l)%*%Sigma_/c(t(l)%*%Sigma_%*%l)
    
    value<-0
    
    
    v_i<-A%*%y_i - x_i%*%beta
    
    value<- value + t(v_i)%*%(Sigma_ - C)%*%v_i
    
    return(c(value))   
    
  }
  
  
  # fun_ineq<-function(par){
  #   beta<-par[1]
  #   phi<-par[2]
  #   psi<-par[3]
  #   
  #   return(c(phi+psi))
  # }
  
  estml<-solnp(c(0, 0,0,0), fun = fun_OLS,
               # ineqfun = fun_ineq, 
               # ineqLB = c(-0.999), 
               # ineqUB = c(0.999),
               LB = c(-0.999,-0.999,-Inf,-Inf), 
               UB = c(0.999,0.999,Inf,Inf)) 
  
  #round(2 * (1 - pnorm(abs(estml$pars/sqrt(diag(solve(estml$hessian)))))),3)   
  
  est_gamma<-estml$pars
  
  est_A<-fun_A(est_gamma[1],T_used)
  
  est_Sigma<-fun_Sigma(est_gamma[2],T_used)
  
  est_Sigma_<-solve(est_Sigma)
  
  c_Sigma<-c(t(l)%*%est_Sigma_%*%l)       
  
  
  est_V<-est_A%*%y_i - x_i%*%matrix(est_gamma[-c(1,2)])
  
  est_mu<- 1/c_Sigma*t(l_T_used)%*%est_Sigma_%*%est_V
 
  
  #######################
  est_V_pred<- fun_A(est_gamma[1],T-T_used)%*%y_pred_i - x_pred_i%*%matrix(est_gamma[-c(1,2)])
  
  
  
  est_pred_Sigma<-fun_Sigma(est_gamma[2],T-T_used)
  
  l_T_pred<-matrix(rep(1,T-T_used))
  
  c_Sigma<-c(t(l_T_pred)%*%solve(est_pred_Sigma)%*%l_T_pred)
  
  est_pred_mu<- 1/c_Sigma*t(l_T_pred)%*%solve(est_pred_Sigma)%*%est_V_pred
  
  est_pred_B<-fun_B(est_gamma[2],T-T_used)
  
  est_pred_U <-solve(est_pred_B)%*%(est_V_pred - est_pred_mu%x%l_T_pred)
  
  est_U_pred_2<- est_pred_U
  
  
  mean_mse_one<- c(mean_mse_one, sqrt(mean(est_pred_U^2)))
  
  
  ############################################
  ###########################################
  
  est_B<-fun_B(est_gamma[2],T_used)
  
  est_U <-solve(est_B)%*%(est_V - est_mu%x%l_T_used)
  
  
  ##########################
  est_omega<-mean((est_U^2))
  est_mu4<-mean((est_U^4))
  est_mu3<-mean((est_U^3))
  
  
  #c(omega)
  
  
  
  ml_GARCH<-function(par){
    
    tau<-par[1]
    nu<-par[2]
    
    
    value<-0
    h0<- est_omega*(1-nu-tau) + nu*est_omega
    value<-value + log(dnorm(est_U[1],sd=sqrt(h0)))
    
    for(t in 2:T_used){
      
      h0<- est_omega*(1-nu-tau) + nu*h0 + tau*est_U[t-1]^2
      
      
      value<-value+log(dnorm(est_U[t],sd=sqrt(h0)))
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
  
  
  print(i)
  
  
  
  
  
}
 


library(Rsolnp)


T<-dim(ymat)[2]

n<-N<-dim(ymat)[1]

# T_used<-96
T_pred<-T - T_used 



l_N<-matrix(rep(1,N))

l_T_used<-l<-l_T<-matrix(rep(1,T_used))

J_N<-diag(N) - l_N%*%t(l_N)/N

J_T_used<-diag(T_used) - l_T_used%*%t(l_T_used)/T_used
 

 
fun_A<-function(phi,T_used = T_used){
  
  A <- diag(T_used)
  diag(A[-1,-T_used]) <- -phi
  
  return(A)
}
############

fun_A_dot1<-function(phi,T_used = T_used){
  
  A_dot1 <- diag(T_used)*0
  diag(A_dot1[-1,-T_used]) <- -1
  return(A_dot1)
}
######
######
fun_B<-function(psi,T_used= T_used){
  
  B <- diag(T_used)
  diag(B[-1,-T_used]) <- psi
  
  return(B)
}
########
fun_B_dot1<-function(psi,T_used= T_used){
  
  B_dot1 <- diag(T_used)*0
  diag(B_dot1[-1,-T_used]) <- 1
  return(B_dot1)
}
###############
fun_Sigma<-function(psi,T_used= T_used){
  B<-fun_B(psi,T_used)
  return(B%*%t(B))
}


#######################################
#######################################


ymat_used <- as.matrix(ymat[,1:T_used])

ymat_pred <- as.matrix(ymat[,(T_used+1):T])

x1<-xmat1[,1:T_used]
# x2<-xmat2[,1:T_used]
# x3<-xmat3[,1:T_used]
x4<-xmat4[,1:T_used]


x1_pred<-xmat1[,(T_used+1):T]
# x2_pred<-xmat2[,(T_used+1):T]
# x3_pred<-xmat3[,(T_used+1):T]
x4_pred<-xmat4[,(T_used+1):T]


Y<-as.vector(t(ymat_used))

X<-cbind(as.vector(t(x1)), as.vector(t(x4)) ) 


Y_pred<-as.vector(t(ymat_pred))

X_pred<-cbind(as.vector(t(x1_pred)), as.vector(t(x4_pred)) ) 


fun_OLS <- function(par){
  
  phi<-par[1]
  psi<-par[2]
  
  
  beta1<-par[3]
  beta2<-par[4]
  # beta3<-par[4]
  # beta4<-par[5]
  # beta5<-par[6]
  # beta6<-par[7]
  
  
  beta<-rbind(beta1,beta2)
  
  A<-fun_A(phi,T_used)
  
  Sigma<-fun_Sigma(psi,T_used)
  
  Sigma_<-solve(Sigma)
  
  C<-Sigma_%*%l%*%t(l)%*%Sigma_/c(t(l)%*%Sigma_%*%l)
  
  
  value<-0
  for(i in 1:N){
    
    y_i<-ymat_used[i,]
    
    
    x_i <-X[((i-1)*T_used+1):(i*T_used),]
    
    
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

estml<-solnp(c(0,0, 0,0), fun = fun_OLS,
             # ineqfun = fun_ineq, 
             # ineqLB = c(-0.999), 
             # ineqUB = c(0.999),
             LB = c(-0.999,-0.999,-Inf,-Inf ), 
             UB = c(0.999,0.999,Inf,Inf )) 

#round(2 * (1 - pnorm(abs(estml$pars/sqrt(diag(solve(estml$hessian)))))),3)   

est_gamma<-estml$pars

est_A<-fun_A(est_gamma[1],T_used)

est_Sigma<-fun_Sigma(est_gamma[2],T_used)

est_Sigma_<-solve(est_Sigma)

c_Sigma<-c(t(l)%*%est_Sigma_%*%l)

est_V<-(diag(N)%x%est_A)%*%Y - X%*%matrix(est_gamma[-c(1,2)])

est_mu<-(diag(N)%x%(1/c_Sigma*t(l)%*%est_Sigma_))%*%est_V


est_V_pred<- (diag(N)%x%fun_A(est_gamma[1],T-T_used))%*%Y_pred - X_pred%*%matrix(est_gamma[-c(1,2)]) 


est_pred_Sigma<-fun_Sigma(est_gamma[2],T-T_used)

l_T_pred<-matrix(rep(1,T-T_used))

c_Sigma<-c(t(l_T_pred)%*%solve(est_pred_Sigma)%*%l_T_pred) 




est_pred_mu<- (diag(N)%x% (1/c_Sigma*t(l_T_pred)%*%solve(est_pred_Sigma)))%*%est_V_pred

est_pred_B<-fun_B(est_gamma[2],T-T_used)

est_pred_U <- (diag(N)%x%solve(est_pred_B))%*% (est_V_pred - est_pred_mu%x%l_T_pred) 

est_U_pred_1 <- est_pred_U[((i-1)*T_pred+1):(i*T_pred)]

sqrt(mean(est_pred_U^2))

mean_mse_our<-c()

for(i in 1:dim(ymat)[1]){
  
  mean_mse_our<-c(mean_mse_our, sqrt(mean((est_pred_U[((i-1)*T_pred+1):(i*T_pred)])^2 )))
  
} 

cbind(mean_mse, mean_mse_our)
############################################
###########################################
 


# dm.test(est_U_pred_1,est_U_pred_2,alternative="less") 
# 
# 
# # dm.test(error_pred_1,error_pred_2,alternative="less") 
# 
# dm.test(var_error_1,var_error_2,alternative="less")  

 

# 
# rm(list = ls())
