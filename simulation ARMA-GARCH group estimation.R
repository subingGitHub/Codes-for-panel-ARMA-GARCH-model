n<-N<-100

T<-100

n_ratio <- c(5, 5)
n_ratio <- n_ratio / sum(n_ratio)

N_1 <- N * n_ratio[1]
N_2 <- N * n_ratio[2]


beta_all <- c(1, 3)


phi_all <- c(0.1, 0.3)


psi_all <- c(0.1, 0.3)

tau_all <- c(0.1, 0.2)
nu_all <- c(0.25, 0.5)



sample_size <- 500



# beta_all <- c(3, 3)

# phi_all <- c(0.3, 0.3)

# psi_all <- c(0.3, 0.3)

# tau_all <- c(0.2, 0.2)
# nu_all <- c(0.5, 0.5) 


mu <- rnorm(N)
omega_all <- matrix(runif(N,1,3))

mu_all <- matrix(mu,nrow=1)
omega_all <- matrix(omega_all,nrow=1)


# se_mu_all <- c()

# se_omega_all <- c()

fun_A<-function(phi,T_used = T){
  
  A <- diag(T_used)
  diag(A[-1,-T_used]) <- -phi
  
  
  return(A)
}
#####

fun_A_dot1<-function(phi, T_used){
  
  A_dot1 <- diag(T_used)*0
  diag(A_dot1[-1,-T_used]) <- -1
  return(A_dot1)
} 

######
fun_B<-function(psi,T_used = T){
  
  B <- diag(T_used)
  diag(B[-1,-T_used]) <- psi
  
  
  return(B)
}
########
fun_B_dot1<-function(psi,T_used = T){
  
  B_dot1 <- diag(T_used)*0
  diag(B_dot1[-1,-T_used]) <- 1
  return(B_dot1)
}
###############
fun_Sigma<-function(psi,T_used = T){
  B<-fun_B(psi,T_used)
  return(B%*%t(B))
}



fun_generate_data <- function(N_, T_, tau_, nu_, phi_, psi_, beta_, mu_, omega_){
    # Generate data for a single subject
       
    U<-c()
    for(i in 1:N_){
        eps0 <- 0
        h0 <- omega_[i]
        eps <- c()
        for(t in 1:(2*T_)){
        h0 <- omega_[i] * (1 - nu_ - tau_) + nu_ * h0 + tau_ * eps0^2
        eps0 <- sqrt(h0) * rnorm(1)
        eps <- rbind(eps, eps0)
        }

    U<-rbind(U,eps)
    }
    
    X_all <- matrix(rnorm(N_ * T_ * 2))
    
    A_phi <- fun_A(phi_, 2*T_) 
    B_psi <- fun_B(psi_, 2*T_) 
    
    Y_all <- (diag(N_) %x% solve(A_phi)) %*% (X_all * beta_ + mu_ %x% matrix(rep(1, 2*T_)) + (diag(N_) %x% B_psi) %*% U)
    
    Y<-c()
    X<-c()
    for(i in 1:N_){
        Y<-cbind(Y,matrix(Y_all[((i-1)*2*T_+1):((i-1)*2*T_+T_)]))
        X<-cbind(X,matrix(X_all[((i-1)*2*T_+1):((i-1)*2*T_+T_)]))
    
    }
    return(list(Y = Y, X = X))
    }





parm_left = c(beta_all[1],phi_all[1],psi_all[1], tau_all[1],nu_all[1])
parm_right = c(beta_all[2],phi_all[2],psi_all[2], tau_all[2],nu_all[2])
estim_left<-matrix(parm_left,nrow=1)
estim_right<-matrix(parm_left,nrow=1)

estim_left_unbiased <- matrix(c(beta_all[1],phi_all[1],psi_all[1], beta_all[1],phi_all[1],psi_all[1], tau_all[1],nu_all[1]), nrow = 1)
estim_right_unbiased <- matrix(c(beta_all[2],phi_all[2],psi_all[2], beta_all[2],phi_all[2],psi_all[2], tau_all[2],nu_all[2]), nrow = 1)

se_left <- c()
se_right <- c()

se_left_unbiased <- c()
se_right_unbiased <- c()


NMI_all <- c()





while(dim(estim_left)[1] < sample_size)
{

res1 <- fun_generate_data(N * n_ratio[1], T, tau_all[1], nu_all[1], phi_all[1], psi_all[1], beta_all[1], mu_all[1: N_1], omega_all[1: N_2])

res2 <- fun_generate_data(N * n_ratio[2], T, tau_all[2], nu_all[2], phi_all[2], psi_all[2], beta_all[2], mu_all[(N_1+1): (N_1 + N_2)], omega_all[(N_1+1): (N_1 + N_2)])


Y <- cbind(res1$Y, res2$Y)
X <- cbind(res1$X, res2$X)


fun_ml_ARMA_with_fixed_effect_onebyone <-  function(par) {
        mu_i<- par[1]  
        beta <- par[2]
        phi <- par[3]
        psi <- par[4]
        omega_i <- par[5]

      
        T <- length(y_i)
        
        A <- fun_A(phi, T)
        B <- fun_B(psi, T)
        
        u_i <- solve(B) %*% (A %*% y_i - x_i * beta - mu_i)
        
        value <- 0
        for (t in 1:T) {
          value <- value + log(dnorm(u_i[t], sd = sqrt(omega_i)))
        }
        
        return(-value)
    }

 
library(parallel)
library(Rsolnp)

n_cores <- 20
cl <- makeCluster(n_cores, type = "PSOCK")


# 每个 worker 加载包
clusterEvalQ(cl, { library(Rsolnp) })

# 导出需要的对象
clusterExport(
  cl,
  varlist = c("Y","X","N","fun_A","fun_B",
              "mu_all","beta_all","phi_all","psi_all","omega_all"),
  envir = environment()
)


# 可选：设置随机数流（重现实验）
# parallel::clusterSetRNGStream(cl, iseed = 123)

fit_one <- function(i) {
  # 取第 i 个体的数据
  y_i <- matrix(Y[, i], ncol = 1)
  x_i <- matrix(X[, i], ncol = 1)
  Tlen <- nrow(y_i)

  # 目标函数（给 solnp）
  obj <- function(par) {
    mu_i   <- par[1]
    beta   <- par[2]
    phi    <- par[3]
    psi    <- par[4]
    omega  <- par[5]

    # 数值保护
    if (!is.finite(omega) || omega <= 0) return(1e12)

    A <- fun_A(phi, Tlen)
    B <- fun_B(psi, Tlen)

    # rhs = A %*% y - x * beta - mu
    rhs <- A %*% y_i - x_i * beta - mu_i

    # u = B^{-1} rhs ；用 solve(B, rhs) 更稳
    u_i <- tryCatch(solve(B, rhs), error = function(e) return(matrix(Inf, nrow = Tlen, ncol = 1)))
    if (!all(is.finite(u_i))) return(1e12)

    # 对数似然（正态）
    # dnorm 的 log=TRUE 更稳；返回负对数似然供最小化
    -sum(dnorm(u_i, sd = sqrt(omega), log = TRUE))
  }

  init  <- c(mu_all[i], beta_all[1], phi_all[1], psi_all[1], omega_all[i])
  lower <- c(-Inf, -Inf, 0, 0, 0.1)
  upper <- c( Inf,  Inf, 1, 1, Inf)

  res <- try(suppressWarnings(
    solnp(pars = init, fun = obj, LB = lower, UB = upper,
          control = list(trace = 0))
  ), silent = TRUE)

  if (inherits(res, "try-error")) {
    return(rep(NA_real_, 5L))  # 失败返回 NA 向量（长度与参数一致）
  } else {
    return(as.numeric(res$pars)) # 成功返回 5 维参数
  }
}
 
res_list <- parLapply(cl, seq_len(N), fit_one)

stopCluster(cl)
 
estimation_result <- do.call(cbind, res_list)




# estimation_result_mean <- apply(estimation_result, 1, mean)

# estimation_result_mean

# beta_estimation_all <- estimation_result[2,]
# phi_estimation_all <- estimation_result[3,]
# psi_estimation_all <- estimation_result[4,]
# 
# df <- data.frame(
#   index = seq_along(beta_estimation_all),  # 生成 1,2,3,... 的序号
#   value = beta_estimation_all
# )
# 
# # 按数值排序
# df_sorted <- df[order(df$value), ]   
#  
# rownames(df_sorted) <- NULL
# sorted_values <- df_sorted$value
# 
# 
# begin_index <- 1
# end_index <- (nrow(df_sorted))
# 
# 
#  
# # 计算给定 mid_index 的 Δ
# fun_get_delta_mid <- function(mid_index, begin_index = 1, end_index = length(sorted_values), sorted_values) {
#   n_seg <- end_index - begin_index + 1
#   nL <- mid_index - begin_index + 1
#   nR <- end_index - mid_index
# 
#   mL <- mean(sorted_values[begin_index:mid_index])
#   mR <- mean(sorted_values[(mid_index+1):end_index])
# 
#   delta_mid <- sqrt(nL * nR / n_seg) * abs(mR - mL)
#   return(delta_mid)
# }
# 
# # 在区间 [begin_index, end_index] 扫一遍所有可能的 mid_index
# find_best_k_loop <- function(sorted_values, begin_index = 1, end_index = length(sorted_values), min_size = 1) {
#   results <- data.frame(mid_index = integer(0), dval = numeric(0))
#   
#   for (mid_index in (begin_index + min_size - 1):(end_index - min_size)) {
#     dval <- fun_get_delta_mid(mid_index, begin_index, end_index, sorted_values)
#     results <- rbind(results, data.frame(mid_index = mid_index, dval = dval))
#   }
#   
#   # 按 dval 从大到小排序
#   results <- results[order(-results$dval), ]
#   
#   # 返回最优解 + 排序表
#   list(
#     k_hat      = results$mid_index[1],
#     delta_hat  = results$dval[1],
#     all_sorted = results
#   )
# }
#  
# # 假设 df_sorted 已经有列 index / value
# b <- df_sorted$value
# res <- find_best_k_loop(b, begin_index = 1, end_index = length(b), min_size = 2)
# 
# k_hat <- res$k_hat
# 
# left_group  <- df_sorted$index[1:k_hat]
# right_group <- df_sorted$index[(k_hat+1):nrow(df_sorted)]
# 
# 
# # get NMI
# # here A1/2 is estimation, and B1/2 is truth
# num_A1 = length(left_group)
# num_A2 = length(right_group)
# 
# num_B1 = N_1
# num_B2 = N_2
# 
# num_A1_B1 = sum(left_group <= N_1)
# num_A1_B2 = sum(left_group > N_1)
# num_A2_B1 = sum(right_group <= N_1)
# num_A2_B2 = sum(right_group > N_1)
# 
# I_AB = 0
# for (count in c(num_A1_B1, num_A1_B2, num_A2_B1, num_A2_B2)) {
#   if (count == 0) next  # 避免 log(0)
#   # 计算每一项并累加
#   I_AB <- I_AB + (count / N) * log((N * count) / ( (if (count %in% c(num_A1_B1, num_A1_B2)) num_A1 else num_A2) * (if (count %in% c(num_A1_B1, num_A2_B1)) num_B1 else num_B2) ))
# }
#  
# H_A = - (num_A1 / N) * log(num_A1 / N) - (num_A2 / N) * log(num_A2 / N)
# H_B = - (num_B1 / N) * log(num_B1 / N) - (num_B2 / N) * log(num_B2 / N)
# NMI = I_AB / ((H_A + H_B) / 2)
# 
# 
# 
# Y_left  <- Y[, left_group, drop = FALSE]  
# Y_right <- Y[, right_group, drop = FALSE]  
# 
# X_left  <- X[, left_group, drop = FALSE]  
# 
# X_right <- X[, right_group, drop = FALSE]


##########################



# fun_get_delta_mid <- function(mid_index, begin_index = 1, end_index = length(sorted_values), sorted_values) {
#   n_seg <- end_index - begin_index + 1
#   nL <- mid_index - begin_index + 1
#   nR <- end_index - mid_index
#   mL <- mean(sorted_values[begin_index:mid_index])
#   mR <- mean(sorted_values[(mid_index+1):end_index])
#   sqrt(nL * nR / n_seg) * abs(mR - mL)
# }
# 
# find_best_k_loop <- function(sorted_values, begin_index = 1, end_index = length(sorted_values), min_size = 1) {
#   if ((end_index - begin_index + 1) < 2 * min_size)
#     return(list(k_hat = NA_integer_, delta_hat = -Inf, all_sorted = NULL))
#   results <- data.frame(mid_index = integer(0), dval = numeric(0))
#   for (mid_index in (begin_index + min_size - 1):(end_index - min_size)) {
#     dval <- fun_get_delta_mid(mid_index, begin_index, end_index, sorted_values)
#     results <- rbind(results, data.frame(mid_index = mid_index, dval = dval))
#   }
#   if (nrow(results) == 0)
#     return(list(k_hat = NA_integer_, delta_hat = -Inf, all_sorted = NULL))
#   results <- results[order(-results$dval), ]
#   list(k_hat = results$mid_index[1],
#        delta_hat = results$dval[1],
#        all_sorted = results)
# }
# 
# # ---- 多变点检测 + 返回每个终端区间的“原始 index 列表” ----
# detect_change_points <- function(sorted_values,
#                                  orig_index,        # 与 sorted_values 对齐的原始 index 向量（如 df_sorted$index）
#                                  delta_thresh,
#                                  min_size = 2,
#                                  begin_index = 1,
#                                  end_index   = length(sorted_values),
#                                  max_cps = Inf) {
#   stopifnot(length(sorted_values) == length(orig_index))
#   
#   cps <- integer(0)
#   deltas <- numeric(0)
#   segs <- list()
#   leaf_groups <- list()  # 每个叶子段的“原始 index”列表
#   
#   # 递归
#   rec <- function(L, R) {
#     # 段太短：作为叶子段收集
#     if ((R - L + 1) < 2 * min_size) {
#       leaf_groups[[length(leaf_groups) + 1]] <<- orig_index[L:R]
#       return(invisible(NULL))
#     }
#     best <- find_best_k_loop(sorted_values, L, R, min_size)
#     k <- best$k_hat; d <- best$delta_hat
#     
#     # 无可切点/不显著/达到上限：作为叶子段收集
#     if (is.na(k) || d < delta_thresh || length(cps) >= max_cps) {
#       leaf_groups[[length(leaf_groups) + 1]] <<- orig_index[L:R]
#       return(invisible(NULL))
#     }
#     
#     # 记录一次分割
#     cps <<- c(cps, k)
#     deltas <<- c(deltas, d)
#     segs[[length(segs) + 1]] <<- data.frame(L = L, K = k, R = R, delta = d)
#     
#     # 分别递归
#     rec(L, k)
#     rec(k + 1, R)
#   }
#   
#   rec(begin_index, end_index)
#   
#   # 输出
#   o <- order(cps)
#   cps <- cps[o]; deltas <- deltas[o]
#   list(
#     cps_indices   = cps,                         # 变点位置（在已排序序列中的下标）
#     cps_values    = if (length(cps)) sorted_values[cps] else numeric(0),
#     cps_deltas    = deltas,
#     segments_info = if (length(segs)) do.call(rbind, segs) else NULL,
#     leaf_groups   = leaf_groups                  # << 每个终端段对应的“原始 index”向量的 list
#   )
# }
# 
# 
# get_BIC_res <- function(Y_mat, X_mat){
#   
#   T_used <- nrow(Y_mat)
#   N_used <- ncol(Y_mat)
#   l <- matrix(1, nrow = T_used, ncol = 1)
#   # 生成 OLS 目标函数：给定 (Y_mat, X_mat, l)，返回一个以 par = (beta, phi, psi) 为自变量的函数
#   make_obj <- function(Y_mat, X_mat) {
#     
#     stopifnot(nrow(Y_mat) == nrow(X_mat))
#     
#     
#     function(par) {
#       beta <- par[1]; phi <- par[2]; psi <- par[3]
#       
#       A     <- fun_A(phi, T_used)           # 维度 T_used×T_used
#       Sigma <- fun_Sigma(psi, T_used)       # 维度 T_used×T_used（对称正定假设）
#       # 数值更稳的“逆”
#       Sigma_inv <- chol2inv(chol(Sigma))
#       
#       denom <- as.numeric(t(l) %*% Sigma_inv %*% l)
#       C     <- Sigma_inv %*% l %*% t(l) %*% Sigma_inv / denom
#       M     <- Sigma_inv - C        # 公共矩阵
#       
#       # V 的列就是各 i 的 v_i：V = A Y_mat - beta X_mat   (T_used×N_used)
#       V  <- A %*% Y_mat - beta * X_mat
#       MV <- M %*% V
#       
#       # sum_i v_i' M v_i = trace(V' M V) = sum(colSums(V * (M V)))
#       value <- sum(colSums(V * MV))
#       return(value)
#     }
#   }
#   
#   obj <- make_obj(Y_mat, X_mat)
#   mean_fit <- try(
#     suppressWarnings(
#       solnp(pars = c(beta_all[1], phi_all[1], psi_all[1]),
#             fun  = obj,
#             LB   = c(-Inf, -0.999, -0.999),
#             UB   = c( Inf,  0.999,  0.999),
#             control = list(trace = 0)),  # 关闭迭代输出)
#     ),
#     silent = TRUE
#   )
#   if (inherits(mean_fit, "try-error")) return(NULL)
#   
#   return(tail(mean_fit$value, 1))
# }
# 
# 
# 
# 
# delta_thresh_list <- seq(0.5, 3.0, by = 0.5)
# 
# final_used <- list()
# BIC_min <- 1e15
# 
# for (delta_thresh in delta_thresh_list) {
#   print(delta_thresh)
#   
#   fit <- detect_change_points(
#     sorted_values = df_sorted$value,
#     orig_index    = df_sorted$index,
#     delta_thresh  = delta_thresh,
#     min_size      = 2
#   )
#   
#   BIC_value <- 0
#   for (i in 1:length(fit$leaf_groups)) {
#     Y_i <- Y[, fit$leaf_groups[[i]], drop = FALSE]
#     X_i <- X[, fit$leaf_groups[[i]], drop = FALSE]
#     BIC_value <- BIC_value + get_BIC_res(Y_i, X_i)
#   }
#   BIC_ <- log(BIC_value / (ncol(Y) * nrow(Y))) + log(ncol(Y) * nrow(Y)) * length(fit$leaf_groups) * 3
#   
#   if (BIC_ < BIC_min) {
#     BIC_min <- BIC_
#     final_used <- list(fit = fit, delta_thresh = delta_thresh, BIC = BIC_)
#   }
# }
# 
# 
# final_used$fit
# 
# 
# # 例如把每段取出来构成数据框列表：
# leaf_dfs <- lapply(final_used$fit$leaf_groups, function(ix) df[df$index %in% ix, ])


############################################
############################################
# Y_mat <- Y_left
# X_mat <- X_left


get_estimation_res <- function(Y_mat, X_mat){

T_used <- nrow(Y_mat)
N_used <- ncol(Y_mat)
l <- matrix(1, nrow = T_used, ncol = 1)
# 生成 OLS 目标函数：给定 (Y_mat, X_mat, l)，返回一个以 par = (beta, phi, psi) 为自变量的函数
make_obj <- function(Y_mat, X_mat) {

  stopifnot(nrow(Y_mat) == nrow(X_mat))


  function(par) {
    beta <- par[1]; phi <- par[2]; psi <- par[3]

    A     <- fun_A(phi, T_used)           # 维度 T_used×T_used
    Sigma <- fun_Sigma(psi, T_used)       # 维度 T_used×T_used（对称正定假设）
    # 数值更稳的“逆”
    Sigma_inv <- chol2inv(chol(Sigma))

    denom <- as.numeric(t(l) %*% Sigma_inv %*% l)
    C     <- Sigma_inv %*% l %*% t(l) %*% Sigma_inv / denom
    M     <- Sigma_inv - C        # 公共矩阵

    # V 的列就是各 i 的 v_i：V = A Y_mat - beta X_mat   (T_used×N_used)
    V  <- A %*% Y_mat - beta * X_mat
    MV <- M %*% V

    # sum_i v_i' M v_i = trace(V' M V) = sum(colSums(V * (M V)))
    value <- sum(colSums(V * MV))
    return(value)
  }
}

obj <- make_obj(Y_mat, X_mat)
mean_fit <- try(
  suppressWarnings(
    solnp(pars = c(beta_all[1], phi_all[1], psi_all[1]),
          fun  = obj,
          LB   = c(-Inf, -0.999, -0.999),
          UB   = c( Inf,  0.999,  0.999),
        control = list(trace = 0)),  # 关闭迭代输出)
  ),
  silent = TRUE
)
if (inherits(mean_fit, "try-error")) return(NULL)


est_gamma <- mean_fit$pars



fun_get_estimates_U_omega <- function(Y_mat, X_mat, estml) {
 
 
  est_gamma<-estml$pars

  est_A<-fun_A(est_gamma[2], T_used)

  est_Sigma<-fun_Sigma(est_gamma[3], T_used)

  est_Sigma_<-solve(est_Sigma)

  c_Sigma<-c(t(l)%*%est_Sigma_%*%l)

  est_V<-(diag(N_used)%x%est_A)%*%Y_mat - X_mat*est_gamma[1]

  est_mu<-(diag(N_used)%x%(1/c_Sigma*t(l)%*%est_Sigma_))%*%est_V

  ############################################
  ###########################################

  est_B<-fun_B(est_gamma[3])

  est_U <- (diag(N_used)%x%solve(est_B)) %*% (est_V - est_mu%x%l) 

  est_U <- matrix(est_U, ncol = N_used)
  ###########################################
  est_omega<- apply(est_U, 2, function(u_i) mean(u_i^2))
  est_mu4 <- apply(est_U, 2, function(u_i) mean(u_i^4))
  est_mu3 <- apply(est_U, 2, function(u_i) mean(u_i^3))
  est_cov1 <- apply(est_U, 2, function(u_i) {
    u_i2 <- u_i[2:T_used]^2
    u_i2_lag <- u_i[1:(T_used-1)]^2
    mean(u_i2 * u_i2_lag)
  })
  est_cov2 <- apply(est_U, 2, function(u_i) {
    u_i2 <- u_i[3:T_used]^2
    u_i2_lag <- u_i[1:(T_used-2)]^2
    mean(u_i2 * u_i2_lag)
  })


 
return(list(est_U = est_U, est_omega = est_omega, est_mu = est_mu, est_mu4 = est_mu4, est_mu3 = est_mu3, est_cov1 = est_cov1, est_cov2 = est_cov2))

}

 
estim_U_omega <- fun_get_estimates_U_omega(as.vector(Y_mat), as.vector(X_mat),mean_fit)
est_U <- estim_U_omega$est_U
est_omega <- estim_U_omega$est_omega
est_mu3 <- estim_U_omega$est_mu3
est_mu4 <- estim_U_omega$est_mu4
est_cov1 <- estim_U_omega$est_cov1
est_cov2 <- estim_U_omega$est_cov2
est_mu <- estim_U_omega$est_mu



# ------- 生成 GARCH 目标函数 --------
make_ml_garch <- function(est_U, est_omega) {
 
 
  # 目标函数：返回 -logLik
  obj <- function(par) {
    tau <- par[1]; nu <- par[2]
    val <- 0

    for (i in 1:N_used) {
      # 初始条件 h0
      h0 <- est_omega[i] * (1 - nu - tau) + nu * est_omega[i] 
      h0 <- max(h0, 1e-12)

      # t = 1
      val <- val + dnorm(est_U[1L, i], sd = sqrt(h0), log = TRUE)

      # t = 2..T_used
      for (t in 2: dim(est_U)[1]) {
        h0 <- est_omega[i] * (1 - nu - tau) + nu * h0 + tau * est_U[t - 1L, i]^2
        h0 <- max(h0, 1e-12)
        val <- val + dnorm(est_U[t, i], sd = sqrt(h0), log = TRUE)
      }
    }
    return(-val)  # solnp 最小化
  }

  # 不等式约束：nu + tau ∈ [0.001, 0.999]
  ineq <- function(par) c(par[1] + par[2])

  list(obj = obj, ineq = ineq)
}



var_fit <- try(
    solnp(pars = c(tau_all[1], nu_all[1]),
          fun  = make_ml_garch(est_U, est_omega)$obj,
          ineqfun = make_ml_garch(est_U, est_omega)$ineq,
          ineqLB  = c(0.001),
          ineqUB  = c(0.999),
          LB = c(0.001, 0.001),
          UB = c(0.999, 0.999),
          control = list(trace = 0)   # 不输出迭代过程
    ),
    silent = TRUE
) 
if (inherits(var_fit, "try-error")) return(NULL)


est_zeta <- var_fit$pars


################################################################
################################################################

# mean part
Y_used <- as.vector(Y_mat)
X_used <- as.vector(X_mat)


fir_deri_function<-function(par){

beta<-par[1]
phi<-par[2]
psi<-par[3]


A<-fun_A(phi, T_used)

B<-fun_B(psi, T_used)

A_dot1<-fun_A_dot1(phi, T_used)
B_dot1<-fun_B_dot1(psi, T_used)

Sigma<-fun_Sigma(psi, T_used)

Sigma_<-solve(Sigma)

Sigma_dot1<-B%*%t(B_dot1)+B_dot1%*%t(B)

C<-Sigma_%*%l%*%t(l)%*%Sigma_/c(t(l)%*%Sigma_%*%l)


exp_beta <- 0

exp_phi <- - 2* sum(diag( (diag(est_omega)%x%(t(B)%*%C%*%A_dot1%*%solve(A)%*%B))))

exp_psi <- 2* sum(diag( (diag(est_omega)%x%(t(B)%*%Sigma_%*%Sigma_dot1%*%C%*%B)))) - sum(diag(Sigma_dot1%*%C)) * sum(diag( (diag(est_omega)%x% (t(B)%*%C%*%B))))


return(matrix(c(exp_beta,exp_phi,exp_psi)))

}


 
sec_deri_function<-function(par){

beta<-par[1]
phi<-par[2]
psi<-par[3]


A<-fun_A(phi, T_used)

B<-fun_B(psi, T_used)

A_dot1<-fun_A_dot1(phi, T_used)
B_dot1<-fun_B_dot1(psi, T_used)

Sigma<-fun_Sigma(psi, T_used)

Sigma_ <- solve(Sigma)

Sigma_dot1<-B%*%t(B_dot1)+B_dot1%*%t(B)

B11<-B12<-B22<-diag(T_used)*0
Sigma_ddot11 <- B_dot1%*%t(B_dot1) + B%*%t(B11)+B11%*%t(B) + B_dot1%*%t(B_dot1)

C<-Sigma_%*%l%*%t(l)%*%Sigma_/c(t(l)%*%Sigma_%*%l)



D <- 2*Sigma_%*%Sigma_dot1%*%C - Sigma_%*%Sigma_dot1%*%Sigma_ - sum(diag(Sigma_dot1%*%C))*C


M<- t(solve(A))%*%t(A_dot1)%*%(Sigma_ - C)%*%A_dot1%*%solve(A)


E<- Sigma_%*%Sigma_ddot11%*%Sigma_ - 2*Sigma_%*%Sigma_dot1%*%Sigma_%*%Sigma_dot1%*%Sigma_ + sum(diag(Sigma_ddot11%*%C))*C +
    2*(sum(diag(Sigma_dot1%*%C))*sum(diag(Sigma_dot1%*%C))*diag(T_used) - sum(diag(Sigma_dot1%*%Sigma_%*%Sigma_dot1%*%C))*diag(T_used)  - sum(diag(Sigma_dot1%*%C))*Sigma_%*%Sigma_dot1 + Sigma_%*%Sigma_dot1%*%Sigma_%*%Sigma_dot1 -
        Sigma_%*%Sigma_ddot11 + 2* Sigma_%*%Sigma_dot1%*%Sigma_%*%Sigma_dot1 -  sum(diag(Sigma_dot1%*%C))*Sigma_%*%Sigma_dot1)%*%C


exp_beta_beta <- 2*t(X_used)%*%( diag(N_used) %x% (Sigma_ - C) )%*%X_used

exp_beta_phi <- -2* t( X_used*beta + est_mu%x%l)%*% ( diag(N_used)%x% (t(solve(A))%*%t(A_dot1)%*%(Sigma_ - C)) )%*%X_used

exp_beta_psi <- -2* ( (t(est_mu)%*%diag(N_used))%x% (t(l)%*%D))%*%X_used



exp_phi_phi <- 2* t(X_used*beta + est_mu%x%l) %*% (diag(N_used) %x% M)%*%(X_used*beta + est_mu%x%l) + 2*sum(diag((diag(est_omega)%x%(t(B)%*%M%*%B))))


exp_phi_psi <- 2* ( (t(est_mu)%*%diag(N_used))%x% (t(l)%*%D%*%A_dot1%*%solve(A)))%*%(X_used*beta + est_mu%x%l) + 2*sum(diag( (diag(est_omega)%x%(t(B)%*%D%*%A_dot1%*%solve(A)%*%B))))


exp_psi_psi <- - (t(est_mu)%*%est_mu)*(t(l)%*%E%*%l) - sum(est_omega)*sum(diag(t(B)%*%E%*%B))

# exp_phi_phi <- 2* t(X_used*beta + est_mu%x%l) %*% (diag(N_used) %x% M)%*%(X_used*beta + est_mu%x%l) + 2*sum(diag((diag(N_used)%x%(t(B)%*%M%*%B))%*%Sigma_U))
#
#
# exp_phi_psi <- 2* ( (t(est_mu)%*%diag(N_used))%x% (t(l)%*%D%*%A_dot1%*%solve(A)))%*%(X_used*beta + est_mu%x%l) + 2*sum(diag( (diag(N_used)%x%(t(B)%*%D%*%A_dot1%*%solve(A)%*%B))%*%Sigma_U))
#
#
# exp_psi_psi <- (t(est_mu)%*%diag(N_used)%*%mu)*(t(l)%*%E%*%l) + sum(diag( (diag(N_used)%x%(t(B)%*%E%*%B))%*%Sigma_U))
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
for(t in 1:T_used){
    for(t_star in 1:T_used){
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
for(t in 1:T_used){
    for(t_star in 1:T_used){
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

lag = 4


se_lambda<-function(par){

beta<-par[1]
phi<-par[2]
psi<-par[3]
tau<-par[4]
nu<-par[5]

A<-fun_A(phi, T_used)
B<-fun_B(psi, T_used)

A_dot1<-fun_A_dot1(phi, T_used)
B_dot1<-fun_B_dot1(psi, T_used)

Sigma_dot1<-B%*%t(B_dot1)+B_dot1%*%t(B)

Sigma<-fun_Sigma(psi, T_used)

Sigma_<-solve(Sigma)

C<-Sigma_%*%l%*%t(l)%*%Sigma_/c(t(l)%*%Sigma_%*%l)

D <- 2*Sigma_%*%Sigma_dot1%*%C - Sigma_%*%Sigma_dot1%*%Sigma_ - sum(diag(Sigma_dot1%*%C))*C

V<-(diag(N_used)%x%A)%*%Y_used - X_used*beta

se_beta_beta<- 4* t(X_used)%*%(diag(est_omega)%x% (Sigma_ - C))%*%X_used

se_phi_phi<-0

se_psi_psi<-0

se_beta_phi<-0

se_beta_psi<-0

se_phi_psi<-0

for(i in 1:N_used){
    est_ui<-est_U[((i-1)*T_used+1):(i*T_used)]
    
    e2<-c()
    for(t in 1:T_used){
    if(t<lag){
        e2<-c(e2, e2_fun(est_ui, t))
    }
    else{
        e2<-c(e2,0)
    }
    
    }
    
    
    x<-X_used[((i-1)*T_used+1):(i*T_used)]
    M_phi<- t(B)%*%(Sigma_ - C)%*%A_dot1%*%solve(A)%*%B
    b_phi<- t(B)%*%(Sigma_ - C)%*%A_dot1%*%solve(A)%*%(est_mu[i]*l + x*beta)
    
    M_psi<- t(B)%*%D%*%B
    b_psi<- t(B)%*%(D+t(D))%*%l*est_mu[i]
    
    
    se_phi_phi<-se_phi_phi + 4 * fun_cov(c(est_omega[i],est_mu4[i]), e2, M_phi, b_phi)
    
    mu_mu_ldl<-0#est_mu[i]^2*c(t(l)%*%D%*%l)
    se_psi_psi<-se_psi_psi + fun_cov(c(est_omega[i],est_mu4[i]), e2, M_psi, b_psi) #+ # mu_mu_ldl^2 + 2*est_omega[i]*sum(diag(M_psi))*mu_mu_ldl
    
    ########################
    ########################
    
    b_beta_phi_1<-t(B)%*%(Sigma_  - C)%*%x
    b_beta_phi_2<-t(B)%*%(Sigma_  - C)%*%A_dot1%*%solve(A)%*%(est_mu[i]*l + x*beta)
    
    se_beta_phi<-se_beta_phi - 4*sum(b_beta_phi_1*b_beta_phi_2)*est_omega[i]
    
    
    
    b_beta_psi_1<- t(B)%*%(Sigma_  - C)%*%x
    b_beta_psi_2<- t(B)%*%(D+t(D))%*%l*est_mu[i]
    
    se_beta_psi<-se_beta_psi - 2*sum(b_beta_psi_1*b_beta_psi_2)*est_omega[i]
    
    
    ########################
    ########################
    
    b_phi_psi_1<- t(B)%*%(Sigma_ - C)%*%A_dot1%*%solve(A)%*%(est_mu[i]*l + x*beta)
    b_phi_psi_2<-  t(B)%*%(D+t(D))%*%l*est_mu[i]
    
    M_phi_psi_1<-t(B)%*%(Sigma_)%*%A_dot1%*%solve(A)%*%B
    
    M_phi_psi_2<-t(B)%*%D%*%B
    UMUUMU<-(est_mu4[i] - 3*est_omega[i]^2)*sum(diag(M_phi_psi_1)*diag(M_phi_psi_2)) + fun_mix_a1(e2, M_phi_psi_1, M_phi_psi_2)  + est_omega[i]^2*(sum(diag(M_phi_psi_1%*%t(M_phi_psi_2)))+sum(diag(M_phi_psi_1%*%M_phi_psi_2)) )
    
    se_phi_psi<- se_phi_psi  + 2*sum(b_phi_psi_1*b_phi_psi_2)*est_omega[i] + 2* UMUUMU #+ 2*sum(diag(M_phi_psi_1))*est_omega[i]*est_mu[i]^2*c(t(l)%*%D%*%l)
    
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
for(t in 2:T_used){

    h0<- omega*(1-nu-tau) + nu*h0 + tau*y[t-1]^2

    h<-c(h,h0)
}

return(h)

}


eps_est<-c()

for(r in 1:N_used){

u_est1_<-est_U[((r-1)*T_used+1):(r*T_used)]


eps_est<-c(eps_est,diag(1/sqrt(fun_h(c(est_omega[r],est_zeta),u_est1_)))%*%u_est1_)

}

mu4<-mean(eps_est^4)



fun_h_nu<-function(par,y){
omega<-par[1]
tau<-par[2]
nu<-par[3]


h0<-0
h<-c(h0)
for(t in 2:T_used){

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
for(t in 2:T_used){

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

for(r in 1:N_used){

    u_est1_<-est_U[((r-1)*T_used+1):(r*T_used)]

    x<-X_used[((r-1)*T_used+1):(r*T_used)]

    y<-Y_used[((r-1)*T_used+1):(r*T_used)]

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
    h0_phi<- -(1-tau-nu)* mean(2*(u_est1_star[-1])*(y_star[-T_used]))
    h0_psi<- -(1-tau-nu)* mean(2*(u_est1_star[-1])*(u_est1_star[-T_used]))

    h_beta<-c(h0_beta)
    h_phi<-c(h0_phi)
    h_psi<-c(h0_psi)

    h0_beta <- - (1-tau-nu)* mean(2*(u_est1_star)*x_star) + nu*h0_beta - tau*2*u_est1_star[1]*x_star[1]

    h0_phi <- -(1-tau-nu)* mean(2*(u_est1_star[-1])*(y_star[-T_used])) + nu*h0_phi

    h0_psi <- -(1-tau-nu)* mean(2*(u_est1_star[-1])*(u_est1_star[-T_used])) + nu*h0_psi

    h_beta<-c(h_beta,h0_beta)
    h_phi<-c(h_phi,h0_phi)
    h_psi<-c(h_psi,h0_psi)


    for(t in 3:T_used){

    h0_beta <- - (1-tau-nu)* mean(2*(u_est1_star)*x_star) + nu*h0_beta - tau*2*u_est1_star[t-1]*x_star[t-1]


    h0_phi <- -(1-tau-nu)* mean(2*(u_est1_star[-1])*(y_star[-T_used])) + nu*h0_phi - tau*2*u_est1_star[t-1]*y_star[t-2]

    h0_psi <- -(1-tau-nu)* mean(2*(u_est1_star[-1])*(u_est1_star[-T_used])) + nu*h0_psi - tau*2*u_est1_star[t-1]*u_est1_star[t-2]

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
for(r in 1:N_used){

    u_est1_<-est_U[((r-1)*T_used+1):(r*T_used)]

    h_est<-fun_h(c(est_omega[r],zeta),u_est1_)

    h_tau<-fun_h_tau(c(est_omega[r],zeta),u_est1_)
    h_nu<-fun_h_nu(c(est_omega[r],zeta),h_est)

    Pi_1<-1/h_est*h_tau - mean(1/h_est^2*h_tau)*h_est
    Pi_2<-1/h_est*h_nu - mean(1/h_est^2*h_nu)*h_est



    O_11<-O_11 + sum(Pi_1*Pi_1)
    O_12<-O_12 + sum(Pi_1*Pi_2)

    O_22<-O_22 + sum(Pi_2*Pi_2)

    # J_1_11<-J_1_11+ sum(  ( 2*u_est1_[-1]^2/(h_est[-1])^3 - 1/(h_est[-1])^2 ) *   (u_est1_[-T_used]^2  )^2)
    # J_1_12<-J_1_12+ sum(( 2*u_est1_[-1]^2/(h_est[-1])^3 - 1/(h_est[-1])^2 ) * (u_est1_[-T_used]^2 )*(h_est[-T_used] ) )
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
for(r in 1:N_used){

    u_est1_<-est_U[((r-1)*T_used+1):(r*T_used)]

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

for(r in 1:N_used){

    u_est1_<-est_U[((r-1)*T_used+1):(r*T_used)]

    h_est<-fun_h(c(est_omega[r],zeta),u_est1_)

    h_tau_est<-fun_h_tau(c(est_omega[r],zeta),u_est1_)
    h_nu_est<-fun_h_nu(c(est_omega[r],zeta),h_est)

    exp_uh<- mu4-1

    #mean((u_est1_^2 - h_est)^2/h_est^2


    omega_est_omega<--1/(1-tau-nu)*((1-nu)*mean(u_est1_^2-h_est)  - 1/T_used*(tau+nu)*u_est1_[T_used]^2)

    sum_1<-sum( (h_tau_est) / (h_est)^2 )

    L_2<- 1/(1-tau-nu)*(tau+nu)/(1-nu)*mean((u_est1_^2-h_est)/h_est^2)*u_est1_[T_used]^2

    #

    c_1<- mean(h_tau_est / h_est^2)

    M_25<- -tau*mean(u_est1_^2)*( mean(u_est1_[-T_used]^2/(h_est[-1])^2)
                                + nu^2*mean(u_est1_[-c(T_used,T_used-1)]^2/(h_est[-c(1,2)])^2)
                                + nu^4*mean(u_est1_[-c(T_used,T_used-1,T_used-2)]^2/(h_est[-c(1,2,3)])^2)
                                + nu^6*mean(u_est1_[-c(T_used,T_used-1,T_used-2,T_used-3)]^2/(h_est[-c(1,2,3,4)])^2))

    M_34<- mean(h_tau_est / (h_est)^2 )*mean(u_est1_^2)
    #mean( 1/tau* (h_est - est_omega[r]) / (h_est)^2)*mean(u_est1_^2)


    J_1<-J_1 -  1/(1-tau-nu) * exp_uh + L_2 + (1-tau-nu)/(1-nu)*omega_est_omega*sum_1 - c_1*sum(u_est1_^2-h_est) + M_25 +  M_34


    L_2_2 <- 1/(1-tau-nu)*tau/(1-nu)^2*(tau+nu)*mean((u_est1_^2-h_est)/h_est^2)*u_est1_[T_used]^2

    sum_2<-sum(h_nu_est / (h_est)^2 )


    c_2<- mean(h_nu_est / (h_est)^2 )

    M_25_2 <- -tau^2*mean(u_est1_^2)*(nu*mean(u_est1_[-c(T_used,T_used-1)]^2/(h_est[-c(1,2)])^2)
                                    + 2*nu^2*mean(u_est1_[-c(T_used,T_used-1,T_used-2)]^2/(h_est[-c(1,2,3)])^2)
                                    + 3*nu^3*mean(u_est1_[-c(T_used,T_used-1,T_used-2,T_used-3)]^2/(h_est[-c(1,2,3,4)])^2)
                                    + 4*nu^4*mean(u_est1_[-c(T_used,T_used-1,T_used-2,T_used-3,T_used-4)]^2/(h_est[-c(1,2,3,4,5)])^2)
                                    + 5*nu^5*mean(u_est1_[-c(T_used,T_used-1,T_used-2,T_used-3,T_used-4,T_used-5)]^2/(h_est[-c(1,2,3,4,5,6)])^2))


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
A<-fun_A(phi, T_used)

B<-fun_B(psi, T_used)

A_dot1<-fun_A_dot1(phi, T_used)
B_dot1<-fun_B_dot1(psi, T_used)

Sigma<-fun_Sigma(psi, T_used)

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
for(i in 1:N_used){
    
    x<-X_used[((i-1)*T_used+1):(i*T_used)]
    y<-Y_used[((i-1)*T_used+1):(i*T_used)]
    
    
    v<-A%*%y - x*beta
    
    u <- solve(B)%*%(v - est_mu[i]) 
    
    
    tildeQ_beta<- -2*diag(c(u))%*%(Sigma_ - C)%*%x
    
    tildeQ_phi<- 2*diag(c(u))%*%lower.tri(t(B)%*%(Sigma_ - C)%*%A_dot1%*%solve(A)%*%B)%*%u + 2*diag(c(u))%*%t(B)%*%(Sigma_ - C)%*%A_dot1%*%solve(A)%*%B%*%(est_mu[i] + x*beta)
    
    tildeQ_psi<- diag(c(u))%*%lower.tri(t(B)%*%D%*%B)%*%u + 2 * diag(c(u))%*%t(B)%*%D%*%l* est_mu[i]
    
    # diag(c(v))%*%D%*%v - 2 * est_omega[i]*(sum(diag((t(B)%*%Sigma_%*%Sigma_dot1%*%C%*%B))) - sum(diag(Sigma_dot1%*%C)) * sum(diag((t(B)%*%C%*%B))))
    
    u_est1_<-est_U[((i-1)*T_used+1):(i*T_used)]
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




return(M_est/N_used/T_used)
}
 
##########################
#SE-mean

Gamma_<--solve(sec_deri_function(est_gamma))
Lambda<-se_lambda(c(est_gamma,est_zeta))

se_gamma_ <- sqrt(diag(Gamma_%*%Lambda%*%Gamma_))

 
Gamma_zeta_<--solve(sec_deri_function_zeta(est_zeta))

Omega_zeta<-Omega_function(est_zeta)

Pi_zeta <- Pi_function(est_zeta)


Lambda_zeta<- Lambda_function(c(est_gamma,est_zeta))
#
Xi<-(mu4-1)*Omega_zeta + Pi_zeta%*%Gamma_%*%Lambda%*%Gamma_%*%t(Pi_zeta) + Pi_zeta%*%Gamma_%*%Lambda_zeta + t(Pi_zeta%*%Gamma_%*%Lambda_zeta)
#
#
#
 

se_zeta_ <- sqrt(diag(Gamma_zeta_%*%Xi%*%Gamma_zeta_))


##################################
##################################
# analytical bias correction
est_gamma_unbias<-est_gamma+solve(sec_deri_function(est_gamma))%*%fir_deri_function(est_gamma)


#################################
# Jacknife_correction

 

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
    for(i in 1:dim(Y_mat)[2]){
      
      y_i<-Y_mat[(1:(T/2)),i]
      
      x_i<-X_mat[(1:(T/2)),i]
      v_i<-A%*%y_i - x_i*beta
      
      value<- value + t(v_i)%*%(Sigma_ - C)%*%v_i
    }
    
    return(value)   
    
  }


  estml<-try(suppressWarnings(solnp(c(3,0.3,0.3), fun = fun_OLS_1,
                                    # ineqfun = fun_ineq, 
                                    # ineqLB = c(-0.999), 
                                    # ineqUB = c(0.999),
                                    LB = c(-Inf,-0.999,-0.999), 
                                    UB = c(Inf,0.999,0.999)))
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
    for(i in 1:dim(Y_mat)[2]){
      
      y_i<-Y_mat[(T/2+1):(T),i]
      
      x_i<-X_mat[(T/2+1):(T),i]
      v_i<-A%*%y_i - x_i*beta
      
      value<- value + t(v_i)%*%(Sigma_ - C)%*%v_i
    }
    
    return(value)   
    
  }
 
  estml<-try(suppressWarnings(solnp(c(3,0.3,0.3), fun = fun_OLS_2,
                                    # ineqfun = fun_ineq, 
                                    # ineqLB = c(-0.999), 
                                    # ineqUB = c(0.999),
                                    LB = c(-Inf,-0.999,-0.999), 
                                    UB = c(Inf,0.999,0.999)))
             ,silent = TRUE)
  
  if('try-error' %in% class(estml)){next} 
  
  est_gamma_2<-estml$pars
     
  ###########################
  est_gamma_unbias_J<- 2*est_gamma - (est_gamma_1+est_gamma_2)/2













est_A<-fun_A(est_gamma_unbias[2])

est_Sigma<-fun_Sigma(est_gamma_unbias[3])

est_Sigma_<-solve(est_Sigma)

c_Sigma<-c(t(l)%*%est_Sigma_%*%l)

est_V<-(diag(N_used)%x%est_A)%*%Y_used - X_used*est_gamma_unbias[1]

est_mu<-(diag(N_used)%x%(1/c_Sigma*t(l)%*%est_Sigma_))%*%est_V

############################################
###########################################

est_B<-fun_B(est_gamma_unbias[3])

est_U <- (diag(N_used)%x%solve(est_B)) %*% (est_V - est_mu%x%l)


##########################
est_omega<-c()
est_mu4<-c()
est_mu3<-c()
for(i in 1:N_used){
est_omega<-c(est_omega,mean((est_U[((i-1)*T_used+1):(i*T_used)])^2))
est_mu4<-c(est_mu4,mean((est_U[((i-1)*T_used+1):(i*T_used)])^4))
est_mu3<-c(est_mu3,mean((est_U[((i-1)*T_used+1):(i*T_used)])^3))

}

ml_GARCH<-function(par){
    
    tau<-par[1]
    nu<-par[2]
    
    
    value<-0
    for(i in 1:N_used){
      h0<- est_omega[i]*(1-nu-tau) + nu*est_omega[i]
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

estml_GARCH_0<-try(suppressWarnings(solnp(c(tau_all[1],nu_all[1]), fun = ml_GARCH,
                                            ineqfun = fun_ineq,
                                            ineqLB = c(0.001),
                                            ineqUB = c(0.999),
                                            LB = c(0.001,0.001),
                                            UB = c(0.999,0.999),
                                            control = list(trace=0)))
                     ,silent = TRUE)
  if('try-error' %in% class(estml_GARCH_0)){return(NULL)}
  
  est_zeta_0<-estml_GARCH_0$pars
  
  
  
  
  
  ##########################
  
  est_omega_1<-c()
  for(i in 1:N_used){
    est_omega_1<-c(est_omega_1,mean((est_U[((i-1)*T_used+1):((i-1)*T_used+T_used/2)])^2))
    
  }
  
  ml_GARCH_1<-function(par){
    
    tau<-par[1]
    nu<-par[2]
    
    value<-0
    for(i in 1:N_used){
      h0<- est_omega_1[i]*(1-nu-tau) + nu*est_omega_1[i]
      value<-value + log(dnorm(est_U[(i-1)*T_used+1],sd=sqrt(h0)))
      
      for(t in 2:(T_used/2)){

        h0<- est_omega_1[i]*(1-nu-tau) + nu*h0 + tau*est_U[(i-1)*T_used+t-1]^2


        value<-value+log(dnorm(est_U[(i-1)*T_used+t],sd=sqrt(h0)))
      }
      
      
    }
    
    return(-value)
  }
  
  
  estml_GARCH<-try(suppressWarnings(solnp(c(tau_all[1],nu_all[1]), fun = ml_GARCH_1,
                                          ineqfun = fun_ineq,
                                          ineqLB = c(0.001),
                                          ineqUB = c(0.999),
                                          LB = c(0.001,0.001),
                                          UB = c(0.999,0.999),
                                          control = list(trace=0)))
                   ,silent = TRUE)
  if('try-error' %in% class(estml_GARCH)){return(NULL)}
  
  est_zeta_1<-estml_GARCH$pars
  
  est_omega_2<-c()
  for(i in 1:N_used){
    est_omega_2<-c(est_omega_2,mean((est_U[((i-1)*T_used+T_used/2+1):((i-1)*T_used+T_used)])^2))

  }
  
  ml_GARCH_2<-function(par){
    
    tau<-par[1]
    nu<-par[2]
    
    
    value<-0
    for(i in 1:N_used){
      h0<- est_omega_2[i]*(1-nu-tau) + nu*est_omega_2[i]
      value<-value + log(dnorm(est_U[(i-1)*T_used+T_used/2+1],sd=sqrt(h0)))
      
      for(t in (T_used/2 + 2) : T_used){

        h0<- est_omega_2[i]*(1-nu-tau) + nu*h0 + tau*est_U[(i-1)*T_used+t-1]^2


        value<-value+log(dnorm(est_U[(i-1)*T_used+t],sd=sqrt(h0)))
      }
      
      
    }
    
    return(-value)
  }
  

  estml_GARCH<-try(suppressWarnings(solnp(c(tau_all[1],nu_all[1]), fun = ml_GARCH_2,
                                          ineqfun = fun_ineq,
                                          ineqLB = c(0.001),
                                          ineqUB = c(0.999),
                                          LB = c(0.001,0.001),
                                          UB = c(0.999,0.999),
                                          control = list(trace=0)))
                   ,silent = TRUE)
  if('try-error' %in% class(estml_GARCH)){return(NULL)}
  
  est_zeta_2 <- estml_GARCH$pars
  
 
  eps_est<-c()
  
  for(r in 1:N_used){
    
    u_est1_<-est_U[((r-1)*T_used+1):(r*T_used)]
    
    
    eps_est<-c(eps_est,diag(1/sqrt(fun_h(c(est_omega[r],est_zeta_1),u_est1_)))%*%u_est1_)
    
  }
  
  mu4<-mean(eps_est^4)
  
  est_zeta_unbias <- 2*est_zeta_0 - (est_zeta_1+est_zeta_2)/2
  
  if(sum(est_zeta_unbias)>1){
    return(NULL)
  }
   
  Gamma_unbias_<--solve(sec_deri_function(est_gamma_unbias))
  Lambda_unbias<-se_lambda(c(est_gamma_unbias,est_zeta_unbias))

  se_gamma_unbias_ <- sqrt(diag(Gamma_unbias_%*%Lambda_unbias%*%Gamma_unbias_))
 
  Gamma_zeta_unbias_<--solve(sec_deri_function_zeta(est_zeta_unbias))
  
  Omega_zeta_unbias<-Omega_function(est_zeta_unbias)
  
  Pi_zeta_unbias <- Pi_function(est_zeta_unbias)
  
  Lambda_zeta_unbias<- Lambda_function(c(est_gamma_unbias,est_zeta_unbias))
  
  
  Xi_unbias<-(mu4-1)*Omega_zeta_unbias + Pi_zeta_unbias%*%Gamma_unbias_%*%Lambda_unbias%*%Gamma_unbias_%*%t(Pi_zeta_unbias) + Pi_zeta_unbias%*%Gamma_unbias_%*%Lambda_zeta_unbias + t(Pi_zeta_unbias%*%Gamma_unbias_%*%Lambda_zeta_unbias)
 
  se_zeta_unbias_ <- sqrt(diag(Gamma_zeta_unbias_%*%Xi_unbias%*%Gamma_zeta_unbias_))
  

  Gamma_unbias_J<--solve(sec_deri_function(est_gamma_unbias_J))
  Lambda_unbias_J<-se_lambda(c(est_gamma_unbias_J,est_zeta_unbias)) 
 

 se_gamma_unbias_J <- sqrt(diag(Gamma_unbias_J%*%Lambda_unbias_J%*%Gamma_unbias_J))



return(list(est_gamma=est_gamma, se_gamma_=se_gamma_,
     est_zeta=est_zeta, se_zeta_=se_zeta_,
     est_gamma_unbias=est_gamma_unbias, se_gamma_unbias=se_gamma_unbias_,
     est_gamma_unbias_J=est_gamma_unbias_J, se_gamma_unbias_J=se_gamma_unbias_J,
     est_zeta_unbias=est_zeta_unbias, se_zeta_unbias_=se_zeta_unbias_))

}


NMI_all <- c(NMI_all, NMI)



if (dim(estim_left)[1] %% 2 == 0) {
  cat("Reached", dim(estim_left)[1], "iterations\n")
  flush.console()
}

estimation_left <- get_estimation_res(Y_left, X_left)
estimation_right <- get_estimation_res(Y_right, X_right)

if (is.null(estimation_left) | is.null(estimation_right)) {
  next
}


estim_left <- rbind(estim_left, c(estimation_left$est_gamma, estimation_left$est_zeta))
estim_right <- rbind(estim_right, c(estimation_right$est_gamma, estimation_right$est_zeta))

se_left <- rbind(se_left, c(estimation_left$se_gamma_, estimation_left$se_zeta_))
se_right <- rbind(se_right, c(estimation_right$se_gamma_, estimation_right$se_zeta_))

estim_left_unbiased <- rbind(estim_left_unbiased, c(estimation_left$est_gamma_unbias,estimation_left$est_gamma_unbias_J, estimation_left$est_zeta_unbias))
estim_right_unbiased <- rbind(estim_right_unbiased, c(estimation_right$est_gamma_unbias, estimation_right$est_gamma_unbias_J, estimation_right$est_zeta_unbias))

se_left_unbiased <- rbind(se_left_unbiased, c(estimation_left$se_gamma_unbias, estimation_left$se_gamma_unbias_J, estimation_left$se_zeta_unbias))
se_right_unbiased <- rbind(se_right_unbiased, c(estimation_right$se_gamma_unbias, estimation_right$se_gamma_unbias_J, estimation_right$se_zeta_unbias))


 


}




apply(estim_left, 2, mean) - parm_left
apply(estim_left_unbiased, 2, mean) - parm_left


apply(estim_right, 2, mean) - parm_right
apply(estim_right_unbiased, 2, mean) - parm_right

estim_left_unbiased




