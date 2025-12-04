 
apply(est_result,2, var) 

df <- data.frame(
  index = seq_along(est_result[,5]),  # 生成 1,2,3,... 的序号
  value = est_result[,5]
)

# 按数值排序
df_sorted <- df[order(df$value), ]    
plot(df_sorted[,2])


rownames(df_sorted) <- NULL
sorted_values <- df_sorted$value


begin_index <- 1
end_index <- (nrow(df_sorted))



delta_thresh_list <- seq(0.01, 0.5, by = 0.1)

# delta_thresh_list <- seq(1, 20, by = 1)

# delta_thresh_list <- c(100)
final_used <- list()
BIC_min <- 1e15

delta_thresh = 0.1
for (delta_thresh in delta_thresh_list) {
  print(delta_thresh)
  
  fit <- detect_change_points(
    sorted_values = df_sorted$value,
    orig_index    = df_sorted$index,
    delta_thresh  = delta_thresh,
    min_size      = 50
  )
  
  BIC_value <- 0
  for (i in 1:length(fit$leaf_groups)) {
    Y_i <- ymat_used[fit$leaf_groups[[i]],, drop = FALSE]
    # x1_i <- x1[fit$leaf_groups[[i]],, drop = FALSE]
    # x4_i <- x4[fit$leaf_groups[[i]],, drop = FALSE]
    # 
    # print(get_BIC_res(Y_i))
    
    BIC_value <- BIC_value + get_BIC_res(Y_i)
    
  }
  BIC_ <- BIC_value + log(ncol(ymat_used) * nrow(ymat_used)) * length(fit$leaf_groups) * 3
  
  if (BIC_ < BIC_min) {
    BIC_min <- BIC_
    final_used <- list(fit = fit, delta_thresh = delta_thresh, BIC = BIC_)
  }
}

final_used <- list(fit = fit, delta_thresh = delta_thresh, BIC = BIC_)

final_used$fit

final_used$delta_thresh

