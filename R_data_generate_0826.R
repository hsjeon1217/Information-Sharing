########################################################################################

sub.data.generate.0826 = function(m, m1, n_i, v.delta, intercept, v.sigma, sigma_ratio){
  
  # n is the sample size for each treatment level :
  v.sigma = v.sigma*sigma_ratio;
  
  Y.ctrl  <- (matrix( rnorm((n_i*m), 0, 1), nrow = m, ncol = n_i ) * v.sigma) + intercept
  Y.trt   <- (matrix( rnorm((n_i*m), 0, 1), nrow = m, ncol = n_i ) * v.sigma) + intercept
  
  # For the first m1 genes, add fixed effect :
  # if g is in H_A (g in 1:m1)
  Y.trt[1:m1,] = Y.trt[1:m1,] + v.delta;
  
  # m x n dimension :
  # The order of Y and X is trt -> ctrl.
  Y = cbind(Y.trt, Y.ctrl)
  X = model.matrix(~ 0 + as.factor(c(rep(1,n_i), rep(2,n_i))))
  colnames(X) = c("trt","ctrl")
  
  return(list(Y = Y, X = X))
  
}

sub.data.generate.0521 = function(m, m1, n_i, delta, intercept, v.sigma, sigma_ratio){
  
  # n is the sample size for each treatment level :
  v.sigma = v.sigma*sigma_ratio;
  
  Y.ctrl  <- (matrix( rnorm((n_i*m), 0, 1), nrow = m, ncol = n_i ) * v.sigma) + intercept
  Y.trt   <- (matrix( rnorm((n_i*m), 0, 1), nrow = m, ncol = n_i ) * v.sigma) + intercept
  
  # For the first m1 genes, add fixed effect :
  # if g is in H_A (g in 1:m1)
  Y.trt[1:m1,] = Y.trt[1:m1,] + delta;
  
  # m x n dimension :
  # The order of Y and X is trt -> ctrl.
  Y = cbind(Y.trt, Y.ctrl)
  X = model.matrix(~ 0 + as.factor(c(rep(1,n_i), rep(2,n_i))))
  colnames(X) = c("trt","ctrl")
  
  return(list(Y = Y, X = X))
  
}

sub.data.generate.0521 = function(m, m1, n_i, delta, intercept, v.sigma, sigma_ratio){
  
  # n is the sample size for each treatment level :
  v.sigma = v.sigma*sigma_ratio;
  
  Y.ctrl  <- (matrix( rnorm((n_i*m), 0, 1), nrow = m, ncol = n_i ) * v.sigma) + intercept
  Y.trt   <- (matrix( rnorm((n_i*m), 0, 1), nrow = m, ncol = n_i ) * v.sigma) + intercept
  
  # For the first m1 genes, add fixed effect :
  # if g is in H_A (g in 1:m1)
  Y.trt[1:m1,] = Y.trt[1:m1,] + delta;
  
  # m x n dimension :
  # The order of Y and X is trt -> ctrl.
  Y = cbind(Y.trt, Y.ctrl)
  X = model.matrix(~ 0 + as.factor(c(rep(1,n_i), rep(2,n_i))))
  colnames(X) = c("trt","ctrl")
  
  return(list(Y = Y, X = X))
  
}

data.generate.0521 = function(m, m1, n_pilot, n_main, delta, sigma_ratio){
  
  # Generate sigmas from inverse chisq-distribution
  
  v.sigma = rinvchisq(n = m, df = 5, ncp = 0) # mean = 1/3  
  
  # The intercept can be interpreted as a trial effect : 
  
  data_pilot = sub.data.generate.0521(m, m1, n_i = n_pilot, delta, intercept = 0, v.sigma, sigma_ratio = sigma_ratio)
  data_main = sub.data.generate.0521(m, m1, n_i = n_main, delta, intercept = 1, v.sigma, sigma_ratio = 1)
  
  Y_pilot = data_pilot$Y; X_pilot = data_pilot$X;
  Y_main = data_main$Y; X_main = data_main$X;
  
  data_list = list( Y_pilot = Y_pilot, 
                    X_pilot = X_pilot,
                    Y_main = Y_main, 
                    X_main = X_main )
  
  return(data_list)
  
}

MSE.generate.0521 = function(n_i, X, Y){
  
  # This degree of freedom is only for this case :
  df = 2*n_i - 2
  b = solve(t(X) %*% X) %*% t(X) %*% t(Y)
  Yhat = t(X %*% b)
  SSE = apply((Y - Yhat)^2,1,sum)
  MSE = SSE/df
  
  mu = t(b)
  colnames(mu) = colnames(X)
  
  test_source = list(MSE = MSE, df = df, mu = mu)
  
  return(test_source)
  
}

p_values.generate.0521 = function( Y_pilot, X_pilot, 
                                   Y_main, X_main, 
                                   n_pilot, n_main ){
  
  # test sources:
  TS_pilot = MSE.generate.0521(n_i = n_pilot, X = X_pilot, Y = Y_pilot)
  TS_main = MSE.generate.0521(n_i = n_main, X = X_main, Y = Y_main)
  
  MSE = data.frame( MSE_pilot = TS_pilot$MSE, 
                    MSE_main = TS_main$MSE );
  
  mu_hat = cbind(TS_pilot$mu, TS_main$mu) %>% data.frame(.) %>% 
    setNames(c( paste0("mu_pilot_",colnames(TS_pilot$mu)),
                paste0("mu_main_",colnames(TS_main$mu)) ))
  
  # p-value sources :
  df_sources = cbind(MSE, mu_hat) %>% data.frame(.);
  
  # p-values generation :
  # Two-Tailed Tests :
  # 2 * pt( abs(tx ̄), df = n-1, lower.tail=FALSE)
  # For T_comb, CS - degree of freedom is used :
  
  p_df = df_sources %>% mutate( T_pilot = (mu_pilot_trt - mu_pilot_ctrl)/sqrt(2*MSE_pilot/n_pilot),
                                T_main = (mu_main_trt - mu_main_ctrl)/sqrt(2*MSE_main/n_main),
                                T_comb = ((mu_pilot_trt + mu_main_trt)/2 - (mu_pilot_ctrl + mu_main_ctrl)/2)/sqrt((MSE_pilot/n_pilot + MSE_main/n_main)/2),
                                df_comb = (MSE_pilot/n_pilot + MSE_main/n_main)^2/(MSE_pilot^2/(n_pilot^2*TS_pilot$df) + MSE_main^2/(n_main^2*TS_main$df)) 
  ) %>%
    mutate( p_pilot = 2*pt(abs(T_pilot), df = TS_pilot$df, lower.tail = F),
            p_main = 2*pt(abs(T_main), df = TS_main$df, lower.tail = F),
            p_comb = 2*pt(abs(T_comb), df = df_comb, lower.tail = F)) %>% 
    dplyr::select(p_pilot, p_main, p_comb)
  
  return(p_df)
  
}

########################################################################################