library(tidyverse)
library(dplyr)
load("cycle_parameter.RData")
#model
estimate_ts_model <- function(Y_obs, S_matrix, cycle, D, K, freq_inital, N){
  
  calculate_likelihood <- function(params){
    
    ## 1. Parameter disassembly
    N_i <- params[1:K]
    alpha <- params[K+1]
    
    ## 2. relative abundance
    r_i <- N_i / sum(N_i)
    
    ## 3. inital
    p_ij <- matrix(0, nrow = K, ncol = K)
    q_i <- numeric(K)
    q_i_other <- numeric(K)
    X_i <- numeric(K)
    Y_ij <- matrix(0, nrow = K, ncol = K)
    
    ## 4. p_ij = alpha * S_ij * r_j
    for(i in seq_len(K)){
      for(j in seq_len(K)){
        p_ij[i, j] <- alpha * S_matrix[i, j] * r_i[j]
      }
      q_i[i] <- sum(p_ij[i, ])
      q_i_other[i] <- q_i[i] - p_ij[i, i]  
    }
    
    q_i_other[q_i_other > 0.99] <- 0.99
    q_i_other[q_i_other < 1e-8] <- 1e-8
    
    ## 5. 
    ##    X_i = N_i * (2 - q_i_other)^c
    for(i in seq_len(K)){
      X_i[i] <- N_i[i] * (2 - q_i_other[i])^cycle
    }
    
    ## 6. Y_ij (i != j)
    ##    Y_ij = (p_ij * N_i / q_i_other) * [2^c - (2 - q_i_other)^c]
    for(i in seq_len(K)){
      for(j in seq_len(K)){
       
          term1 <- (p_ij[i, j] * N_i[i] / q_i_other[i]) * (2^cycle)
          term2 <- (p_ij[i, j] * N_i[i] / q_i_other[i]) * ((2 - q_i_other[i])^cycle)
          Y_ij[i, j] <- term1 - term2
      
      }
    }
    
    ## 7.  M_total = sum(X_i) + sum_{i!=j} Y_ij
    M_total <- sum(X_i) + sum(Y_ij)
    
    
    if(M_total <= 0) M_total <- 1e-8
    
    ## 8. log-likelihood
    log_lik <- 0
    
    ## 8.1 
    for(i in seq_len(K)){
      for(j in seq_len(K)){
        
          expected <- D * Y_ij[i, j] / M_total
          if(expected > 0 && Y_obs[i, j] > 0){
            log_lik <- log_lik + Y_obs[i, j] * log(expected) - expected
          } else if(expected > 0){
            log_lik <- log_lik - expected  
          }
      }
    }
    
    
    if(!is.finite(log_lik)) log_lik <- -1e10
    
    return(log_lik)
  }
  
  ## 9. initial value
  N_initial <- N * freq_inital
  
  N_initial[N_initial < 1e-6] <- 1e-6
  alpha_initial <- 1e-4
  
  initial_params <- c(N_initial, alpha_initial)
  
  ## 10. L-BFGS-B 
  result <- optim(
    par = initial_params,
    fn = calculate_likelihood,
    method = "L-BFGS-B",
    lower = c(rep(1e-6, K), 1e-8),      
    upper = c(rep(1e8, K), 0.5),         
    control = list(fnscale = -1, maxit = 1000)
  )
  
  ## 11. output
  N_estimated <- result$par[1:K]
  alpha_estimated <- result$par[K+1]
  r_estimated <- N_estimated / sum(N_estimated)
  
  
  final_lik <- calculate_likelihood(result$par)
  
  return(list(
    N_estimated = N_estimated,
    alpha_estimated = alpha_estimated,
    r_estimated = r_estimated,
    log_likelihood = final_lik,
    convergence = result$convergence,
    message = result$message
  ))
}

result <- estimate_ts_model(Y_obs, S_matrix,effective_cycle, D, K,freq_inital,N)


# print result
print(paste("Alpha:", round(result$alpha_estimated, 6)))
print("r_i:")
print(round(result$r_estimated, 4))
print(paste("log_likehood:", round(result$log_likelihood, 2)))
