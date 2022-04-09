##  Approx. Gibbs sampler for HBMs for grouped count data
library('invgamma')

## functions to sample from conditional posterior distributions

cal_sum_each_group = function(y_arr){
  
  # the quantity to be calculated for each group
  # a vector with length = n_data
  
  # sum of log_lambda of each group
  y_arr_each_group = rep(0, n_group)
  y_arr_each_group[1] = sum(y_arr[1:len_each_group[1]])
  
  for (i in 2:n_group){
    index_temp = seq(sum(len_each_group[1:(i-1)])+1, sum(len_each_group[1:i]), by=1)
    y_arr_each_group[i] = sum(y_arr[index_temp])
  } 
  
  return(y_arr_each_group)
}


# hyper parameters mu and sigma2. Conditional posterior for hyper parameters. Conjugate.

draw_mu_x = function(x, sigma2_x) {
  # input: x - an array, e.g. vector of alpha_j; sigma2_x: a positive 'double'
  # output: a 'double'
  mu_mu_x = tau2_all*sum(x)/(n_group*tau2_all+sigma2_x)
  sigma2_mu_x = tau2_all*sigma2_x/(n_group*tau2_all+sigma2_x)

  return(rnorm(1, mu_mu_x, sqrt(sigma2_mu_x)))
}

draw_sigma2_x = function(x, mu_x) {
  b_sigma2_x = b_all + sum((x-mu_x)^2)/2
  # note that invgamma pdf is defined over the support x > 0
  return(rinvgamma(1, a_sigma2_x, b_sigma2_x))
}

# proposal function
prop_norm = function(coeff_curr, coeff_curr_id, hyper_param_curr, ...){
  # input: coeff_curr - a vector; 
  # e.g.: beta1, beta2
  # coeff_curr_id - an integer.
  # hyper_param_curr is redundant here, but is needed in the proposal using approximate likelihood
  # output: a vector
  sd = c(0.02/34, 3e-6/13)
  
  coeff_curr_1d = coeff_curr[coeff_curr_id, ]
  
  coeff_curr_new_1d = rnorm(n_group, mean=coeff_curr_1d, sd=rep(sd[coeff_curr_id], n_group))
  
  return(coeff_curr_new_1d)
}


cal_exact_mean = function(y){
  if(y<1){
    stop('Error: the input must be >=1 in calculating the exact mean of the approximate distribution')
  }
  euler_con = 0.57721566
  if (y==1) {
    exact_mean = -euler_con
  }else{
    k = seq(1,y-1)
    harm_num = sum(1/k)    # harmonic number
    exact_mean = -euler_con + harm_num
  }
  return(exact_mean)
}

cal_exact_var = function(y){
  if(y<1){
    stop('Error: the input must be >=1 in calculating the exact mean of the approximate distribution')
  }
  if(y==1){
    exact_var = pi^2/6
  }else{
    seq_neg_part = seq(1,y-1)
    exact_var = (pi^2/6-sum(1/seq_neg_part^2))^2
  }
  return(exact_var)
}


prop_approx_likeli = function(coeff_curr, coeff_curr_id, hyper_param_curr, ...){
  # coeff_curr: a 2d array, dim = n_coeff, n_group
  # hyper param curr: a vector of mu and sigma2

  coeff_curr_rep = apply(coeff_curr, 1, function(x) rep(x, len_each_group))

  prod_part_temp = exact_mean - (rowSums(coeff_curr_rep*X) - coeff_curr_rep[, coeff_curr_id]*X[, coeff_curr_id])
  numerator = hyper_param_curr[1] + hyper_param_curr[2]*
              cal_sum_each_group(1/exact_var*X[, coeff_curr_id]*prod_part_temp)
  denominator = hyper_param_curr[2]*cal_sum_each_group(1/exact_var*X[,coeff_curr_id]^2) + 1
  
  mu_hat = numerator/denominator
  sigma_hat = sqrt(hyper_param_curr[2]/denominator)
  
  coeff_curr_new_1d = rnorm(n_group, mean=mu_hat, sd=sigma_hat)
  
  # cat('\nmu', mu_hat)
  # cat('\nsigma', sigma_hat)
  # cat('\ncoeff_new_1d', coeff_curr_new_1d)
  
  return(coeff_curr_new_1d)
}


log_post_orig = function(coeff_curr, coeff_curr_id, hyper_param_curr) {
  # prior: p(beta_j|mu_beta, sigma_beta)
  # output: a vector of conditional posterior of coeff, dim =c(n_group,1).
  
  # log of prior
  coeff_curr_subset = coeff_curr[coeff_curr_id,]
  log_prior = dnorm(coeff_curr_subset, hyper_param_curr[1], hyper_param_curr[2], log=T)
  
  # log of likeli: the prod over i of p(y_ij|log_lambda_ij). log_lambda_ij = alpha_j + beta_j*x1_ij + gamma_j*x2_j
  
  coeff_curr_rep = apply(coeff_curr, 1, function(x) rep(x, len_each_group))
  log_lambda = rowSums(coeff_curr_rep*X)
  
  # the constant term, - log_fac_approx(y), is removed. Cancelled out in log_post_new - log_post_current
  log_likeli = -exp(log_lambda) + y*log_lambda
  log_likeli_each_group = cal_sum_each_group(log_likeli)  # sum of log_lambda of each group
  
  # log of posteriror
  log_post = log_prior + log_likeli_each_group
  
  return(log_post)
}


# metropolis draw
mh_draw = function(coeff_curr, coeff_curr_id, hyper_param_curr, targ_dens_mh=log_post_orig, prop_mh) {
  # input: coeff_curr - a 2d array of samples of coeffs;
  # coeff_curr_id - id of group-level or indiv-level coeff to be updated
  # hyper_param_curr = hyper_param_curr[i_coeff, ]. A vector. e.g. (mu, sigma2)
  # output: 
  # generate new x_ii, evaluate with x, return new x_ii or reject it
  coeff_curr = as.matrix(coeff_curr)  # keep the dimension
  coeff_curr_new = coeff_curr
  coeff_curr_new[coeff_curr_id, ] = prop_mh(coeff_curr, coeff_curr_id, hyper_param_curr) # draw a new sample from current sample
  
  # new coeff values, but evaluating the dens requires the current values of other coeffs and hyper params.
  log_rate = targ_dens_mh(coeff_curr_new, coeff_curr_id, hyper_param_curr) - 
             targ_dens_mh(coeff_curr, coeff_curr_id, hyper_param_curr)
  is_accept = log(runif(n_group)) < log_rate
  
  # print(coeff_curr_new)
  
  coeff_curr[coeff_curr_id, is_accept] = coeff_curr_new[coeff_curr_id, is_accept]
  
  accept_count[coeff_curr_id, is_accept] <<- accept_count[coeff_curr_id, is_accept] + 1
  
  return(coeff_curr)
}

gibbs_with_mh = function(coeff_curr, hyper_param_curr, targ_dens, prop_fun){
  
  # update hyperpriors
  for (i_coeff in 1:n_coeff) {
    hyper_param_curr[i_coeff, 1] = draw_mu_x(coeff_curr[i_coeff, ], hyper_param_curr[i_coeff, 2])
    hyper_param_curr[i_coeff, 2] = draw_sigma2_x(coeff_curr[i_coeff, ], hyper_param_curr[i_coeff, 1])
  }
  
  
  for (i_coeff in 1:n_coeff){
    coeff_curr = mh_draw(coeff_curr, i_coeff, hyper_param_curr[i_coeff, ],
                         targ_dens_mh = targ_dens,  prop_mh = prop_fun)
  }
  
  return(list(coeff_curr, hyper_param_curr))
}


gibbs_without_mh = function(coeff_curr, hyper_param_curr, approx_fun = prop_approx_likeli){

  # update hyperpriors
  for (i_coeff in 1:n_coeff) {
    hyper_param_curr[i_coeff, 1] = draw_mu_x(coeff_curr[i_coeff, ], hyper_param_curr[i_coeff, 2])
    hyper_param_curr[i_coeff, 2] = draw_sigma2_x(coeff_curr[i_coeff, ], hyper_param_curr[i_coeff, 1])
  }

  for (i_coeff in 1:n_coeff) {
    # cat('\ni coeff', i_coeff,'\n')
    # print(coeff_curr[i_coeff, ])
    coeff_curr[i_coeff, ] = prop_approx_likeli(coeff_curr, i_coeff, hyper_param_curr[i_coeff,])
  }
  
  return(list(coeff_curr, hyper_param_curr))
}


## function to run the mh-in-gibbs

fit_gibbs = function(data_gibbs, n_keep, n_warmup, n_chain){
  
  y <<- data_gibbs$y
  X <<- data_gibbs$X 
  n_data <<- data_gibbs$n_data 
  n_group <<- data_gibbs$n_group 
  n_coeff <<- data_gibbs$n_coeff
  group_id <<- data_gibbs$group_id
  len_each_group <<- data_gibbs$len_each_group

  
  # values of hyper params 
  tau2_all <<- 1
  a_all <<- 2
  b_all <<- 2
  a_sigma2_x <<- a_all + n_group
  
  n_total_sample = n_warmup + n_keep
  n_hyper_param = 2
  
  
  meth_name_all = c('MH-in-Gibbs, exact conditional posterior','MH-in-Gibbs, approximate conditional posterior',
                'Gibbs, approximate conditional posterior')
  meth_name = meth_name_all[c(3)]
  
  dens_fun_list = c(prop_approx_likeli)    #c(prop_norm, prop_approx_likeli)
  n_meth <<- length(meth_name)
  
  hyper_param = array(NA, dim = c(n_coeff, n_total_sample, n_hyper_param, n_chain, n_meth))

  coeff = array(NA, dim = c(n_coeff, n_total_sample, n_group, n_chain, n_meth))
  
  accept_count_all = array(0, dim = c(n_coeff, n_group, n_chain, n_meth))
  
  t_run = array(0, dim = c(n_meth))
  
  
  for(i_meth in 1:n_meth){
    
    t_start  = Sys.time()
    
    for (i_chain in 1:n_chain){
      cat('\n\nMethod:', meth_name[i_meth])
      cat('\nChain: ', i_chain, '\n')
      
      # initialization
      hyper_param[ , 1, , i_chain, i_meth] = array(abs(rnorm(1, 1, 1/2)), dim = c(n_coeff, n_hyper_param))
      coeff[, 1, , i_chain, i_meth] = array(rnorm(1, 0.001, 0.001/2), dim = c(n_coeff, n_group))
      if (grepl('MH', meth_name[i_meth])==T){
        accept_count <<- array(0, dim = c(n_coeff, n_group))
      }
      
      for (i_samp in 2:n_total_sample){
        if (i_samp%%1000==0){
        cat('\nsample: ', i_samp)
        }
        
        if (grepl('MH', meth_name[i_meth])==T){
          coeff_hyper_param_return = gibbs_with_mh(
              coeff[, i_samp-1, , i_chain, i_meth], hyper_param[ ,i_samp-1, ,i_chain, i_meth],
              log_post_orig, dens_fun_list[[i_meth]])

        }else{
          coeff_hyper_param_return = gibbs_without_mh(
              coeff[, i_samp-1, , i_chain, i_meth], hyper_param[ ,i_samp-1, ,i_chain, i_meth])
        }
        
        coeff[, i_samp, , i_chain, i_meth] = as.matrix(coeff_hyper_param_return[[1]])
        hyper_param[, i_samp, , i_chain, i_meth] = as.matrix(coeff_hyper_param_return[[2]])
        
      }
      
      if (grepl('MH', meth_name[i_meth])==T){
        accept_count_all[,,i_chain, i_meth] = accept_count
      }
      
    }
    
    t_end = Sys.time()
    t_run[i_meth] = as.numeric(as.POSIXct(t_end)-as.POSIXct(t_start ), units="secs")
  }
  cat('\n\n Total sampling time for approx. Gibbs:', t_run, '\n')
  
  return(list(coeff, hyper_param, accept_count_all, t_run))
}



