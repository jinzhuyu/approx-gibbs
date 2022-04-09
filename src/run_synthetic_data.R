# library(dplyr)  # convert factor to numeric
library('tidyr')
library('rstan')
setwd('C:/GitHub/approx_gibbs_for_HBM')

# dependent source code
source('data_preprocessing.R')
source('gibbs_all_vectorized.R')
source('diagnostics_and_plotting.R')
options(digits = 5)


## generate simulated data
generate_data = function(n_point_group, n_group, n_cov, normalize=1, add_noise=1){
  
  n_total = n_point_group*n_group
  
  # generate covariates
  x_group = runif(n_group, 1e3, 1e5)  
  x_group = matrix(rep(x_group, each = n_point_group), n_total,1)
  
  x_indiv = matrix(NA, n_total, n_cov-1)
  lb = c(0.1,0.1,0.1,1,0.5,10)
  ub = c(2,1,0.5,10,5,100)
  for (i in 1:ncol(x_indiv)){
    for (j in 1:n_group){
      row_index = seq(((j-1)*n_point_group+1), j*n_point_group, by=1)
      x_indiv[row_index,i] = sort(runif(n_point_group,lb[i],ub[i]))
    }
  }
  
  # generate coeff
  mu = 1e-3
  sd = mu
  alpha = abs(rnorm(n_group, mu, sd))
  alpha_rep = rep(alpha, each = n_point_group)
  coeff_indiv = matrix(abs(rnorm((n_cov-1)*n_group, mu, sd)), n_group,n_cov-1)
  coeff_indiv_rep = matrix(NA, n_total, n_cov-1)
  coeff_indiv = matrix(abs(rnorm((n_cov-1)*n_group, mu, sd)), n_group,n_cov-1)
  coeff_indiv_rep = matrix(NA, n_total, n_cov-1)
  for (j in 1:n_group){
    row_index = seq(((j-1)*n_point_group+1),j*n_point_group,by=1)
    coeff_indiv_rep[row_index,] = rep_row(coeff_indiv[j,], n_point_group)
  }
  
  # calculate y
  rate = exp(alpha_rep+rowSums(coeff_indiv_rep*x_indiv))
  rate_norm = rate
  for (j in 1:n_group){
    row_index = seq(((j-1)*n_point_group+1),j*n_point_group,by=1)
    lb_norm = max(abs(rnorm(1,5e-3,5e-3)), 1e-3)
    rate_norm[row_index] = min_max(rate[row_index],lb_norm)
  }
  
  y = round(rate_norm*x_group, 0)
  
  x = as.matrix(cbind(x_indiv, x_group))
  x[,n_cov] = log(x[,n_cov])
  
  # normalize x
  if (normalize==1){
    for (i in 1:n_cov){
      x[,i] = min_max(x[,i],lb_norm=1e-4)
    }
  }
  
  # add noise
  if (add_noise==1){
    cv = 0.50
    mu_y_noise = y/100
    sigma_y_noise = mu_y_noise*cv
    len_y = length(y)
    y = y + (-1)^sample(1:len_y, replace=T)*round(rnorm(len_y, mu_y_noise, sigma_y_noise), 0)
  }
  
  # handle values of x and y that are too small
  y[y<5] = 5
  x[x==0] = 1e-5
  
  data_df = as.data.frame(cbind(x, y))
  
  return(data_df)
}


n_cov = 6
n_group = 10
n_point_group = 20
group_attr_id = n_cov

syn_data = generate_data(n_point_group, n_group, n_cov, normalize=1, add_noise=0)
data_all = prepare_data(syn_data, group_attr_id)

# # sampling parameters
n_keep = 10000L
n_warmup = 5000L
n_chain = 4L

fit_model(data_all, group_attr_id)

