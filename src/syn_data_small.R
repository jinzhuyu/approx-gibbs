# library(dplyr)  # convert factor to numeric
library('tidyr')

setwd('C:/code/approx-gibbs/src')
library('rstan')
# import dependent functions
source('data_preprocess.R')
source('approx_gibbs.R')
source('diagnose_and_plot.R')
options(digits = 5)



pdf_exp = function(x, lambda){
  return( lambda*exp(-lambda*(max(x) - x + 1)) )
}

## generate simulated data
generate_data = function(n_point_group, n_group, n_cov, normalize=1){
  
  n_total = n_point_group*n_group
  
  # generate covariates
  x_group_max = 5
  x_group_min = 1
  # choose seed to get the desired proportion of small counts
    # the samples are somewhat close to uniformly distributed
  if (x_group_max == 10) {
    seed = 4  # <=5 rate: 0.496
  }
  if (x_group_max == 5) {
    seed = 4  # <=3 rate: 0.6625
  }
  if (x_group_max == 15) {
    seed = 7  # <=5 rate: 0.30312 
  }
  set.seed(seed)
  print(seed)
  x_group_values = seq(x_group_min, x_group_max, by=1)
  prob = pdf_exp(x_group_values, lambda=0.7)
  x_group = sample(x_group_values, n_group, replace=T, prob=prob)  #runif(n_group, 1, 5)  
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
  mu = 1e-1
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
    lb_norm = max(abs(rnorm(1,5e-2,5e-2)), 1e-1)
    rate_norm[row_index] = min_max(rate[row_index],lb_norm)
  }
  
  # y = round(rate_norm*x_group, 0)
  y = ceiling(rate_norm*x_group)
  
  x = as.matrix(cbind(x_indiv, x_group))
  x[,n_cov] = log(x[,n_cov])
  
  # normalize x
  if (normalize==1){
    for (i in 1:n_cov){
      x[,i] = min_max(x[,i], lb_norm=1e-4)
    }
  }
  
  # handle values of x and y that are too small
  
  print(length(y[y<=5])/length(y))
  png(filename=paste0("hist_count_max_count_",x_group_max, ".png", sep=''))
  hist(y)
  dev.off()
  
  x[x==0] = 1e-5
  
  data_df = as.data.frame(cbind(x, y))
  
  return(data_df)
}


n_cov = 3
n_group = 4
n_point_group = 60
group_attr_id = n_cov

syn_data = generate_data(n_point_group, n_group, n_cov, normalize=1)
# print(syn_data)
data_all = prepare_data(syn_data, group_attr_id)

# # sampling parameters
n_keep = 10000
n_warmup = 10000
n_chain = 4

fit_model(data_all, group_attr_id)





