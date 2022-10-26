################
# covid_19 test data with small counts (<=5) only 
################

setwd('C:/code/approx-gibbs/src')
library('rstan')
# import dependent functions
source('data_preprocess.R')
source('approx_gibbs.R')
source('diagnose_and_plot.R')
options(digits = 5)

## import and adapt data
# import data
data_list = read.csv("../test_data/covid19_test_data/covid19_data_raw.csv", header = T)

# select useful columns
data_list = data_list[, c(1,3,6,7)] 
t_exp_before = 5 + 1
data_list[,2] = data_list[,2] + t_exp_before

# keep only small positive count to see the performance
test_pos = data_list[, 4] 
data_list_new = data_list[which(test_pos >=1 & test_pos <= 5),]
# run to test if increasing counts by 5 will lead to better performance
# data_list_new[, 4] = data_list_new[, 4] + 5

# add a group-level covariate~norm(0,1)
# find group id
n_data = length(data_list_new[,1])
# n_data
group_id = rep(1,n_data)
group_id_temp = 1
cov_group_rep = rep(rnorm(1,0,0.25), n_data)
for(i in 2:n_data){
  study_name_curr = data_list_new[i,1]
  if(data_list_new[i,1]==data_list_new[i-1,1]){
    group_id[i] = group_id_temp
    cov_group_rep[i] = cov_group_rep[i-1]
    i = i+1
  }else{
    group_id_temp = group_id_temp+1
    group_id[i] = group_id_temp
    cov_group_rep[i] = rnorm(1,0,0.25)
  }
}

# generate covariates based on the number of days after exposure, according to the original model
log_t = log(data_list_new$day)
log_t2 = log_t^2
log_t3 = log_t^3

# add covariates to data
data_list_new = cbind(data_list_new,log_t, log_t2, log_t3, cov_group_rep)

# remove study name and reorder the columns
data_list_select = data_list_new[,c(5,6,7,3,8,4)]

# data features
n_cov = length(data_list_select[1,]) - 1
group_attr_id = n_cov
n_group = max(group_id)
len_each_group = find_len_each_group(data_list_select, group_attr_id, n_data, n_group)

# normalize individual level covariates
for (i in 1:(n_cov-1)){
  data_list_select[,i] = min_max(data_list_select[,i])
}

# format data for each model
data_all = prepare_data(data_list_select, group_attr_id)
data_stan = data_all[[1]]
data_gibbs = data_all[[2]]

# sampling parameters
n_keep = 10000
n_warmup = 8000
n_chain = 4

# fit model
# stan
# stan_fit = stan(file='1d_HBM_multi_attri.stan', data = data_stan,
#                 iter = n_keep+n_warmup,warmup=n_warmup, chains = n_chain)
# gibbs
y = data_gibbs$y
exact_mean = sapply(y, cal_exact_mean)
exact_var = sapply(y, cal_exact_var)

gibbs_fit = fit_gibbs(data_gibbs, n_keep, n_warmup, n_chain)
index_good = (n_warmup+1):(n_warmup+n_keep)

# value of performance metrics
metric_all = cal_metric(stan_fit, gibbs_fit, y, print_out=TRUE)
metric_all


