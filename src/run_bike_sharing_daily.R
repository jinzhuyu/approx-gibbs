################
# Daily bike sharing data of each day


################
# library(dplyr)  # convert factor to numeric
library('tidyr')
library('rstan')

# dependent source code
source('data_preprocessing.R')
source('gibbs_all_vectorized.R')
source('diagnostics_and_plotting.R')
options(digits = 5)

################
# import data
data_list = read.csv("../test_data/bike_sharing_data/bike_sharing_day.csv", header = T)

# select important features temp, hum, working day, windspeed, cnt
  # refs.: https://github.com/pgebert/bike-sharing-dataset
data_list_new = data_list[, c(8,10,12,13,14,16)]

# scale the unregistered counts to avoid crashing R when running stan
data_list_new[, 5] = data_list_new[, 5]/100

# add a group-level covariate~norm(0,1)
# find group id
n_data = length(data_list_new[,1])
group_id = data_list_new[,1] + 1  # working day is a binary factor covariate
n_group = max(group_id)
cov_group_rep = rep(0, n_data)
cov_group_rep[group_id==1] = abs(rnorm(1,0,0.5))
cov_group_rep[group_id==2] = abs(rnorm(1,0,0.5))

# add covariates to data
data_list_new = cbind(data_list_new, cov_group_rep)

# reorder the columns
data_list_select = data_list_new[,c(2,3,4,5,7,6)]

# data features
n_cov = length(data_list_select[1,])-1
group_attr_id = n_cov

# sort data based on group_id first.
len_each_group = find_len_each_group(data_list_select, group_attr_id, n_data, n_group)

# normalize individual level covariates
for (i in 1:(n_cov-1)){
  data_list_select[,i] = min_max(data_list_select[,i])
}

# prepare data for the two samplers
data_all = prepare_data(data_list_select, group_attr_id)

# # sampling parameters
n_keep = 10000L
n_warmup = 5000L
n_chain = 2L

# fit model
data_stan <<- data_all[[1]]
data_gibbs <<- data_all[[2]]

# stan
stan_fit <<- stan(file='1d_HBM_multi_attri.stan', data = data_stan, iter = n_keep+n_warmup, warmup=n_warmup, chains = n_chain)

# # gibbs
y <<- data_gibbs$y
exact_mean <<- sapply(y, cal_exact_mean)
exact_var <<- sapply(y, cal_exact_var)
gibbs_fit <<- fit_gibbs(data_gibbs, n_keep, n_warmup, n_chain)
index_good <<- (n_warmup+1):(n_warmup+n_keep)

# value of performance metrics
metric_all = cal_metric(stan_fit, gibbs_fit, y, print_out=1)

metric_all
