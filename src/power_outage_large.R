################
# power outage counts after weather events accross the US

################
# dependent source code and packages
setwd('C:/code/approx-gibbs/src')
library('rstan')
# import dependent functions
source('data_preprocess.R')
source('approx_gibbs.R')
source('diagnose_and_plot.R')
options(digits = 5)

################
## import and select data
data_list = read.csv("../test_data/power_outage_data_large/outage_data_all.csv", header = T)

# select data using event_id
event_id =  18
# data_list_select = data_list[data_list$event_id==event_id, 6:12]
data_list_select = data_list[data_list$event_id==event_id, 2:8]

# data covariates
n_cov = length(data_list_select[1,])-1
group_attr_id = n_cov
n_data = length(data_list_select[,1])

  
# scale individual-level covariate group by group
group_attr_id_temp = 6
group_id = find_group_id(data_list_select, group_attr_id_temp, n_data)

for (i in 1:max(group_id)){
  for (j in 1:(length(data_list_select[1,])-2)) {
    data_list_select[group_id==i, j] = min_max(data_list_select[group_id==i, j])
  }
}

# # scale group covariate
data_list_select[, n_cov] = min_max(data_list_select[, n_cov])

ub_approx = 3
data_list_select$outage_count[data_list_select$outage_count<ub_approx] = ub_approx


# prepare data for the two samplers
data_all = prepare_data(data_list_select, group_attr_id)

# # sampling parameters
n_keep = 10000
n_warmup = 5000
n_chain = 2

# fit model
data_stan = data_all[[1]]
data_gibbs = data_all[[2]]


# stan
stan_fit = stan(file='1d_HBM_multi_attri.stan', data = data_stan, iter = n_keep+n_warmup, warmup=n_warmup, chains = n_chain)

# gibbs
y = data_gibbs$y    # mean=4514, var=3724356. Check out mean and variance of y in each group.
exact_mean = sapply(y, cal_exact_mean)
exact_var = sapply(y, cal_exact_var)
gibbs_fit = fit_gibbs(data_gibbs, n_keep, n_warmup, n_chain)
index_good = (n_warmup+1):(n_warmup+n_keep)

# value of performance metrics
metric_all = cal_metric(stan_fit, gibbs_fit, y, print_out=1)

metric_all


