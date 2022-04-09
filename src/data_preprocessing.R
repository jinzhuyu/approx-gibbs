
rep_row = function(x_vec,n){
  
  return(matrix(rep(x_vec,each=n),nrow=n))
}

min_max = function(x_vec,lb_norm=0.005){
  
  range = max(x_vec) - min(x_vec)
  x_norm = (x_vec-min(x_vec)+lb_norm*range)/((1+lb_norm)*range)
  
  return(x_norm)
}

find_group_id = function(data_list, group_attr_id, n_data){
  
  group_attr = data_list[[group_attr_id]]
  group_id = rep(1, n_data)
  id_temp = 1
  for (i in 2:n_data){
    if(group_attr[i]!=group_attr[i-1]){
      id_temp = id_temp + 1
    }
    group_id[i] = id_temp
  }
  
  return(group_id)
}


find_len_each_group = function(data_list, group_attr_id, n_data, n_group){
  
  group_attr = data_list[[group_attr_id]]
  len_each_group = rep(1, n_group)
  
  group_id_temp = 1
  for (i in 2:n_data) {
    if(group_attr[i]==group_attr[i-1]){
      len_each_group[group_id_temp] = len_each_group[group_id_temp] + 1
    }else{
      group_id_temp = group_id_temp + 1
      # len_each_group[group_id_temp] = 1
    }
  }
  
  return(len_each_group)
}


prepare_data = function(data_list, group_attr_id){
  
  len_list = length(data_list)
  y = as.numeric(unlist(data_list[len_list])) # convert to numeric values, otherwise R will crash at the start of Stan sampling.
  X_list = data_list[-len_list]
  # X_list[group_attr_id] = log(X_list[group_attr_id])
  n_data = length(y)
  n_coeff = length(X_list)
  X = matrix(as.numeric(unlist(X_list)), n_data, n_coeff)
  
  group_id = find_group_id(data_list, group_attr_id, n_data)
  n_group = max(group_id)
  
  len_each_group = find_len_each_group(data_list, group_attr_id, n_data, n_group)
  
  data_stan = list(y = y,  X = X, n_data = n_data, n_group=n_group, n_coeff=n_coeff, group_id=group_id,
                   len_each_group=len_each_group)
  data_gibbs = list(y = y,  X = X, n_data = n_data, n_group=n_group, n_coeff=n_coeff, group_id=group_id,
                    len_each_group=len_each_group)
  
  return(list(data_stan, data_gibbs))
}


fit_model = function(data_all, group_attr_id){
  
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
  
  return(metric_all)
}








  