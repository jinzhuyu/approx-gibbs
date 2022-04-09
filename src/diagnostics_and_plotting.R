
library('bayesplot')
library('ggplot2')
library('LaplacesDemon')    # ESS function
library('mcmcse')    # ess function

# functions to generate plots

# calculate effective sample size
cal_rho = function(m, n, phi, t, var_phi_hat){
  
  # calculate variogram at each lag t, denoted by V_t
  V_t = 1/(m*(n-t))*sum((phi[(t+1):n, ] - phi[1:(n-t), ])^2)
  rho_t_hat = 1 - V_t/(2*var_phi_hat)
  
  return(rho_t_hat)
}


cal_n_eff = function(phi){
  
  # Ref.: BDA, Gelman et al. pp 284 to 287.
  # calculate the effective sample size given m chains with each chain having n samples.
  
  # calculate variance
  dim_phi = dim(phi)
  n = dim_phi[1]
  m = dim_phi[2]
  
  # between sequence variance
  phi_j_bar = colMeans(phi)
  phi_bar = mean(phi_j_bar)
  var_B = n/(m-1)*sum(phi_j_bar - phi_bar)
  
  # within sequence variance
  var_W = mean(apply(phi, 2, var))
  
  # marginal posterior variance (total variance)
  var_phi_hat = (n-1)/n*var_W + var_B/n
  

   
  # identify T: from 0 to first odd positive integer for which rho_T_hat +rho_T+1_hat is negative
  rho_hat_list = rep(0, n)
  T = n  
  t=1

  while (t<=n-3){
    for (i in 1:2){
      rho_hat_list[t+i] = cal_rho(m,n,phi,t+i,var_phi_hat)
    }
    if (sum(rho_hat_list[(t+1):(t+2)]) < 0){
      T = t
      break
    } else{
      t = t + 2 
    }
  }
    
    
  # effective sample size
  n_eff_hat = round(m*n/(1+2*sum(rho_hat_list[1:T])),0)
  
  return(n_eff_hat)
}

# error metric

cal_y_hat_stan = function(coeff_stan_list, data_stan){
  # usage: cal_y_hat_stan(coeff_stan_list,data_stan)
  x = data_stan$X
  len_each_group = data_stan$len_each_group
  coeff_mean = matrix(NA, n_group, n_cov)
  for (i in 1:n_cov) {
    for (j in 1:n_group){
      coeff_mean[j,i] = mean(coeff_stan_list[,j,i])
    }
  }
  
  coeff_mean_rep = matrix(NA,nrow(x), n_cov)
  for(i in 1:n_cov){
    for (j in 1:n_group) {
      if (j==1){
        row_id_curr = seq(1,len_each_group[1])
        coeff_mean_rep[row_id_curr,i] = rep(coeff_mean[j,i], each=len_each_group[j])
      } else {
        row_id_curr = seq(sum(len_each_group[1:(j-1)])+1, sum(len_each_group[1:j]))
        coeff_mean_rep[row_id_curr,i] = rep(coeff_mean[j,i], each=len_each_group[j])
      }
    }
  }
  
  y_hat_stan = exp(rowSums(coeff_mean_rep*x))
  
  return(y_hat_stan)
}

cal_y_hat_gibbs = function(coeff_gibbs, data_gibbs){
  
  x = data_gibbs$X
  len_each_group = data_gibbs$len_each_group
  
  n_coeff = dim(coeff_gibbs)[1]
  n_group = dim(coeff_gibbs)[3]
  n_meth = dim(coeff_gibbs)[5]
  coeff_mean = array(NA, dim=c(n_coeff, n_group, n_meth))
  for (i_coeff in 1:n_coeff) {
    for (i_group in 1:n_group) {
      for (i_meth in 1:n_meth) {
        coeff_mean[i_coeff, i_group, i_meth] = mean(coeff_gibbs[i_coeff, index_good, i_group, ,i_meth]) 
      }
    }
  }
  
  coeff_mean_rep = matrix(NA, nrow(x), n_cov)
  for(i in 1:dim(coeff_mean)[1]){
    for (j in 1:n_group) {
      if (j==1){
        row_id_curr = seq(1,len_each_group[1])
        coeff_mean_rep[row_id_curr,i] = rep(coeff_mean[i,j,], each=len_each_group[j])
      } else {
        row_id_curr = seq(sum(len_each_group[1:(j-1)])+1, sum(len_each_group[1:j]))
        coeff_mean_rep[row_id_curr,i] = rep(coeff_mean[i,j,], each=len_each_group[j])
      }
    }
  }
  
  y_hat_gibbs = exp(rowSums(coeff_mean_rep*x))
  
  return(y_hat_gibbs)
}


r_squared = function(y, y_hat){
  total_var = sum((y-mean(y))^2)
  resi_sum_squared_err = sum((y-y_hat)^2)
  r_squared = 1 - resi_sum_squared_err/total_var
  
  return(r_squared)
}


rmse = function(y, y_hat){
  rmse = sqrt(sum((y-y_hat)^2)/length(y))
  
  return(rmse)
}


cal_n_eff_stan = function(stan_fit_samples){
  # usage: cal_n_eff_stan(stan_fit_samples=coeff_stan)
  n_param = length(stan_fit_samples)
  n_eff_list = vector(mode = "list", length = n_param)
  
  for (i in 1:n_param){
    sample_temp = stan_fit_samples[[i]]
    dim_sample = dim(sample_temp)
    if (length(dim_sample)==3){
      n_eff_list[[i]] = matrix(0,dim_sample[2],dim_sample[3])
      for (j in 1:dim_sample[2]) {
        for (k in 1:dim_sample[3]){
          sample_temp_curr = matrix(as.vector(sample_temp[,j,k]), n_keep, n_chain)
          n_eff_list[[i]][j,k] = cal_n_eff(sample_temp_curr)
        }
      }
    }
    if (length(dim_sample)==2){
      n_eff_list[[i]] = rep(0,dim_sample[2])
      for (j in 1:dim_sample[2]) {
          sample_temp_curr = matrix(as.vector(sample_temp[,j]), n_keep, n_chain)
          n_eff_list[[i]][j] = cal_n_eff(sample_temp_curr)
      }
    }
  }
  
  return(n_eff_list)
}


cal_n_eff_gibbs = function(coeff_gibbs){
  n_coeff = dim(coeff_gibbs)[1]
  e_eff_gibbs = array(NA, dim = c(n_group, n_meth, n_coeff))
  for (i_meth in 1:n_meth){
    for (i_coeff in 1:n_coeff){
      for (i_group in 1:n_group){
        # cat('\n group id:', i_group, '\n')
        sample_temp_curr = matrix(as.vector(coeff_gibbs[i_coeff,index_good,i_group,,i_meth]),
                                  n_keep, n_chain)
        e_eff_gibbs[i_group,i_meth,i_coeff] = cal_n_eff(sample_temp_curr)
      }
    }
  }
  
  return(e_eff_gibbs)
}


cal_ESS_gibbs = function(coeff_gibbs){
  
  ESS_gibbs = array(NA, dim = c(n_group, n_meth, n_coeff))
  for (i_meth in 1:n_meth){
    for (i_coeff in 1:n_coeff){
      for (i_group in 1:n_group){
        # cat('\n group id:', i_group, '\n')
        sample_temp_curr = matrix(as.vector(coeff_gibbs[i_coeff,index_good,i_group,,i_meth]),
                                  n_keep, n_chain)
        ESS_gibbs[i_group,i_meth,i_coeff] = sum(ESS(sample_temp_curr))
      }
    }
  }
  
  return(ESS_gibbs)
}

cal_metric = function(stan_fit, gibbs_fit, y, print_out=1){
  ## execution time
  # stan
  n_iter = 1000    # time per 1000 iterations
  
  time_stan_all = get_elapsed_time(stan_fit)
  time_stan_sample = sum(time_stan_all[,2])
  mean_time_stan = time_stan_sample/(n_keep*n_chain/n_iter)
  
  # Approximate Gibbs
  time_gibbs_all = gibbs_fit[[4]]
  time_gibbs_sample = time_gibbs_all*n_keep/(n_keep+n_warmup)
  mean_time_gibbs = time_gibbs_sample/(n_keep*n_chain/n_iter)
  
  
  ## sampling efficiency
  # stan
  coeff_stan = rstan::extract(stan_fit, pars=c('coeff'))
  n_eff_stan = cal_n_eff_stan(coeff_stan)
  eff_stan = mean(n_eff_stan[[1]])/time_stan_sample
  
  # Approx. Gibbs
  coeff_gibbs = gibbs_fit[[1]]
  n_eff_gibbs = cal_n_eff_gibbs(coeff_gibbs)
  eff_gibbs = mean(n_eff_gibbs)/time_gibbs_sample
  
  # error measure
  coeff_stan_list = coeff_stan$coeff
  y_hat_stan = cal_y_hat_stan(coeff_stan_list,data_stan)
  r2_stan = r_squared(y, y_hat_stan)
  rmse_stan = rmse(y, y_hat_stan)
  
  y_hat_gibbs = cal_y_hat_gibbs(coeff_gibbs, data_gibbs)
  r2_gibbs = r_squared(y, y_hat_gibbs)
  rmse_gibbs = rmse(y, y_hat_gibbs)
  
  # print
  if(print_out==1){
    # execution time
    cat('Time for each chain of samples with Stan: ', time_stan_all[,2],'\n')
    # cat('Total time for 4 chain of samples with Stan: ', time_stan_sample,'\n')
    cat('Average time per 1000 iterations with Stan: ', mean_time_stan,'\n')
    cat('Average time per 1000 iterations with Approx. Gibbs: ', mean_time_gibbs,'\n')
    
    # effeciency
    cat('Sampling efficiency for stan and Approx. Gibbs: ', eff_stan, eff_gibbs,'\n')
    
    # error measure
    cat('r2 for stan and approx. Gibbs: ', r2_stan, r2_gibbs,'\n')
    cat('rmse for stan and approx. Gibbs: ', rmse_stan, rmse_gibbs,'\n')
  }
  
  return(c(mean_time_stan, mean_time_gibbs, eff_stan, eff_gibbs,
              r2_stan, r2_gibbs,rmse_stan, rmse_gibbs))
}


# calculated the R2 of gibbs only
cal_r2_gibbs = function(gibbs_fit){
  coeff_gibbs = gibbs_fit[[1]]
  y_hat_gibbs = cal_y_hat_gibbs(coeff_gibbs, data_gibbs)
  r2_gibbs = r_squared(y, y_hat_gibbs)
  
  return(r2_gibbs)
}


###########################3
# time complexity

plot_time = function(x, y, x_lab, y_lab){
  
  # fit linear model
  lm_fit = lm(y~x)
  y_fit = predict(lm_fit)
  par(mar =  c(6, 6, 4, 4)) 
  plot(x, y, type='p', ylab=y_lab, font.lab = 2,
       xlab=x_lab, col='black', pch=16, cex=1)
  lines(x, y_fit, col='darkred', lty=1, lwd=2, font.lab = 2)
  legend("bottomright", legend=c("True value", "Fitted value"), bty = "n",
         col=c("black", "darkred"), pch=c(16,NA), lty=c(NA,1), cex=rep(1,2), lwd=c(NA,2))
  
  # add R2
  modsum = summary(lm_fit)
  r2 = modsum$adj.r.squared
  mylabel = bquote(italic(R)^2 == .(format(r2, digits = 4)))
  text(x = max(x)*0.925, y = (max(y))*0.725, labels = mylabel, cex=1.0)
}


################################################################################
##### 3 - Summarize and visualize posterior distributions 
################################################################################


plot_colors <<- c('blue','red', 'pink', "orange")
legend_names = sapply(1:n_chain, toString)
n_total_sample <<- n_warmup+n_keep
index_good = (n_total_sample*0.5):n_total_sample
# if(max(coeff)<=1e-4){
#   axis(2,at=marks,labels=format(marks,scientific=T))
# }
plot_sample = function(coeff, coeff_name) {
  for (i_meth in 1:n_meth) {
    for (i_coeff in 1:length(coeff_name)){
      par(mfrow=c(2, ceiling(n_group/2)))
      for (i_group in 1:n_group){
        for(i_chain in 1:n_chain){
          if (i_chain==1){
            plot(index_good, coeff[i_coeff, index_good, i_group,i_chain, i_meth], type='l',
                 ylab=paste(coeff_name[i_coeff],i_group), xlab='', col=plot_colors[i_chain])
          } else{
            lines(index_good, coeff[i_coeff, index_good, i_group,i_chain, i_meth],
                  ylab=paste(coeff_name[i_coeff],i_group), xlab= '', col=plot_colors[i_chain])
          }
        }
      }
      legend("bottom", legend=legend_names, col=plot_colors, lwd=1.0,
             cex=1.0, bty='n', xpd =NA, horiz=T, inset=c(0,-0.5))
    }
  }
}


plot_density = function(coeff, coeff_name) {
  for (i_meth in 1:n_meth) {
    for (i_coeff in 1:length(coeff_name)){
      par(mfrow=c(2, ceiling(n_group/2)))
      for (i_group in 1:n_group){
        for(i_chain in 1:n_chain){
          if (i_chain==1){
            plot(density(coeff[i_coeff, index_good, i_group,i_chain, i_meth]), type='l', main = '',
                 xlab=paste(coeff_name[i_coeff],i_group), ylab='Density', col=plot_colors[i_chain])
          } else{
            lines(density(coeff[i_coeff, index_good, i_group,i_chain, i_meth]), main = '',
                  xlab=paste(coeff_name[i_coeff],i_group), ylab='Density', col=plot_colors[i_chain])
          }
        }
      }
      legend("bottom", legend=legend_names, col=plot_colors, lwd=1.0,
             cex=1.0, bty='n', xpd = NA, horiz = T, inset = c(0,-0.5))
    }
  }
}


