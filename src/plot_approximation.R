# show quality of Gaussian approximation to the log-gamma distribution
library(tidyr)
library(plyr)
library(ggplot2)
theme_set(theme_bw())
source('diagnose_and_plot.R')



cal_phi_0 = function(y){
  if(y<1){
    stop('Error: the input must be >=1 in calculating the exact mean of the approximate distribution')
  }
  euler_con = 0.57721566
  if (y==1) {
    phi_0 = -euler_con
  }else{
    k = seq(1, y-1)
    harm_num = sum(1/k)    # harmonic number
    phi_0 = -euler_con + harm_num
  }
  return(phi_0)
}


cal_phi_1 = function(y){
  if(y<1){
    stop('Error: the input must be no smaller than 1 in calculating the exact var of the approximate distribution')
  }
  if(y==1){
    phi_1 = pi^2/6
  }else{
    seq_neg_part = seq(1,y-1)
    phi_1 = pi^2/6 - sum(1/seq_neg_part^2)
  }
  return(phi_1)
}


gen_log_gamma = function(v, y){
  return (1/factorial(y-1)*exp(v*y)*exp(-exp(v)))
}


gen_approx_norm = function(v, y){
  phi_1 = cal_phi_1(y)
  phi_0 = cal_phi_0(y)
  return ( 1/(sqrt(2*pi*phi_1))*exp((v-phi_0)^2/(-2*phi_1)) )
}

gen_df = function(y){
  intv = 0.025
  if (y == 1){
    v = seq(-3.25, 3.25, by=intv)
  } else if (y==2){
    v = seq(-2, 3.2, by=intv)  
  } else if (y==3){
    v = seq(-1.15, 3.25, by=intv)  
  } else if (y==5){
    v = seq(0, 3.5, by=intv)
  } else if (y==10){
    v = seq(0.75, 3.8, by=intv)
  } else if (y==20){
    v = seq(1.95, 4.0, by=intv)
  } else if (y==40){
    v = seq(3, 4.35, by=intv)
  } else {
    print('This y value is not needed')
  }
  Approximate = gen_approx_norm(v, y)
  True = gen_log_gamma(v, y)
  df = data.frame(v, Approximate, True)
  return (df)
}

plot_approx = function(df, save_fig_title, y_current){
  # stack vertically to feed into ggplot
  df = df%>%gather(key, value, True, Approximate)  
  
  ggplot(df, aes(x=v, y =value, color=key, group=key)) + 
    geom_line(aes(linetype=key), size=0.65) +
    labs(x = 'v', y = 'Probability density', size=2.2) +
    scale_linetype_manual(values = c("dashed", "solid")) +
    scale_color_manual(values=c('dark blue', 'dark red')) +
    annotation_compass(paste('y = ', y_current, spe=''),'E')+
    theme(axis.text = element_text(size=12), 
          axis.title = element_text(size=13),
          axis.title.y = element_text(margin = margin(t=0, r=3, b=0, l=0)),
          plot.margin = unit(c(0.2,0.4,0.2,0.2), "cm"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position=c(0.78,0.88),
          legend.background=element_rect(fill="transparent",colour=NA),
          legend.spacing.y = unit(0, "mm"),
          legend.title = element_blank(),
          legend.text = element_text(size=11),
          legend.key = element_rect(colour = NA, fill = NA))
  ggsave(paste(save_fig_title, '.pdf', sep=''), plot = last_plot(), width=4.25, height=3, dpi=300)
  ggsave(paste(save_fig_title, '.png', sep=''), plot = last_plot(), width=4.25, height=3, dpi=600)
}


main_plot_approx = function(){
  # setwd('C:/code/approx-gibbs/results')
  y_list = c(1, 2, 3, 5, 10, 20,40)
  for (y in y_list){
    df = gen_df(y)
    save_fig_title = paste("../results/approx_vs_real_y_", y, sep='')
    plot_approx(df, save_fig_title, y)
  }
}


cal_dist = function(y){
  n = 1e7
  
  approx_norm = rnorm(n, cal_phi_0(y), sqrt(cal_phi_1(y)))
  log_gamma = log(rgamma(n, shape=y, scale = 1))
  
  ks_dist =  ks.test(approx_norm, log_gamma)$statistic[['D']]
  abs_diff_mean = abs(mean(approx_norm) - mean(log_gamma))
  return (c(ks_dist, abs_diff_mean))
}


main_plot_dist = function(){
  y_list = seq(1, 15, by=1)
  dist_list = sapply(y_list, cal_dist)
  ks_value_scale = 100
  ks_list = dist_list[1,]*ks_value_scale
  diff_mean_value_scale = ks_value_scale * 10
  abs_diff_mean_list = dist_list[2,]*diff_mean_value_scale
  
  df = data.frame(y_list, ks_list, abs_diff_mean_list)
  
  y_scale_diff_factor = 3 # transform data on right y axis
  ggplot(df, aes(x=y_list)) +
    geom_line( aes(y=ks_list), size=0.7, color='darkred') + 
    geom_line( aes(y=abs_diff_mean_list * y_scale_diff_factor), size=0.7, color='darkblue',linetype =2) +
    scale_y_continuous(
      # Features of the first axis
      name = bquote("KS distance (x"~10^{-2}~')'),
      # Add a second axis and specify its features
      sec.axis = sec_axis(~./y_scale_diff_factor, name=bquote("Abs. error of the mean (x"~10^{-3}~')') )
    ) +
    xlab('y') +
    scale_x_continuous(breaks = round(seq(min(y_list), max(y_list), by = 2),1)) + 
    theme(
      axis.text.x = element_text(size=13),
      axis.title.x = element_text(size=13),
      axis.title.y.left = element_text(color = 'darkred', size=13,
                                       margin = margin(t=0, r=4, b=0, l=1)),
      axis.title.y.right = element_text(color = 'darkblue', size=13,
                                        margin = margin(t=0, r=-6.5, b=0, l=4)),
      axis.text.y.left = element_text(color = 'darkred',size=13),
      axis.text.y.right = element_text(color = 'darkblue',size=13),
      plot.margin = unit(c(0.2,0.4,0.2,0.2), "cm"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),    )
  save_fig_title = '../results/dist'
  ggsave(paste(save_fig_title, '.pdf', sep=''), plot = last_plot(), width=4.25, height=3, dpi=300)
  ggsave(paste(save_fig_title, '.png', sep=''), plot = last_plot(), width=4.25, height=3, dpi=600)
}


main_plot_approx()

main_plot_dist()
