# show time complexity of nuts and ags
library(tidyr)
library(plyr)
library(ggplot2)
theme_set(theme_bw())


plot_scalability = function(df_select, lgd_x, lgd_y_vec,save_fig_title){

  # tack vertically to feed into ggplot
  df_select =  df_select%>%
    gather(key, value, NUTS, AGS)  
  
  # calculate r2
  r2<-ddply(df_select,.(key), function(x) summary(lm(x$value ~ 0+x$Nd.Nc))$r.squared)
  r2$V1 = round(r2$V1, 4)
  names(r2) = c('key','r2')
  
  # use ggplot2 to make plots
  pdf(file=save_fig_title, width=4.2, height=3)
  ggplot(df_select, aes(x=Nd.Nc, y =value, color=key, group=key)) + 
    # geom_line(aes(linetype=key), size=0.75) + 
    geom_point(aes(shape=key), size=1.65) +
    geom_smooth(formula=y~0+x, method = "lm", fill = NA, se = FALSE, aes(linetype=key), size=0.8) +
    labs(x = "Size of dataset", y = 'Average sampling time\nper 1000 iter. (s)', size=2) +
    scale_color_manual(values=c('dark red', 'dark blue')) +
    geom_text(data=r2,aes(label = paste("R^2: ", r2, sep="")),parse=T,x=lgd_x,y=lgd_y_vec, show.legend = F)+
    theme(axis.text = element_text(size=12), 
          axis.title = element_text(size=13, face="bold"),
          axis.title.y = element_text(margin = margin(t = 0, r = 5, b = 0, l = 0)),
          plot.margin = unit(c(0.2,0.7,0.2,0.2), "cm"),
          legend.position=c(0.175,0.829),
          legend.background=element_rect(fill="transparent",colour=NA),
          legend.spacing.y = unit(0, "mm"),
          legend.title = element_blank(),
          legend.text = element_text(size=11),
          legend.key = element_rect(colour = NA, fill = NA)) #,
  # legend.box.background = element_rect(colour = "black"))
  
  dev.off()
}


# real data
# import data
data = read.csv("../results/scalability_real_data.csv", header = T)

df_select = data[,c(2,6,7)]

# define plot parameters
lgd_x = 16950
lgd_y_vec = c(110,615)
save_fig_title = "../results/scalability_real_data.pdf"
plot_scalability(df_select, lgd_x, lgd_y_vec, save_fig_title)


# simualted data
# import data
data = read.csv("../results/scalability_simulated_data.csv", header = T)
df_select = data[,c(2,7,8)]

# define plot parameters
lgd_x = 3900
lgd_y_vec = c(12,77)
save_fig_title = "../results/scalability_simulated_data.pdf"
plot_scalability(df_select, lgd_x, lgd_y_vec, save_fig_title)

