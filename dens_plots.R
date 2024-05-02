setwd("~/GitHub/environment-modelling")

library(cond.extremes)
library(gridExtra)
library(patchwork)
source("~/GitHub/env-contours/FORM_functions_revised.R")

## read in data
data <- read.csv("data/cnsTS.txt")

## isolate peaks

hs_peak <- c()
t2_peak <- c()
for (i in 1:max(data$StrIdn)){
  hs_peak[i] <- max(data$Hs[data$StrIdn == i])
  t2_peak[i] <- max(data$T2[data$StrIdn == i])
}

stp_peak <- (2 * pi * hs_peak) / (9.81 * t2_peak^2)

# load dens ---------------------------------------------------

dens = read.csv("data/env_probs.csv")

# plotting conditioned density --------------------------------------------

cond_dens <- read.csv("data/cond_dens_fixed.csv", header=F)

cond_dens_df = data.frame(
  hs = dens$x,
  stp = dens$y,
  cd = cond_dens
)

addCntr=T
justCntr=F
lwd=1.2

plot1 = ggplot(cond_dens_df, aes(x=hs, y=stp)) + geom_raster(aes(fill=((V1)))) +
  scale_fill_gradientn(colours=c("white", "yellow", "orange", "red", "black"),
                       limits = c(0,  37.8)) +
  labs(fill="CDE")

plot1 = plot1 +
  ylim(layer_scales(plot1)$y$range$range[1], layer_scales(plot1)$y$range$range[2]) +
  theme_linedraw() +
  scale_x_continuous(limits = c(7.5, 22.5), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0.04, 0.08), expand = c(0, 0)) +
  theme(legend.position = "none", panel.background = element_rect(fill = NA),
        panel.ontop = TRUE, legend.title = element_text(size=20), legend.text = element_text(size=15), legend.key.size = unit(1.5, 'cm'),
        axis.text=element_text(size=15), axis.title=element_text(size=20)) + xlab("Hs") + ylab("S2")

if(addCntr){
  
  if(justCntr){
    plot1 = ggplot() +
      ylim(layer_scales(plot1)$y$range$range[1], layer_scales(plot1)$y$range$range[2]) +
      theme_linedraw() +
      scale_x_continuous(limits = c(0, 20), expand = c(0, 0)) +
      scale_y_continuous(limits = c(0.001, 0.09), expand = c(0, 0)) +
      theme( panel.background = element_rect(fill = NA),
             panel.ontop = TRUE, legend.title = element_text(size=20), legend.text = element_text(size=15), legend.key.size = unit(1.5, 'cm'),
             axis.text=element_text(size=15), axis.title=element_text(size=20)) +
      xlab("Hs") + ylab("S2")
  }

  # data --------------------------------------------------------------------
  sample.df = data.frame(
    x = hs_peak,
    y = stp_peak
  )
  
  plot1 = plot1 + geom_point(data=sample.df, aes(x=x, y=y), col="grey") 
  # read in contours ---------------------------------------------------------
  cntr_path = "~/GitHub/env-contours/FORM fits/"

  contours=c("neggev_exp_exp_form_p1.38888888888889e-05", "posweibull_lin_qua_form_p1.38888888888889e-05", "neglnorm_lin_lin_form_p1.38888888888889e-05" ,  "neggamma_qua_lin_form_p1.38888888888889e-05"   ,
             "negweibull_exp_exp_form_p1.38888888888889e-05" ,  "posgamma_lin_lin_form_p1.38888888888889e-05" ,    "poslnorm_qua_exp_form_p1.38888888888889e-05" ,  "posgev_lin_qua_form_p1.38888888888889e-05") 
  # AICs = c()
  # for (i in 1:length(contours)){
  #   #cntr_name = substr(contours[i], 1, nchar(contours[i])-12)
  #   cntr_name = substr(contours[i], 1,  nchar(contours[i])-12)
  #   load(file = paste("~/GitHub/env-contours/FORM_AICs/", cntr_name, "_AIC", sep=""))
  #   AICs[i] = AIC
  # }

  #contours_w_aic =c("S'~GEV" , "S~Weibull" , "S'~Lognormal", "S'~Gamma", "S'~Weibull", "S~Gamma", "S~Lognormal", "S~GEV")
  contours_w_aic =c("C1" , "C2" , "C3", "C4", "C5", "C6", "C7", "C8")
  
  for (cntr in contours){
    cntr_file = load(paste(cntr_path, cntr, sep=""))
    assign(cntr, get(cntr_file))
    assign(cntr, data.frame(x= get(cntr)$x, y=get(cntr)$y))
  }
  
  # adding contours ---------------------------------------------------------
  # pallette = c("red", "orange", "yellow4", "lightgreen",
  #  "green", "lightblue", "blue", "purple",
  # "violet", "black")
  pallette = c("red", "orange", "yellow4")

  plot1 = plot1 + geom_path(data=get(contours[1]), aes(x=x, y=y, col = contours_w_aic[1]), linewidth=lwd) +
    geom_path(data=get(contours[2]), aes(x=x, y=y, col = contours_w_aic[2]), linewidth=lwd) +
    geom_path(data=get(contours[3]), aes(x=x, y=y, col = contours_w_aic[3]), linewidth=lwd) #+ 
    # geom_path(data=get(contours[4]), aes(x=x, y=y, col = contours_w_aic[4])) +
    # geom_path(data=get(contours[5]), aes(x=x, y=y, col = contours_w_aic[5])) +
    # geom_path(data=get(contours[6]), aes(x=x, y=y, col = contours_w_aic[6])) +
    # geom_path(data=get(contours[7]), aes(x=x, y=y, col = contours_w_aic[7])) +
    # geom_path(data=get(contours[8]), aes(x=x, y=y, col = contours_w_aic[8]))

  plot1 = plot1 + scale_color_manual(name="", breaks=contours_w_aic,
                                   values=pallette) 
}

# plotting second conditioned density -------------------------------------

cond_dens <- read.csv("data/cond_dens_pos.csv", header=F)
# area =  0.0001296296
# cond_dens = cond_dens/sum(cond_dens*area)

cond_dens_df = data.frame(
  hs = dens$x,
  stp = dens$y,
  cd = cond_dens
)

addCntr=T
justCntr=F
lwd=1.2

plot2 = ggplot(cond_dens_df, aes(x=hs, y=stp)) + geom_raster(aes(fill=((V1)))) +
  scale_fill_gradientn(colours=c("white", "yellow", "orange", "red", "black"),
                       limits = c(0,  37.8)) +
  labs(fill="CDE")

plot2 = plot2 +
  ylim(layer_scales(plot2)$y$range$range[1], layer_scales(plot2)$y$range$range[2]) +
  theme_linedraw() +
  scale_x_continuous(limits = c(7.5, 22.5), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0.04, 0.08), expand = c(0, 0)) +
  theme(legend.position = "none", panel.background = element_rect(fill = NA),
        panel.ontop = TRUE, legend.title = element_text(size=20), legend.text = element_text(size=15), legend.key.size = unit(1.5, 'cm'),
        axis.text=element_text(size=15), axis.title=element_text(size=20)) +
  
  xlab("Hs") + ylab("S2")


if(addCntr){
  
  if(justCntr){
    plot2 = ggplot() +
      ylim(layer_scales(plot2)$y$range$range[1], layer_scales(plot2)$y$range$range[2]) +
      theme_linedraw() +
      scale_x_continuous(limits = c(0, 20), expand = c(0, 0)) +
      scale_y_continuous(limits = c(0.001, 0.09), expand = c(0, 0)) +
      theme( panel.background = element_rect(fill = NA),
             panel.ontop = TRUE, legend.title = element_text(size=20), legend.text = element_text(size=15), legend.key.size = unit(1.5, 'cm'),
             axis.text=element_text(size=15), axis.title=element_text(size=20)) +
      xlab("Hs") + ylab("S2")
  }
  
  # data --------------------------------------------------------------------
  sample.df = data.frame(
    x = hs_peak,
    y = stp_peak
  )
  
  plot2 = plot2 + geom_point(data=sample.df, aes(x=x, y=y), col="grey") 
  # read in contours ---------------------------------------------------------
  cntr_path = "~/GitHub/env-contours/FORM fits/"
  
  contours=c("neggev_exp_exp_form_p1.38888888888889e-05", "posweibull_lin_qua_form_p1.38888888888889e-05", "neglnorm_lin_lin_form_p1.38888888888889e-05" ,  "neggamma_qua_lin_form_p1.38888888888889e-05"   ,
             "negweibull_exp_exp_form_p1.38888888888889e-05" ,  "posgamma_lin_lin_form_p1.38888888888889e-05" ,    "poslnorm_qua_exp_form_p1.38888888888889e-05" ,  "posgev_lin_qua_form_p1.38888888888889e-05") 
  # AICs = c()
  # for (i in 1:length(contours)){
  #   #cntr_name = substr(contours[i], 1, nchar(contours[i])-12)
  #   cntr_name = substr(contours[i], 1,  nchar(contours[i])-12)
  #   load(file = paste("~/GitHub/env-contours/FORM_AICs/", cntr_name, "_AIC", sep=""))
  #   AICs[i] = AIC
  # }
  
  #contours_w_aic =c("S'~GEV" , "S~Weibull" , "S'~Lognormal", "S'~Gamma", "S'~Weibull", "S~Gamma", "S~Lognormal", "S~GEV")
  contours_w_aic =c("C1", "C2" , "C3", "C4", "C5", "C6", "C7", "C8")
  
  
  for (cntr in contours){
    cntr_file = load(paste(cntr_path, cntr, sep=""))
    assign(cntr, get(cntr_file))
    assign(cntr, data.frame(x= get(cntr)$x, y=get(cntr)$y))
  }
  
  # adding contours ---------------------------------------------------------
  # pallette = c("red", "orange", "yellow4", "lightgreen",
  #  "green", "lightblue", "blue", "purple",
  # "violet", "black")
  pallette = c("red", "orange", "yellow4")
  
  plot2 = plot2 + geom_path(data=get(contours[1]), aes(x=x, y=y, col = contours_w_aic[1]), linewidth=lwd) +
    geom_path(data=get(contours[2]), aes(x=x, y=y, col = contours_w_aic[2]), linewidth=lwd) +
    geom_path(data=get(contours[3]), aes(x=x, y=y, col = contours_w_aic[3]), linewidth=lwd) #+ 
  # geom_path(data=get(contours[4]), aes(x=x, y=y, col = contours_w_aic[4])) +
  # geom_path(data=get(contours[5]), aes(x=x, y=y, col = contours_w_aic[5])) +
  # geom_path(data=get(contours[6]), aes(x=x, y=y, col = contours_w_aic[6])) +
  # geom_path(data=get(contours[7]), aes(x=x, y=y, col = contours_w_aic[7])) +
  # geom_path(data=get(contours[8]), aes(x=x, y=y, col = contours_w_aic[8]))
  
  plot2 = plot2 + scale_color_manual(name="", breaks=contours_w_aic,
                                   values=pallette) 
  
}


# plotting third conditioned density -------------------------------------

cond_dens <- read.csv("data/cond_dens_neg.csv", header=F)
# area =  0.0001296296
# cond_dens = cond_dens/sum(cond_dens*area)

cond_dens_df = data.frame(
  hs = dens$x,
  stp = dens$y,
  cd = cond_dens
)

addCntr=T
justCntr=F
lwd=1.2

plot3 = ggplot(cond_dens_df, aes(x=hs, y=stp)) + geom_raster(aes(fill=((V1)))) +
  scale_fill_gradientn(colours=c("white", "yellow", "orange", "red", "black"),
                       limits = c(0,  37.8)) +
  labs(fill="CDE")

plot3 = plot3 +
  ylim(layer_scales(plot3)$y$range$range[1], layer_scales(plot3)$y$range$range[2]) +
  theme_linedraw() +
  scale_x_continuous(limits = c(7.5, 22.5), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0.04, 0.08), expand = c(0, 0)) +
  theme(legend.position = "none", panel.background = element_rect(fill = NA),
        panel.ontop = TRUE, legend.title = element_text(size=20), legend.text = element_text(size=15), legend.key.size = unit(1.5, 'cm'),
        axis.text=element_text(size=15), axis.title=element_text(size=20)) +
  
  xlab("Hs") + ylab("S2")

if(1-addCntr){
  
  print(plot3)
  
}

if(addCntr){
  
  if(justCntr){
    plot3 = ggplot() +
      ylim(layer_scales(plot3)$y$range$range[1], layer_scales(plot3)$y$range$range[2]) +
      theme_linedraw() +
      scale_x_continuous(limits = c(0, 20), expand = c(0, 0)) +
      scale_y_continuous(limits = c(0.001, 0.09), expand = c(0, 0)) +
      theme( panel.background = element_rect(fill = NA),
             panel.ontop = TRUE, legend.title = element_text(size=20), legend.text = element_text(size=15), legend.key.size = unit(1.5, 'cm'),
             axis.text=element_text(size=15), axis.title=element_text(size=20)) +
      xlab("Hs") + ylab("S2")
  }
  
  # data --------------------------------------------------------------------
  sample.df = data.frame(
    x = hs_peak,
    y = stp_peak
  )
  
  plot3 = plot3 + geom_point(data=sample.df, aes(x=x, y=y), col="grey") 
  # read in contours ---------------------------------------------------------
  cntr_path = "~/GitHub/env-contours/FORM fits/"
  
  contours=c("neggev_exp_exp_form_p1.38888888888889e-05", "posweibull_lin_qua_form_p1.38888888888889e-05", "neglnorm_lin_lin_form_p1.38888888888889e-05" ,  "neggamma_qua_lin_form_p1.38888888888889e-05"   ,
             "negweibull_exp_exp_form_p1.38888888888889e-05" ,  "posgamma_lin_lin_form_p1.38888888888889e-05" ,    "poslnorm_qua_exp_form_p1.38888888888889e-05" ,  "posgev_lin_qua_form_p1.38888888888889e-05") 
  # AICs = c()
  # for (i in 1:length(contours)){
  #   #cntr_name = substr(contours[i], 1, nchar(contours[i])-12)
  #   cntr_name = substr(contours[i], 1,  nchar(contours[i])-12)
  #   load(file = paste("~/GitHub/env-contours/FORM_AICs/", cntr_name, "_AIC", sep=""))
  #   AICs[i] = AIC
  # }
  
  #contours_w_aic =c("S'~GEV" , "S~Weibull" , "S'~Lognormal", "S'~Gamma", "S'~Weibull", "S~Gamma", "S~Lognormal", "S~GEV")
  contours_w_aic =c("C1" , "C2" , "C3", "C4", "C5", "C6", "C7", "C8")
  
  for (cntr in contours){
    cntr_file = load(paste(cntr_path, cntr, sep=""))
    assign(cntr, get(cntr_file))
    assign(cntr, data.frame(x= get(cntr)$x, y=get(cntr)$y))
  }
  
  # adding contours ---------------------------------------------------------
  # pallette = c("red", "orange", "yellow4", "lightgreen",
  #  "green", "lightblue", "blue", "purple",
  # "violet", "black")
  pallette = c("red", "orange", "yellow4")
  
  plot3 = plot3 + geom_path(data=get(contours[1]), aes(x=x, y=y, col = contours_w_aic[1]), linewidth=lwd) +
    geom_path(data=get(contours[2]), aes(x=x, y=y, col = contours_w_aic[2]), linewidth=lwd) +
    geom_path(data=get(contours[3]), aes(x=x, y=y, col = contours_w_aic[3]), linewidth=lwd) #+ 
  # geom_path(data=get(contours[4]), aes(x=x, y=y, col = contours_w_aic[4])) +
  # geom_path(data=get(contours[5]), aes(x=x, y=y, col = contours_w_aic[5])) +
  # geom_path(data=get(contours[6]), aes(x=x, y=y, col = contours_w_aic[6])) +
  # geom_path(data=get(contours[7]), aes(x=x, y=y, col = contours_w_aic[7])) +
  # geom_path(data=get(contours[8]), aes(x=x, y=y, col = contours_w_aic[8]))
  
  plot3 = plot3 + scale_color_manual(name="", breaks=contours_w_aic,
                                     values=pallette) 
  print(plot3)
  
}

# combining density plots ---------------------------------------------------------

combined = plot1 + plot2 + plot3 & theme(legend.position="bottom")
print(combined + plot_layout(guides="collect"))

# plotting first failure ps -----------------------------------------------------

fail_ps = read.csv("data/fixed_fail_ps.csv", header=F)

fail_ps_df = data.frame(
  hs=dens$x,
  stp=dens$y,
  p = fail_ps
)

plot4 = ggplot(fail_ps_df, aes(x=hs, y=stp)) + geom_raster(aes(fill=(log(V1)))) +
  scale_fill_gradientn(colours=c("white", "yellow", "orange", "red", "black"), na.value = "white", breaks=c(0, -10, -20, -30), limits=c(-30, 0)) +
  #ggtitle(paste("Probability of exceeding", as.integer(1/p), "Year Response - Weakpoint 5m Above Sea Level", sep=" ")) +
  labs(fill="Log exceedance \n probability") 
plot4 = plot4 + 
  ylim(layer_scales(plot4)$y$range$range[1], layer_scales(plot4)$y$range$range[2]) +
  scale_x_continuous(limits = c(0, 25), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0.01, 0.08), expand = c(0, 0)) +
  theme_linedraw() +
  theme( panel.background = element_rect(fill = NA),
         panel.ontop = TRUE, legend.title = element_text(size=20), legend.text = element_text(size=15), legend.key.size = unit(1.5, 'cm'),
         axis.text=element_text(size=15), axis.title=element_text(size=20)) +

  xlab("Hs") + ylab("S2")

# for adding contours 

if(addCntr){

  plot4 = plot4 + geom_point(data=sample.df, aes(x=x, y=y), col="grey")

  plot4 = plot4 + geom_path(data=get(contours[1]), aes(x=x, y=y, col = contours_w_aic[1]), linewidth=lwd) +
    geom_path(data=get(contours[2]), aes(x=x, y=y, col = contours_w_aic[2]), linewidth=lwd) +
    geom_path(data=get(contours[3]), aes(x=x, y=y, col = contours_w_aic[3]), linewidth=lwd) #+
    # geom_path(data=get(contours[4]), aes(x=x, y=y, col = contours_w_aic[4])) +
    # geom_path(data=get(contours[5]), aes(x=x, y=y, col = contours_w_aic[5])) +
    # geom_path(data=get(contours[6]), aes(x=x, y=y, col = contours_w_aic[6])) +
    # geom_path(data=get(contours[7]), aes(x=x, y=y, col = contours_w_aic[7])) +
    # geom_path(data=get(contours[8]), aes(x=x, y=y, col = contours_w_aic[8]))

  plot4 = plot4 + scale_color_manual(name="", breaks=contours_w_aic, values=pallette)

}




# plotting second failure ps -----------------------------------------------------

fail_ps = read.csv("data/pos_fail_ps.csv", header=F)

fail_ps_df = data.frame(
  hs=dens$x,
  stp=dens$y,
  p = fail_ps
)

plot5 = ggplot(fail_ps_df, aes(x=hs, y=stp)) + geom_raster(aes(fill=(log(V1)))) +
  scale_fill_gradientn(colours=c("white", "yellow", "orange", "red", "black"), na.value = "white", breaks=c(0, -10, -20, -30), limits=c(-30, 0)) +
  #ggtitle(paste("Probability of exceeding", as.integer(1/p), "Year Response - Weakpoint 5m Above Sea Level", sep=" ")) +
  labs(fill="Log exceedance \n probability") 
plot5 = plot5 + 
  ylim(layer_scales(plot5)$y$range$range[1], layer_scales(plot5)$y$range$range[2]) +
  scale_x_continuous(limits = c(0, 25), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0.01, 0.08), expand = c(0, 0)) +
  theme_linedraw() +
  theme( panel.background = element_rect(fill = NA),
         panel.ontop = TRUE, legend.title = element_text(size=20), legend.text = element_text(size=15), legend.key.size = unit(1.5, 'cm'),
         axis.text=element_text(size=15), axis.title=element_text(size=20)) +
  
  xlab("Hs") + ylab("S2")

# for adding contours 

if(addCntr){
  
  plot5 = plot5 + geom_point(data=sample.df, aes(x=x, y=y), col="grey")
  
  plot5 = plot5 + geom_path(data=get(contours[1]), aes(x=x, y=y, col = contours_w_aic[1]), linewidth=lwd) +
    geom_path(data=get(contours[2]), aes(x=x, y=y, col = contours_w_aic[2]), linewidth=lwd) +
    geom_path(data=get(contours[3]), aes(x=x, y=y, col = contours_w_aic[3]), linewidth=lwd) #+
  # geom_path(data=get(contours[4]), aes(x=x, y=y, col = contours_w_aic[4])) +
  # geom_path(data=get(contours[5]), aes(x=x, y=y, col = contours_w_aic[5])) +
  # geom_path(data=get(contours[6]), aes(x=x, y=y, col = contours_w_aic[6])) +
  # geom_path(data=get(contours[7]), aes(x=x, y=y, col = contours_w_aic[7])) +
  # geom_path(data=get(contours[8]), aes(x=x, y=y, col = contours_w_aic[8]))
  
  plot5 = plot5 + scale_color_manual(name="", breaks=contours_w_aic, values=pallette)

}

# plotting third failure ps -----------------------------------------------------

fail_ps = read.csv("data/neg_fail_ps.csv", header=F)

fail_ps_df = data.frame(
  hs=dens$x,
  stp=dens$y,
  p = fail_ps
)

plot6 = ggplot(fail_ps_df, aes(x=hs, y=stp)) + geom_raster(aes(fill=(log(V1)))) +
  scale_fill_gradientn(colours=c("white", "yellow", "orange", "red", "black"), na.value = "white", breaks=c(0, -10, -20, -30), limits=c(-30, 0)) +
  #ggtitle(paste("Probability of exceeding", as.integer(1/p), "Year Response - Weakpoint 5m Above Sea Level", sep=" ")) +
  labs(fill="Log exceedance \n probability") 
plot6 = plot6 + 
  ylim(layer_scales(plot6)$y$range$range[1], layer_scales(plot6)$y$range$range[2]) +
  scale_x_continuous(limits = c(0, 25), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0.01, 0.08), expand = c(0, 0)) +
  theme_linedraw() +
  theme( panel.background = element_rect(fill = NA),
         panel.ontop = TRUE, legend.title = element_text(size=20), legend.text = element_text(size=15), legend.key.size = unit(1.5, 'cm'),
         axis.text=element_text(size=15), axis.title=element_text(size=20)) +
  
  xlab("Hs") + ylab("S2")

# for adding contours 

if(addCntr){
  
  plot6 = plot6 + geom_point(data=sample.df, aes(x=x, y=y), col="grey")
  
  plot6 = plot6 + geom_path(data=get(contours[1]), aes(x=x, y=y, col = contours_w_aic[1]), linewidth=lwd) +
    geom_path(data=get(contours[2]), aes(x=x, y=y, col = contours_w_aic[2]), linewidth=lwd) +
    geom_path(data=get(contours[3]), aes(x=x, y=y, col = contours_w_aic[3]), linewidth=lwd) #+
  # geom_path(data=get(contours[4]), aes(x=x, y=y, col = contours_w_aic[4])) +
  # geom_path(data=get(contours[5]), aes(x=x, y=y, col = contours_w_aic[5])) +
  # geom_path(data=get(contours[6]), aes(x=x, y=y, col = contours_w_aic[6])) +
  # geom_path(data=get(contours[7]), aes(x=x, y=y, col = contours_w_aic[7])) +
  # geom_path(data=get(contours[8]), aes(x=x, y=y, col = contours_w_aic[8]))
  
  plot6 = plot6 + scale_color_manual(name="", breaks=contours_w_aic, values=pallette)
  
}

# combining probability plots ---------------------------------------------------------

combined2 = plot4 + plot5 + plot6 & theme(legend.position="bottom")
print(combined2 + plot_layout(guides="collect"))


# just contours -----------------------------------------------------------

plot7 = ggplot(cond_dens_df, aes(x=hs, y=stp)) + geom_raster(aes(fill=((V1)))) +
  scale_fill_gradientn(colours=c("white", "yellow", "orange", "red", "black"),
                       breaks=c(0, 5e-7, 1e-6), limits=c(0,1e-6)) +
  labs(fill="CDE")

plot7 = ggplot() +
  ylim(layer_scales(plot7)$y$range$range[1], layer_scales(plot7)$y$range$range[2]) +
  theme_linedraw() +
  scale_x_continuous(limits = c(0, 20), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0.001, 0.09), expand = c(0, 0)) +
  theme( panel.background = element_rect(fill = NA),
         panel.ontop = TRUE, legend.title = element_text(size=20), legend.text = element_text(size=15), legend.key.size = unit(1.5, 'cm'),
         axis.text=element_text(size=15), axis.title=element_text(size=20)) +
  xlab("Hs") + ylab("S2")

# data 
sample.df = data.frame(
  x = hs_peak,
  y = stp_peak
)

plot7 = plot7 + geom_point(data=sample.df, aes(x=x, y=y), col="grey") 
# read in contours 
cntr_path = "~/GitHub/env-contours/FORM fits/"

contours=c("neggev_exp_exp_form_p1.38888888888889e-05", "posweibull_lin_qua_form_p1.38888888888889e-05", "neglnorm_lin_lin_form_p1.38888888888889e-05" ,  "neggamma_qua_lin_form_p1.38888888888889e-05"   ,
           "negweibull_exp_exp_form_p1.38888888888889e-05" ,  "posgamma_lin_lin_form_p1.38888888888889e-05" ,    "poslnorm_qua_exp_form_p1.38888888888889e-05" ,  "posgev_exp_qua_form_p1.38888888888889e-05") 

contours_w_aic =c("C1", "C2" , "C3", "C4", "C5", "C6", "C7", "C8")

for (cntr in contours){
  cntr_file = load(paste(cntr_path, cntr, sep=""))
  assign(cntr, get(cntr_file))
  assign(cntr, data.frame(x= get(cntr)$x, y=get(cntr)$y))
}

# adding contours 
pallette = c("red", "orange", "yellow4", "lightgreen",
 "green", "lightblue", "blue", "purple",
"violet", "black")


plot7 = plot7 + geom_path(data=get(contours[1]), aes(x=x, y=y, col = contours_w_aic[1]), linewidth=lwd) +
  geom_path(data=get(contours[2]), aes(x=x, y=y, col = contours_w_aic[2]), linewidth=lwd) +
  geom_path(data=get(contours[3]), aes(x=x, y=y, col = contours_w_aic[3]), linewidth=lwd) + 
geom_path(data=get(contours[4]), aes(x=x, y=y, col = contours_w_aic[4]), linewidth=lwd) +
geom_path(data=get(contours[5]), aes(x=x, y=y, col = contours_w_aic[5]), linewidth=lwd) +
geom_path(data=get(contours[6]), aes(x=x, y=y, col = contours_w_aic[6]), linewidth=lwd) +
geom_path(data=get(contours[7]), aes(x=x, y=y, col = contours_w_aic[7]), linewidth=lwd) +
geom_path(data=get(contours[8]), aes(x=x, y=y, col = contours_w_aic[8]), linewidth=lwd)

plot7 = plot7 + scale_color_manual(name="", breaks=contours_w_aic,
                                   values=pallette) 
print(plot7)




# just data ---------------------------------------------------------------

  plot8 = ggplot(cond_dens_df, aes(x=hs, y=stp)) + geom_raster(aes(fill=((V1)))) +
  scale_fill_gradientn(colours=c("white", "yellow", "orange", "red", "black"),
                       breaks=c(0, 5e-7, 1e-6), limits=c(0,1e-6)) +
  labs(fill="CDE")

plot8 = ggplot() +
  ylim(layer_scales(plot8)$y$range$range[1], layer_scales(plot8)$y$range$range[2]) +
  theme_linedraw() +
  scale_x_continuous(limits = c(0, 13), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0.0125, 0.075), expand = c(0, 0)) +
  theme( panel.background = element_rect(fill = NA),
         panel.ontop = TRUE, legend.title = element_text(size=20), legend.text = element_text(size=15), legend.key.size = unit(1.5, 'cm'),
         axis.text=element_text(size=15), axis.title=element_text(size=20)) +
  xlab("Hs") + ylab("S2")

# data 
sample.df = data.frame(
  x = hs_peak,
  y = stp_peak
)

sample.df2 = data.frame(
  x = hs_peak,
  y = t2_peak
)

plot8.5 = plot8 + geom_point(data=sample.df2, aes(x=x, y=y), col="black") + scale_y_continuous(limits = c(3, 12.5), expand = c(0, 0)) + ylab("T2")
plot8 = plot8 + geom_point(data=sample.df, aes(x=x, y=y), col="black") 


print(plot8)
print(plot8.5)