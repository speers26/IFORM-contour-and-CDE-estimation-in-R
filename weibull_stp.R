# -------------------------------------------------------------------------

setwd("~/GitHub/env-contours")
source("~/GitHub/env-contours/FORM_functions_revised.R")
source("~/GitHub/env-contours/predictive_llh_cross_validation.R")
source("~/GitHub/R_packages/cond.extremes/R/thresh_select.R")
library("stats")
library("cond.extremes")

# overall function --------------------------------------------------------
form_fit_weibull = function(data, q, theta0, plot=F, p=1e-5, x_all){
  
  # fitting hs --------------------------------------------------------------
  hsfit = gpd.fit(hs, quantile(hs, q), show=F)
  
  # fitting stp model -------------------------------------------------------
  s2fit = optim(c(shape_theta0, scale_theta0), dep_weibull_nll, data=data, x_all=x_all)
  
  # making form -------------------------------------------------------------
  beta = qnorm(1-p)
  u_ctr = (beta*exp(2i * pi * (1:10000)/10000))
  
  FORM_x = sapply(pnorm(Re(u_ctr)), FUN=qspliced, x=data[,1], q=q, gpd_par=hsfit$mle)
  uh_matrix = matrix(data=c(Im(u_ctr), FORM_x), ncol=2)
  FORM_y = apply(uh_matrix, 1, weibull_inv_rsblt, shape_theta = s2fit$par[1:3],
                 scale_theta = s2fit$par[4:6]) 
  
  # plotting ----------------------------------------------------------------
  if (plot){
    plot(data[,1], data[,2], cex=0.5, pch=16)
    lines(FORM_x, FORM_y, col="red")
  }
  return(list(hsfit=hsfit, s2fit=s2fit, x=FORM_x, y=FORM_y))
}

# fitting stp modelfunctions -------------------------------------------------

dep_weibull_nll = function(data, all_theta, x_all){
  
  x = data[,1] ; y = data[,2]
  
  shape = shape_func(all_theta[1:3], x)
  scale = scale_func(all_theta[4:6], x)
  
  shape_all = shape_func(all_theta[1:3], x_all)
  scale_all = scale_func(all_theta[4:6], x_all)
  
  if(sum(shape_all<=0)>0){
    return(1e10)
  }
  if(sum(scale_all<=0)>0){
    return(1e10)
  }
  
  return(-sum(dweibull(y, shape=shape, scale=scale, log=TRUE)))
  
}

weibull_inv_rsblt = function(uh, shape_theta, scale_theta){
  
  u = uh[1] ; h = uh[2]
  
  shape = shape_func(shape_theta, h)
  scale = scale_func(scale_theta, h)
  
  p = pnorm(u)
  
  qweibull(p, shape=shape, scale=scale)
  
}


# read in data ------------------------------------------------------------

cnsTS = read.csv("data/cnsTS.txt")
data = cnsTS
hs <- c()
t2 <- c()
for (i in 1:max(data$StrIdn)){
  hs[i] <- max(data$Hs[data$StrIdn == i])
  t2[i] <- max(data$T2[data$StrIdn == i])
}

neg = T
stp = (hs*2*pi)/(t2^2*9.81)
max_stp = max(stp)
neg_stp = -stp + max_stp + 0.001

if (neg){
  data = matrix(data=c(hs, neg_stp), ncol=2)
}
if(1-neg){
  data = matrix(data=c(hs, stp), ncol=2)
}
data[,2] = data[,2][order(data[,1])] ; data[,1] = sort(data[,1])

# parametric forms for weibull parameters -------------------------------------
shape_type = "exponential"
scale_type = "exponential"
shape_func = par_func(shape_type)
scale_func = par_func(scale_type)
shape_abr = substr(shape_type, 1, 3)
scale_abr = substr(scale_type, 1, 3)
if(1-neg){
  shape_theta0 = switch(shape_type,
                      "linear"=c(5, 2.5, 0),
                      "quadratic"=c(0.1, 0, 9),
                      "exponential"=c(10, 1, 0.25))
  scale_theta0 = switch(scale_type,
                       "linear"=c(0.2075, 0.05, 0),
                       "quadratic"=c(-0.01/20, -15, 2.285/20),
                       "exponential"=c(0.03069107, 0.01289684, 0.10101747))
}
if(neg){
  shape_theta0 = switch(shape_type,
                      "exponential"=c(4,2, -1),
                      "linear"=c(6,-0.5, 0),
                      "quadratic"=c(0.005031297, -17.040816241 ,  2.607981848))#c(0.05,-14, 0))
  scale_theta0 = switch(scale_type,
                       "exponential"=c(0.0125, 0.025, -0.25),
                       "quadratic"=c(0.0001, -16, 0),
                       "linear"=c(0.031031587 ,-0.002124255, -7.572580238))# c(0.0375, -0.00075, 0))
}

# fitting form ------------------------------------------------------------
theta0 = c(scale_theta0, shape_theta0)
p = 1e-3/72
form = form_fit_weibull(data, 0.8, theta0, plot=T, p=p, x_all=hs) 
if(neg){
  form$y = -(form$y - max_stp - 0.001)
}
if (neg){
  save(form, file=paste("~/GitHub/env-contours/FORM fits/negweibull_",
                        shape_abr, "_", scale_abr, "_form_p", p, sep=""))
}
if (1-neg){
  save(form, file=paste("~/GitHub/env-contours/FORM fits/posweibull_",
                        shape_abr, "_", scale_abr, "_form_p", p, sep=""))
}

# 
# AIC ---------------------------------------------------------------------
AIC = (2 * form$s2fit$value + 2 * length(form$s2fit$par))

if (neg){
  save(AIC, file=paste("~/GitHub/env-contours/FORM_AICs/negweibull_", shape_abr, "_", scale_abr, "_AIC", sep=""))
}
if (1-neg){
  save(AIC, file=paste("~/GitHub/env-contours/FORM_AICs/posweibull_", shape_abr, "_", scale_abr, "_AIC", sep=""))
}


# checking stp parameter form ---------------------------------------------
weibull_NLL = function(x, theta){
  shape = theta[1] ; scale = theta[2]
  if (shape<=0|scale<=0){return(1e10)}
  else{
    -sum(dweibull(x, shape, scale, log=TRUE))
  }
}

weibull_fit = function(x, theta0){
  optim(theta0, weibull_NLL, x=x)
}

banded_pars_weibull = function(data, k){
  data[,2] = data[,2][order(data[,1])] ; data[,1] = sort(data[,1])
  bands = get_bands(data, k)
  fit = lapply(bands$bandsY, weibull_fit, theta0=c(1,1))
  mles = lapply(fit, '[[', 1)
  
  return(list(x = sapply(bands$bandsX, median),
              shape = as.numeric(lapply(mles, '[[', 1)),
              scale = as.numeric(lapply(mles, '[[', 2)))
  )
}

k=100
banded_pars = banded_pars_weibull(data, k)

par(mfrow=c(1,2))
x = seq(0, 7, length=100)

plot(banded_pars$x, banded_pars$shape)
lines(x, shape_func(shape_theta0, x), col="red")
lines(x, shape_func(form$s2fit$par[1:3], x), col="blue")

plot(banded_pars$x, banded_pars$scale, xlim=c(0, 7))
lines(x, scale_func(scale_theta0, x), col="red")
lines(x, scale_func(form$s2fit$par[4:6], x), col="blue")

par(mfrow=c(1,1))

# tail cross validation ---------------------------------------------------
set.seed(1)
if(1){
  R = 30 ; q = 0.9
  k_set = c(5, 10) ; q_cv_set = c(0.8, 0.9, 0)
  for (q_cv in q_cv_set){
    print(q_cv)
    for (k in k_set){
      # leave = F
      # if (k==10 && q_cv ==0.8){
      #   leave = T
      # }
      # if (leave){break}
      print(k)
      llh_mean_set = c()
      for (r in 1:R){
        print(r)
        folds = split_training_test(get_permutation(data), k, threshold_q=q_cv)
        llh = c()
        for (i in 1:k){
          fold_fit = form_fit_weibull(folds[[i]]$train, q, theta0, x_all=hs)
          llh[i] = -dep_weibull_nll(folds[[i]]$test, fold_fit$s2fit$par, x_all=hs)
        }
        llh_mean_set[r] = (mean(llh))
        print(llh_mean_set)
      }
      
      if (neg){
        save(llh_mean_set, file=paste("~/GitHub/env-contours/FORM_cross_validation_scores/vary_par_forms/negweibullq",
                                      shape_abr, scale_abr, q_cv,"k",k, sep='_'))
      }
      if (1-neg){
        save(llh_mean_set, file=paste("~/GitHub/env-contours/FORM_cross_validation_scores/vary_par_forms/posweibullq",
                                      shape_abr, scale_abr, q_cv,"k",k, sep='_'))
      }
    }
  }
}

