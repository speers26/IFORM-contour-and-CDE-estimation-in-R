# -------------------------------------------------------------------------

setwd("~/GitHub/env-contours")
source("~/GitHub/env-contours/FORM_functions_revised.R")
source("~/GitHub/R_packages/cond.extremes/R/thresh_select.R")
source("~/GitHub/env-contours/predictive_llh_cross_validation.R")
library("stats")
library("cond.extremes")

# overall function --------------------------------------------------------
form_fit_gamma = function(data, q, theta0, plot=F, p=1e-5, x_all){
  
  # fitting hs --------------------------------------------------------------
  hsfit = gpd.fit(hs, quantile(hs, q), show=F)
  
  # fitting stp model -------------------------------------------------------
  s2fit = optim(theta0, dep_gamma_nll, data=data, x_all=x_all)
  
  # making form -------------------------------------------------------------
  beta = qnorm(1-p)
  u_ctr = (beta*exp(2i * pi * (1:10000)/10000))
  
  FORM_x = sapply(pnorm(Re(u_ctr)), FUN=qspliced, x=hs, q=q, gpd_par=hsfit$mle)
  uh_matrix = matrix(data=c(Im(u_ctr), FORM_x), ncol=2)
  FORM_y = apply(uh_matrix, 1, gamma_inv_rsblt, shape_theta = s2fit$par[1:3],
                 rate_theta = s2fit$par[4:6])
  
  # plotting ---------------------------------------------------------------
  if (plot){
    plot(data[,1], data[,2], xlim=c(1, 16), cex=0.5, pch=16) ; lines(FORM_x, FORM_y, col="red")
  }
  
  return(list(hsfit=hsfit, s2fit=s2fit, x=FORM_x, y=FORM_y))
}


# fitting parametric form functions ---------------------------------------

dep_gamma_nll = function(data, all_theta, x_all){
  
  x = data[,1] ; y = data[,2]
  
  shape = shape_func(all_theta[1:3], x)
  rate = rate_func(all_theta[4:6], x)
  
  shape_all = shape_func(all_theta[1:3], x_all)
  rate_all = rate_func(all_theta[4:6], x_all)
  
  if(sum(shape_all<=0)>0){
    #print("1")
    return(1e10)
  }
  if(sum(rate_all<=0)>0){
    #print("2")
    return(1e10)
  }
  
  return(-sum(dgamma(y, shape=shape, rate=rate, log=TRUE)))
  
}

gamma_inv_rsblt = function(uh, shape_theta, rate_theta){
  
  u = uh[1] ; h = uh[2]
  
  shape = shape_func(shape_theta, h)
  rate = rate_func(rate_theta, h)
  
  p = pnorm(u)
  
  qgamma(p, shape=shape, rate=rate)
  
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

# parametric forms for gamma parameters -------------------------------------
shape_type="quadratic"
rate_type="linear"
shape_abr = substr(shape_type, 1, 3)
rate_abr = substr(rate_type, 1, 3)
shape_func = par_func(shape_type)
rate_func = par_func(rate_type)

if(neg){
  shape_theta0 = switch(shape_type,
                        "linear" = c(25, -2, 0),
                        "quadratic" = c(0.01, -12, 25),
                        "exponential" = c(20, 10, -1)
                        )
  rate_theta0 = switch(rate_type,
                       "linear" = c(1000, 1, 0), # c(500, 100, 0)
                       "quadratic" = c(0.01, 0, 500),
                       "exponential" = c(800, 200, -1)
                       )
}

if (1-neg){
  shape_theta0 = switch(shape_type,
                        "linear" = c(25, -2, 0),
                        "quadratic" = c(0.1, 0, 75),
                        "exponential" = c(30, 20, 0.4)
  )
  rate_theta0 = switch(rate_type,
                       "linear" = c(1000, 1, 0),
                       "quadratic" = c(0.1, 0, 2000),
                       "exponential" = c(800, 200, 0.4)
  )
}

# make FORM ---------------------------------------------------------------
p = 1e-3/72
theta0 = c(shape_theta0, rate_theta0)
form = form_fit_gamma(data, 0.8, theta0, plot=T, p=p, x_all=hs)

if(neg){
  form$y = -(form$y - max_stp - 0.001)
}
if (neg){
  save(form, file=paste("~/GitHub/env-contours/FORM fits/neggamma_",
                        shape_abr, "_", rate_abr, "_form_p", p, sep=""))
}
if (1-neg){
  save(form, file=paste("~/GitHub/env-contours/FORM fits/posgamma_",
                        shape_abr, "_", rate_abr, "_form_p", p, sep=""))
}

# 
# AIC ---------------------------------------------------------------------
AIC = (2 * form$s2fit$value + 2 * length(form$s2fit$par))

if (neg){
  save(AIC, file=paste("~/GitHub/env-contours/FORM_AICs/neggamma_", shape_abr, "_", rate_abr, "_AIC", sep=""))
}
if (1-neg){
  save(AIC, file=paste("~/GitHub/env-contours/FORM_AICs/posgamma_", shape_abr, "_", rate_abr, "_AIC", sep=""))
}


# checking stp parameter form ---------------------------------------------
gamma_NLL = function(x, theta){
  shape = theta[1] ; rate = theta[2]
  if (shape<=0|rate<=0){return(1e10)}
  else{
    -sum(dgamma(x, shape, rate, log=TRUE))
  }
}

gamma_fit = function(x, theta0){
  optim(theta0, gamma_NLL, x=x)
}

banded_pars_gamma = function(data, k){
  
  data[,2] = data[,2][order(data[,1])] ; data[,1] = sort(data[,1])
  bands = get_bands(data, k)
  fit = lapply(bands$bandsY, gamma_fit, theta0=c(1,1))
  mles = lapply(fit, '[[', 1)
  
  return(list(x = sapply(bands$bandsX, median),
              shape = as.numeric(lapply(mles, '[[', 1)),
              rate = as.numeric(lapply(mles, '[[', 2)))
  )
}

k=100
banded_pars = banded_pars_gamma(data, k)

par(mfrow=c(1,2))
plot(banded_pars$x, banded_pars$shape)
x = seq(0, 7, length=100)
lines(x, shape_func(shape_theta0, x), col="red")
lines(x, shape_func(form$s2fit$par[1:3], x), col="blue")

plot(banded_pars$x, banded_pars$rate)
x = seq(0, 7, length=100)
lines(x, rate_func(rate_theta0, x), col="red")
lines(x, rate_func(form$s2fit$par[4:6], x), col="blue")
par(mfrow=c(1,1))


# tail cross validation ---------------------------------------------------
set.seed(1)
if(1){
  R = 30 ; q = 0.9
  k_set = c(5, 10) ; q_cv_set = c(0.8, 0.9, 0)
  
  for (q_cv in q_cv_set){
    print(q_cv)
    for (k in k_set){
      leave = F
      if (k==10 && q_cv ==0.8){
        leave = T
      }
      print(k)
      llh_mean_set = c()
      for (r in 1:R){
        print(r)
        folds = split_training_test(get_permutation(data), k, threshold_q=q_cv)
        llh = c()
        for (i in 1:k){
          fold_fit = form_fit_gamma(folds[[i]]$train, q, theta0, x_all=hs)
          llh[i] = -dep_gamma_nll(folds[[i]]$test, fold_fit$s2fit$par, x_all=hs)
        }
        llh_mean_set[r] = (mean(llh))
        print(llh_mean_set)
      }
      if (neg){
        save(llh_mean_set, file=paste("~/GitHub/env-contours/FORM_cross_validation_scores/vary_par_forms/neggammaq",
                                      shape_abr, rate_abr, q_cv,"k",k, sep='_'))
      }
      if (1-neg){
        save(llh_mean_set, file=paste("~/GitHub/env-contours/FORM_cross_validation_scores/vary_par_forms/posgammaq",
                                      shape_abr, rate_abr, q_cv,"k",k, sep='_'))
      }
    }
  }
}

