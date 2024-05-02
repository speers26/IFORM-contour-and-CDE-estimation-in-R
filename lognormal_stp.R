### R script to draw FORM/IFORM bivariate environmental contours ###
setwd("~/GitHub/env-contours")
source("~/GitHub/env-contours/FORM_functions_revised.R")
source("~/GitHub/R_packages/cond.extremes/R/thresh_select.R")
source("~/GitHub/env-contours/predictive_llh_cross_validation.R")
library("stats")
library("cond.extremes")

# overall function --------------------------------------------------------
form_fit_lnorm = function(data, q, theta0, plot=F, p=1e-5, x_all){
  
  # fit hs
  v = quantile(data[,1], q, names=F)
  hsfit = gpd.fit(data[,1], q, show=FALSE)
  
  # fit s2
  s2fit = optim(theta0, lgnrm_NLL, data=data, x_all=x_all)
  
  # making FORM
  beta = qnorm(1-p)
  u_ctr = (beta*exp(2i * pi * (1:10000)/10000))
  
  FORM_x = sapply(pnorm(Re(u_ctr)), FUN=qspliced, x=hs, q=q, gpd_par=hsfit$mle)
  uh_matrix = matrix(data=c(Im(u_ctr), FORM_x), ncol=2)
  FORM_y = apply(uh_matrix, 1, lnorm_inv_rsn, theta=s2fit$par)
  
  # plotting
  if (plot){
    plot(data[,1], data[,2], cex=0.5, pch=16) ; lines(FORM_x, FORM_y, col="red")
  }
  
  return(list(hsfit=hsfit, s2fit=s2fit, x=FORM_x, y=FORM_y))
  
}

# fitting parametric form functions ---------------------------------------


lgnrm_NLL = function(data, theta, x_all){
  
  x = data[,1] ; y = data[,2]
  mu_theta = theta[1:3] ; sd_theta = theta[4:6]
  
  mu_val = mu(mu_theta, x)
  sd_val = sigma(sd_theta, x)
  
  sd_val_all = sigma(sd_theta, x_all)
  
  if(sum(sd_val_all<=0)>0){
    #print("hi")
    return(1e10)
  }
  
  return(-sum(dlnorm(y, meanlog=mu_val, sdlog=sd_val, log=TRUE)))
  
}

lnorm_inv_rsn = function(uh, theta, stp=T){
  
  u = uh[1] ; h = uh[2]
  # extract parameters from theta
  mu = mu(theta[1:3], h) ; sigma = sigma(theta[4:6], h)
  
  # get normal probability
  p = pnorm(u)
  
  # inverse cdf
  x = qlnorm(p, meanlog=mu, sdlog=sigma)
  
  return(x)
  
}

# read in data ------------------------------------------------------------
cnsTS = read.csv("data/cnsTS.txt")
data = cnsTS
hs <- c()
tp <- c()
for (i in 1:max(data$StrIdn)){
  hs[i] <- max(data$Hs[data$StrIdn == i])
  tp[i] <- max(data$T2[data$StrIdn == i])
}

stp = (hs*2*pi)/(tp^2*9.81)
max_stp = max(stp)
neg_stp = -stp + max_stp + 0.001

neg = T
if(neg){
  data = matrix(c(hs, neg_stp), ncol=2)
}
if(1-neg){
  data = matrix(c(hs, stp), ncol=2)
}

# parametric forms --------------------------------------------------------
mu_type = "linear"
sig_type = "linear"
mu = par_func(mu_type)
sigma = par_func(sig_type)
mu_abr = substr(mu_type, 1, 3)
sig_abr = substr(sig_type, 1, 3)
if(1-neg){
  mu_theta0 = switch(mu_type,
                     "exponential" = c(-2.9, -2.8, -1),
                     "linear" = c(-3.35, 0.08, 0),
                     "quadratic" = c(-.01, -12, -2.3) # c(-.01, -8, -2.8)
                     )
  sig_theta0 = switch(sig_type,
                     "exponential" = c(0.075, 1, -1),
                     "quadratic" = c(0.001, -10, 0.2),
                     "linear"= c(0.20, -0.015, 0)
  )
}
if (neg){
  mu_theta0 = switch(mu_type,
                     "exponential" = c(1.6, -1, 0.1),
                     "linear" = c(-3.4, -0.1, 0), # c(-3.45, -0.13, 0)
                     "quadratic" = c(0.01, 15, -8)
  )
  sig_theta0 = switch(sig_type,
                      "exponential" = c(0.2534281, -0.2465539, -1.0508090),
                      "quadratic" =  c(0.005, 0.1, 0.15), 
                      "linear" =  c(0.1, 0.05, 0)
  )
}

theta0 = c(mu_theta0, sig_theta0)

# fitting form ------------------------------------------------------------
p = 1e-3/72
fit = form_fit_lnorm(data, 0.8, theta0, plot=T, p=p, x_all=hs)

if(neg){
  fit$y = -(fit$y - max_stp - 0.001)
}
if (neg){
  save(fit, file=paste("~/GitHub/env-contours/FORM fits/neglnorm_",
                        mu_abr, "_", sig_abr, "_form_p", p, sep=""))
}
if (1-neg){
  save(fit, file=paste("~/GitHub/env-contours/FORM fits/poslnorm_",
                        mu_abr, "_", sig_abr, "_form_p", p, sep=""))
}

# 
# AIC ---------------------------------------------------------------------
AIC = (2 * fit$s2fit$value + 2 * length(fit$s2fit$par))

if (neg){
  save(AIC, file=paste("~/GitHub/env-contours/FORM_AICs/neglnorm_", mu_abr, "_", sig_abr, "_AIC", sep=""))
}
if (1-neg){
  save(AIC, file=paste("~/GitHub/env-contours/FORM_AICs/poslnorm_", mu_abr, "_", sig_abr, "_AIC", sep=""))
}


# check parametric form fit k--------------------------------------------
banded_pars = function(data, k){
  
  data[,2] = data[,2][order(data[,1])] ; data[,1] = sort(data[,1])
  bands = get_bands(data, k)
  log_bands = lapply(bands$bandsY, log)
  mus = lapply(log_bands, mean)
  sds = lapply(log_bands, sd)
  
  return(list(x = sapply(bands$bandsX, median),
              mu = mus,
              sd = sds)
  )
}

k=100
emp_pars = banded_pars(data,k)

par(mfrow=c(1,2))
plot(emp_pars$x, emp_pars$mu, xlim=c(2, 20))
lines(emp_pars$x, mu(fit$s2fit$par[1:3], emp_pars$x), col="red")
lines(emp_pars$x, mu(theta0[1:3], emp_pars$x), col="blue")

plot(emp_pars$x, emp_pars$sd)
lines(emp_pars$x, sigma(fit$s2fit$par[4:6], emp_pars$x), col="red")
lines(emp_pars$x, sigma(theta0[4:6], emp_pars$x), col="blue")
par(mfrow=c(1,1))

# tail cross validation ---------------------------------------------------
set.seed(1)
if(1){
  R = 30 ; q = 0.9
  k_set = c(5, 10) ; q_cv_set = c(0, 0.8, 0.9)
  
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
          fold_fit = form_fit_lnorm(folds[[i]]$train, q, theta0, x_all=hs)
          llh[i] = -lgnrm_NLL(folds[[i]]$test, fold_fit$s2fit$par, x_all=hs)
        }
        llh_mean_set[r] = (mean(llh))
        print(llh_mean_set)
      }
      if (neg){
        save(llh_mean_set, file=paste("~/GitHub/env-contours/FORM_cross_validation_scores/vary_par_forms/neglnormq",
                                      mu_abr, sig_abr, q_cv,"k",k, sep='_'))
      }
      if (1-neg){
        save(llh_mean_set, file=paste("~/GitHub/env-contours/FORM_cross_validation_scores/vary_par_forms/poslnormq",
                                      mu_abr, sig_abr, q_cv,"k",k, sep='_'))
      }
    }
  }
}




