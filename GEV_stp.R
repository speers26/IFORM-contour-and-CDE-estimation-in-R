# -------------------------------------------------------------------------
source("FORM_functions_revised.R")
source("predictive_llh_cross_validation.R")
library("stats")
library("cond.extremes")

# overall function --------------------------------------------------------
form_fit_GEV = function(data, q, theta0, plot=F, p=1e-5){
  
  # fitting hs --------------------------------------------------------------
  hsfit = gpd.fit(hs, quantile(hs, q), show=F)
  
  # fitting stp model -------------------------------------------------------
  xi = gev.fit(data[,2], show=FALSE)$mle[3]
  s2fit = optim(theta0, GEV_nll, data=data, xi=xi, range_y=range(data[,2]))
  # NLL = 1e10
  # while(NLL>=1e10){
  #   s2fit = optim(theta0, GEV_nll, data=data, xi=xi, range_y=range(data[,2]))
  #   NLL = s2fit$value
  #   xi = xi + 0.1
  # }

  # making form -------------------------------------------------------------
  beta = qnorm(1-p)
  u_ctr = (beta*exp(2i * pi * (1:10000)/10000))
  FORM_x = sapply(pnorm(Re(u_ctr)), FUN=qspliced, x=data[,1], q=q, gpd_par=hsfit$mle)
  uh_matrix = matrix(data=c(Im(u_ctr), FORM_x), ncol=2)
  FORM_y = apply(uh_matrix, 1, GEV_inv_rsblt, mu_theta = s2fit$par[1:3],
                 sig_theta = s2fit$par[4:6], xi=xi)

  # plotting ----------------------------------------------------------------
  if (plot){
    plot(data[,1], data[,2], cex=0.5, pch=16)
    lines(FORM_x, FORM_y, col="red")
  }
  
  return(list(hsfit=hsfit, s2fit=s2fit, xi=xi, x=FORM_x, y=FORM_y))
}

# fitting stp model functions  --------------------------------------------

GEV_nll = function(data, all_theta, xi, range_y){
  
  x = data[,1] ; y = data[,2]
  
  mu = mu_func(all_theta[1:3], x)
  sig = sig_func(all_theta[4:6], x)
  
  if(sum(sig<=0)>0){
    #print("1")
    return(1e10)
  }
  if(xi>0 & sum(range_y[1]<(mu - sig/xi))>0){
    #print("2")
    return(1e10)
  }
  if(xi<0 & sum(range_y[2]>(mu-sig/xi))>0){
    #print("3") 
    return(1e10)
  }
  
  return(-sum(dgev(y, location=mu, scale=sig, shape=xi, log=TRUE)))
  
}
GEV_inv_rsblt = function(uh, mu_theta, sig_theta, xi){
  
  u = uh[1] ; h = uh[2]
  
  mu = mu_func(mu_theta, h) 
  sigma = sig_func(sig_theta, h)
  
  p = pnorm(u)
  
  qgev(p, location=mu, scale=sigma, shape=xi)
  
}

# read in data ------------------------------------------------------------

cnsTS = read.csv("cnsTS.txt")
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

mu_theta_0 = c(0.055, -0.04, -0.5)

if (neg){
  data = matrix(data=c(hs, neg_stp), ncol=2)
}
if(1-neg){
  data = matrix(data=c(hs, stp), ncol=2)
}

# parametric forms for GEV parameters -------------------------------------
mu_type = "exponential"
sig_type = "exponential"

mu_abr = substr(mu_type, 1, 3)
sig_abr = substr(sig_type, 1, 3)

mu_func = par_func(mu_type)
sig_func = par_func(sig_type)

if(1-neg){
  mu_theta_0 = switch(mu_type,
                      "exponential"=c(0.055, -0.04, -0.5),
                      "linear"=c(0.03, 0.004, 0),
                      "quadratic"=c(-0.00015, -16, 0.05+0.02))
  sig_theta_0 = switch(sig_type,
                       "exponential"=c(0.004, 0.02, -0.75),
                       "quadratic"=c(1.840716e-05+0.005, -1.412231e+01,  3.919991e-03),
                       "linear" = c(0.006, -0.001/6, 0))
}
if(neg){
  mu_theta_0 = switch(mu_type,
                      "exponential"= c(-0.007616481,  0.036772154, -0.092734397),
                      "linear"=c(0.027049802, -0.002239907, 5.418543004),
                      "quadratic"=c(0.01/100, -16, 0))
  sig_theta_0 = switch(sig_type,
                       "exponential"= c(0.00357652,  0.00313629, -0.57849285),
                       "quadratic"=c(1.833977e-05, -3.340960e+00,  4.084472e-03),
                       "linear"=c(0.003, 0, 0))
}

# fitting form ------------------------------------------------------------

theta0 = c(mu_theta_0, sig_theta_0)
p = 1e-3/72
form = form_fit_GEV(data, 0.8, theta0, plot=T, p=p)
if(neg){
  form$y = -(form$y - max_stp - 0.001)
}
if (neg){
  save(form, file=paste("neggev_",
                        mu_abr, "_", sig_abr, "_form_p", p, sep=""))
}
if (1-neg){
  save(form, file=paste("posgev_",
                        mu_abr, "_", sig_abr, "_form_p", p, sep=""))
}

# 
# AIC ---------------------------------------------------------------------
AIC = (2 * form$s2fit$value + 2 * length(form$s2fit$par) + 2)

if (neg){
  save(AIC, file=paste("neggev_", mu_abr, "_", sig_abr, "_AIC", sep=""))
}
if (1-neg){
  save(AIC, file=paste("posgev_", mu_abr, "_", sig_abr, "_AIC", sep=""))
}

# checking stp parameter form ---------------------------------------------

banded_pars_GEV = function(data, k){
  
  data[,2] = data[,2][order(data[,1])] ; data[,1] = sort(data[,1])
  bands = get_bands(data, k)
  fit = lapply(bands$bandsY, gev.fit, show=F)
  mles = lapply(fit, '[[', 7)
  
  return(list(x = sapply(bands$bandsX, median),
              mu = as.numeric(lapply(mles, '[[', 1)),
              sigma = as.numeric(lapply(mles, '[[', 2)),
              xi = as.numeric(lapply(mles, '[[', 3)))
         )
}

banded_pars = banded_pars_GEV(data, 100)

par(mfrow=c(1, 3))
plot(banded_pars$x, banded_pars$xi, xlab="Hs", ylab="xi")
lines(banded_pars$x, rep(form$xi, length(banded_pars$x)), col="blue")

plot(banded_pars$x, banded_pars$mu, xlab="Hs", ylab="mu")
lines(banded_pars$x, mu_func(form$s2fit$par[1:3], banded_pars$x), col="red")
lines(banded_pars$x, mu_func(mu_theta_0,banded_pars$x), col="blue")

plot(banded_pars$x, banded_pars$sigma, xlab="Hs", ylab="sigma")
lines(banded_pars$x, sig_func(form$s2fit$par[4:6], banded_pars$x), col="red")
lines(banded_pars$x, sig_func(sig_theta_0, banded_pars$x), col="blue")

par(mfrow=c(1,1))

# tail cross validation ---------------------------------------------------
set.seed(1)
if(1){
  R = 30 ; q = 0.9
  k_set = c(5,10) ; q_cv_set = c(0.8, 0.9, 0)
  
  for (q_cv in q_cv_set){
    print(q_cv)
    for (k in k_set){
      print(k)
      # leave = F
      # if (k==10 && q_cv ==0.8){
      #   leave = T
      # }
      llh_mean_set = c()
      for (r in 1:R){
        print(r)
        folds = split_training_test(get_permutation(data), k, threshold_q=q_cv)
        llh = c()
        for (i in 1:k){
          fold_fit = form_fit_GEV(folds[[i]]$train, q, theta0)
          llh[i] = -GEV_nll(folds[[i]]$test, fold_fit$s2fit$par, xi=fold_fit$xi, range_y = range(stp))
        }
        llh_mean_set[r] = (mean(llh))
        print(llh_mean_set)
      }
      if (neg){
        save(llh_mean_set, file=paste("neggevq",
                                      mu_abr, sig_abr, q_cv,"k",k, sep='_'))
      }
      if (1-neg){
        save(llh_mean_set, file=paste("posgevq",
                                      mu_abr, sig_abr,q_cv,"k",k, sep='_'))
      }
    }
  }
}

# AIC ---------------------------------------------------------------------
print(2 * form$s2fit$value + 2 * length(form$s2fit$par) + 2)
print(-form$s2fit$value)



