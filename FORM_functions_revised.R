#source("~/GitHub/env-contours/FORM_functions.R")

# parametric families ------------------------------------------------------

constant = function(theta, x){
  rep(theta[1], length(x))
}

linear = function(theta, x){
  a = theta[1] ; b = theta[2]
  a + b*x
}

exponential = function(theta, x){
  a = theta[1] ; b = theta[2] ; c = theta[3]
  a + b*exp(c*x)
}

quadratic = function(theta, x){
  a = theta[1] ; b = theta[2] ; c = theta[3]
  a*(x+b)^2 + c
}

par_func = function(type){
  switch(type,
         "constant"=constant,
         "quadratic"=quadratic,
         "linear"=linear,
         "exponential"=exponential)
}
