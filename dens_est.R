setwd("~/GitHub/environment-modelling")

library(cond.extremes)
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
plot(hs_peak, stp_peak)


# generate dens -----------------------------------------------------------

dens <- grid.dens(data = matrix(c(hs_peak, stp_peak), ncol = 2),
                  q.marginal = 0.8, q.cond = 0.95, xlim = c(0, 30),
                  ylim = c(0.01, 0.08), nx = 180, ny = 90, log = T, adjust = 0.2)
print(sum(dens$p))

write.csv(dens, "data/env_probs.csv", row.names = FALSE)
write.csv(dens, "~/GitHub/force_calculation_and_wave_sim/scripts/response_sim/env_probs.csv",
          row.names = FALSE)
