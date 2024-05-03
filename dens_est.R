library(cond.extremes)
source("FORM_functions.R")

## read in data
data <- read.csv("data.txt")

## isolate peaks
hs_peak <- c()
t2_peak <- c()
for (i in 1:max(data$Idn)){
  hs_peak[i] <- max(data$Hs[data$Idn == i])
  t2_peak[i] <- max(data$T2[data$Idn == i])
}

stp_peak <- (2 * pi * hs_peak) / (9.81 * t2_peak^2)
plot(hs_peak, stp_peak)


# generate dens -----------------------------------------------------------

dens <- grid.dens(data = matrix(c(hs_peak, stp_peak), ncol = 2),
                  q.marginal = 0.8, q.cond = 0.95, xlim = c(0, 30),
                  ylim = c(0.01, 0.08), nx = 180, ny = 90, log = T, adjust = 0.2)
print(sum(dens$p))

write.csv(dens, "env_probs.csv", row.names = FALSE)
