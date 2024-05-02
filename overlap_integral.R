library("pracma")

# get contours  -----------------------------------------------------------
fits =  c("neggev_exp_exp", "negweibull_exp_exp", "neglnorm_lin_lin", "neggamma_qua_lin",
          "posgev_exp_qua", "posweibull_lin_qua", "poslnorm_qua_exp", "posgamma_lin_lin")

file_names = paste("FORM fits/",sapply((sapply(fits, list.files, path="FORM fits", simplify=FALSE)), '[[', 1), sep="")

for (i in 1:length(fits)){
  load(file_names[i])
  if (exists("form")){
    assign(fits[i], form)
    rm(form)
  }
  if (exists("fit")){
    assign(fits[i], fit)
    rm(fit)
  }
}


# load cond_density -------------------------------------------------------

cond_dens = read.csv("cond_dens_pos.csv") # currently for structure B, can load in different file for others
dens = read.csv("env_probs.csv")
cond_dens = cond_dens/sum(cond_dens)

# calculate integral ------------------------------------------------------

overlap_int = function(cntr, cell_x, cell_y, dens, plot=F){
  if (plot){
    plot(cell_x, cell_y)
    mask = inpolygon(cell_x, cell_y, cntr$x, cntr$y) 
    points(cell_x[mask], cell_y[mask], col="yellow")
  }
  sum(inpolygon(cell_x, cell_y, cntr$x, cntr$y) * dens)
}

contour_data = sapply(fits, get)
overlap_ps = lapply(contour_data, overlap_int, cell_x=dens$x, cell_y=dens$y, dens=cond_dens)

metric = function(p){
  2 * p -1
}

lapply(overlap_ps, metric)
