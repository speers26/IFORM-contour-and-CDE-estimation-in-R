# IFORM-contour-and-CDE-estimation-in-R
R Code used alongside code in https://github.com/speers26/Linear-Waves-and-Base-Shear-Estimation- and https://github.com/speers26/Conditional-Extremes-in-R- to produce results in https://arxiv.org/abs/2404.16775

To reproduce results found in https://arxiv.org/abs/2404.16775

1: Install R package on conditional extremes https://github.com/speers26/Conditional-Extremes-in-R-
2: Estimate extreme region density for example data using dens_est.R
3: Use Python package https://github.com/speers26/Linear-Waves-and-Base-Shear-Estimation- to estimate the conditional density of the environment
4: Use lognormal_stp.R, weibull_stp.R, gamma_stp.R and GEV_stp.R to fit IFORM contours discussed in the paper
5: Use dens_plots.R to generate plots seen in paper

