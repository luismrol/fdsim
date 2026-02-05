### Simulation example
remotes::install_github("luismrol/fdsim")
library(fdsim)
d1<-sim_out(10, 200, 100, moderate = TRUE, perc_out= c(0.05, 0.0))
d2<-sim_out(10, 200, 100, moderate = FALSE, perc_out= c(0.05, 0.05))
d3<-sim_out(10, 200, 100, moderate = TRUE, perc_out= c(0.25, 0.25))


depth<-lapply(d1, function(x) {fdsim::uni_depth(x$data, depth = "maha")})


c_matrix(d1, method = boxplot_outliers, depth = "maha")
