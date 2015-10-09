nObs <- 400
prmtrAR <- 0.5
prmtrX <- c(0.2,0.4)
sigma <- 0.6

cr <- 0
abscl <- 100


lcl <- -abscl
ucl <- abscl

ind_result_dir <- sprintf("./sim_result/n_%04d_cr_%02d/",nObs,as.integer(cr*100))
