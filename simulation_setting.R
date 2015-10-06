#In this file the settings of the simulation studies are written.
#It will be sourced by 
#   simulation_prepare.R, 
#   simulation_single.R,
#   simulation_summarize.Rmd
nRep <- 1000
nObs <- 400
prmtrAR <- 0.5
prmtrX <- c(0.2,0.4)
sigma <- 0.6
cr_idx <- 1
iRep <- 715

cr_all <- c(0.01,
        0.05,
        0.1,
        0.2,
        0.3,
        0.4,
        0.5,
        0.6,
        0.7,
        0.8
				)

abscl_all <-c(
			2.13, #0.01
			1.62, #0.05
			1.35, #0.1
			1.06, #0.2
			0.85, #0.3
			0.69, #0.4
			0.56, #0.5
			0.43, #0.4
			0.32, #0.3
			0.21 #0.2
			)

batch_variable_nObs <- Sys.getenv("carx_batch_variable_nObs")
batch_variable_cr_idx <- Sys.getenv("carx_batch_variable_cr_idx")

if(batch_variable_nObs!="")
	nObs <- as.integer(batch_variable_nObs)
if(batch_variable_cr_idx!="")
	cr_idx <- as.integer(batch_variable_cr_idx)
	

cr <- cr_all[cr_idx]
abscl <- abscl_all[cr_idx]

lcl <- -abscl
ucl <- abscl
ind_result_dir <- sprintf("./sim_result/n_%04d_cr_%02d/",nObs,as.integer(cr*100))
#ind_result_dir <- sprintf("./sim_result/n_%04d_abscl_%4.2f/",nObs,abscl)
