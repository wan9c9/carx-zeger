source("simulation_setting.R")

generateCL <- F
if(generateCL){
  library(carx)
	nObs <- 100000
	simData <- carxSim(nObs=nObs,prmtrAR=prmtrAR,prmtrX=prmtrX,sigmaEps=sigma,lcl=lcl,ucl=ucl,seed=0)

	cr <- c(0.01,0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8)
	ucl <- quantile(simData$y,1-cr/2)
	lcl <- quantile(simData$y,cr/2)
	abscl <- round((ucl+abs(lcl))/2,2)
	save(abscl,file="abscl.RData")
	save(abscl,file="abscl.txt")
}

source("./simulation_setting.R")


n <- 100
crs <- c(1,5,10,20,30,40,50)

rslt <- matrix(nrow=length(crs),ncol=1+4)
timeStats <- matrix(nrow=length(crs),ncol=1+4)
i<- 1
for( i in 1:length(crs))
{
  ind_result_dir <- sprintf("./sim_result/n_%04d_cr_%02d/",nObs,as.integer(crs[i]))
  load(paste0(ind_result_dir,"mseRatio.RData"))
  rslt[i,-1] <- mseRatio
  rslt[i,1] <- crs[i]/100
  if(i == 1) colnames(rslt) <- c("cr",names(mseRatio))

  load(paste0(ind_result_dir,"timeStat.RData"))
  timeStats[i,] <- c(crs[i]/100,c(timeStat))
  if(i == 1) colnames(timeStats) <- c("cr","m_wc","m_zb","sd_wc","sd_zb")
  i <- i + 1
}

setEPS()
postscript(paste0("msePlot_n",nObs,".eps"))

ylim <- range(rslt[,-1])
plot(rslt[,"cr"],rslt[,"X1"],ylim=ylim,xlab="Censoring rate", ylab="MSE ratio (WC/ZB)",type="l",lty=1)
grid()
lines(rslt[,"cr"],rslt[,"X2"],lty=2)
lines(rslt[,"cr"],rslt[,"AR1"],lty=3)
lines(rslt[,"cr"],rslt[,"sigma"],lty=4)
legend("topleft",
			 c(expression(beta[1]), 
				 expression(beta[2]),
				 expression(psi[1]), 
				 expression(sigma)
				 ) ,
			 lty=c(1,2,3,4))

dev.off()

setEPS()
postscript(paste0("timeStats_n",nObs,".eps"))

ylim <- range(timeStats[,-1])
plot(timeStats[,"cr"],timeStats[,"m_zb"],log="y",type="l",ylim=ylim,lty=1,
     xlab="Censoring rate", ylab="Mean time per estimation")
grid()
lines(timeStats[,"cr"],timeStats[,"m_wc"],lty=2)  
#lines(rslt[,"cr"],rslt[,"AR1"],ylim=c(0,1),lty=3)
#lines(rslt[,"cr"],rslt[,"sigma"],ylim=c(0,1),lty=4)
legend("topleft",
			 c("ZB","WC"),
			 lty=c(1,2)
			 )

dev.off()


