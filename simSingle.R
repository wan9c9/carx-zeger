#------------------------------
simSingle <- function(iRep)
{
  #message("rep:",iRep)
  tryCatch({
  simData <- carxSim(nObs=nObs,prmtrAR=prmtrAR,prmtrX=prmtrX,sigmaEps=sigma,lcl=lcl,ucl=ucl,seed=37513*iRep)

  t0 <- proc.time()
  source("./carx_zeger.R")
  cz = carx_zeger(y=simData$y,
                x=as.matrix(simData[,c("X1","X2")]),
                ci=simData$ci,
                lcl=simData$lcl,
                ucl=simData$ucl,
                p=1,
                prmtrX=prmtrX,
                prmtrAR=prmtrAR,
                sigma=sigma,
                verbose=T
                )
  t1 <- proc.time()
  c0 <- carx(y~X1+X2-1,data=simData,p=1,prmtrX=prmtrX,prmtrAR=prmtrAR,sigmaEps=sigma,CI.compute=F)
  t2 <- proc.time()
  ret <- list(zeger=cz,wang=c0,tzeger=t1-t0,twang=t2-t1)
  f <- paste0(sprintf("rep_%05d.RData",iRep))
  save(ret,file=f)
  message(sprintf("Rep: %5d is done, result is saved in %s",iRep,f))
  
  invisible(ret)
  },warning=function(w){message(sprintf("Rep: %5d, warnings.",iRep))
  },error=function(w){message(sprintf("Rep: %5d, errors.",iRep))
  },finally={}
  )
}
