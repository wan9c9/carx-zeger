library(tmvtnorm)
library(matrixcalc)
library(emulator)
library(parallel)
library(matrixStats)
backShift <- function(ts,bsCoef,idx=ifelse(is.matrix(ts),dim(ts)[1],length(ts)))
{
  if(is.matrix(ts))
  {
    as.vector(t(ts[idx:(idx-length(bsCoef)+1),])%*%bsCoef)
  }
  else
  {
    if(is.vector(ts))
      as.numeric(ts[idx:(idx-length(bsCoef)+1)]%*%bsCoef)
  }
}
#test
#backShift(c(1,2),c(2,1))
#v = backShift( matrix(c(1,2,2,1),nrow=2),c(2,1))


covAR2 <- function(prmtrAR, sigma, lags)
{
	roots <- polyroot(c(1,-prmtrAR))
	#print(roots)
	if(any(abs(roots)<=1))
		warning(" prmtrAR is not stationary",paste(prmtrAR,sep=',')," Roots:", roots)
	val <- stats::ARMAacf(ar=prmtrAR, lag.max=max(lags))
	#print(val)
	val <- as.vector(val)
  nlag <- length(lags)

	mat <- matrix(nrow=nlag,ncol=nlag)
	for(i in 1:nlag)
	{
		mat[i,i] <- 1
		if(i > 1){
			for(j in 1:(i-1))
			{
				mat[i,j] <- val[abs(lags[i]-lags[j])+1]
				mat[j,i] <- mat[i,j]
			}
		}
	}

	v <- sigma^2/(1-prmtrAR%*%val[2:(length(prmtrAR)+1)])
	mat <- v[1,1]*mat
	return(mat)
}


#' this function implements the exact MLE of Zeger and Brookmeyer (1986) method for estimating the parameters of regression models with censored autoregressive data.

carx_zeger <- function(y,x,ci=NULL,lcl=NULL,ucl=NULL,
                       p=1,prmtrX=NULL,prmtrAR=NULL,sigma=NULL,
                       tol=1e-4,max.iter=100,
                       verbose=T
                       )
{
  if(p>1)
    stop("function is not suitable for p>1")
  #find the U index s.t. for each t in U, y_t,...y_{t-p} are not censored.
  U <- NULL
  nObs <- length(y)
  if(is.vector(x))
    nX <- 1
  else
    nX <- dim(x)[2]

  for(i in (p+1):nObs)
    if(all(ci[i:(i-p)]==0))
      U <- c(U,i)
  if(is.null(U))
    message("no consecutive p+1 uncensored observations are found.")

  message(paste0("censoring rate: ", mean(1.0*abs(ci))))
  #find indices for censored strings (Y)
  cenStrIdx = list()
  ncs = NULL
  cntU = 0
  inCS = FALSE
  newCenStr = FALSE
  for(i in 1:nObs)
  {
    if(ci[i] != 0) #censored
    {
      if(!inCS)
      {
        inCS = TRUE
        ncs = list(beg = i)
      }
      cntU = 0
    }
    else
    {
      if(inCS)
      {
        cntU = cntU + 1
        if(cntU == p)
        {
          ncs$end = i
          ncs$u = which(ci[ncs$beg:ncs$end]==0)+ncs$beg-1
          ncs$c = which(ci[ncs$beg:ncs$end]!=0)+ncs$beg-1
          ncs$nu = length(ncs$u)
          ncs$nc = length(ncs$c)
          z <- NULL
          if(ncs$beg>1)
            z <- max(1,ncs$beg-p):(ncs$beg-1)
          ncs$z <- z
          ncs$nz <- length(z)
          cenStrIdx[[length(cenStrIdx)+1]] = ncs
          cntU = 0
          inCS = FALSE
        }
      }
    }
  }
  if(inCS) #the last cenStr is still open
  {
    ncs$end = nObs
    ncs$u = which(ci[ncs$beg:ncs$end]==0)+ncs$beg-1
    ncs$c = which(ci[ncs$beg:ncs$end]!=0)+ncs$beg-1
    z <- NULL
    if(ncs$beg>1)
      z <- max(1,ncs$beg-p):(ncs$beg-1)
    ncs$z <- z
    ncs$nz <- length(z)
    ncs$nu = length(ncs$u)
    ncs$nc = length(ncs$c)
    cenStrIdx[[length(cenStrIdx)+1]] = ncs
  }

  if(is.null(prmtrX))  prmtrX <- numeric(nX)
  if(is.null(prmtrAR))  prmtrAR <- numeric(p)
  if(is.null(sigma)) sigma <- 1.0

	trend <- numeric(nObs)
	eps <- numeric(nObs)

  updateTrend <- function()
  {
		trend <<- as.vector(x%*%prmtrX)
    eps <<- y - trend
	}

  updatePrmtrX <- function()
  {
    t_theta <- c(1,-prmtrAR)
    Sxx <- matrix(0,nrow=nX,ncol=nX)
    Sxy <- numeric(nX)
    for(i in U)
    {
      tmp <- backShift(x,t_theta,i)
      Sxx <- Sxx + tmp%o%tmp/sigma^2
      Sxy <- Sxy + tmp*backShift(y,t_theta,i)/sigma^2
    }

    for(csidx in cenStrIdx)
    {
      if(length(csidx$u)>0) #uncensored obs
      {
        tmpIdx <- c(csidx$z,csidx$u)
        varMat <- covAR2(prmtrAR,sigma,tmpIdx-tmpIdx[1]+1)
        meanVec <- trend[tmpIdx]
        if(is.null(csidx$z))
        {
          #message("z is null")
          #browser()
          Xs_ju <- x[csidx$u,,drop=FALSE]
          Ys_ju <- y[csidx$u]
          Sigma_ju <- varMat
          iSigma_ju <- solve(Sigma_ju)
          Sxx <- Sxx + t(Xs_ju)%*%iSigma_ju%*%Xs_ju
          Sxy <- Sxy + t(Xs_ju)%*%iSigma_ju%*%Ys_ju
        }
        else
        {
          nz <- length(csidx$z)
          cdist <- conditionalDistMvnorm(y[csidx$z], 1:nz,
                                         meanVec, varMat)
          Sigma_ju <- cdist$var
          iSigma_ju <- solve(Sigma_ju)
          M_ju <- varMat[-(1:nz),1:nz]%*%solve(varMat[(1:nz),(1:nz)])
          Xs_ju <- x[csidx$u,,drop=FALSE] - M_ju%*%x[csidx$z,,drop=FALSE]
          Ys_ju <- y[csidx$u] - M_ju%*%y[csidx$z]
          Sxx <- Sxx + t(Xs_ju)%*%iSigma_ju%*%Xs_ju
          Sxy <- Sxy + t(Xs_ju)%*%iSigma_ju%*%Ys_ju
        }
      }

      #censored obs
      tmpIdx <- c(csidx$z,c(csidx$beg:csidx$end))
      varMat <- covAR2(prmtrAR,sigma,tmpIdx-tmpIdx[1]+1)
      meanVec <- trend[tmpIdx]
      tmpIdxU <- c(csidx$z, csidx$u)
      cdist <- conditionalDistMvnorm(y[tmpIdxU],
                                     tmpIdxU-tmpIdx[1]+1,
                                     meanVec,
                                     varMat)
      Sigma_jc <- cdist$var
      iSigma_jc <- solve(Sigma_jc)

      M_jc <- varMat[csidx$c-tmpIdx[1]+1,tmpIdxU-tmpIdx[1]+1]%*%
              solve(varMat[tmpIdxU-tmpIdx[1]+1,tmpIdxU-tmpIdx[1]+1])
      Xs_jc <- x[csidx$c,] - M_jc%*%x[tmpIdxU,]
      ##compute Yh_jc
      leftCensored <- ci[csidx$c] < 0 #data are left censored
      rightCensored <- ci[csidx$c] > 0 #right censored
      lower <- rep(-Inf,length(csidx$c))
      lower[rightCensored] <- ucl[csidx$c][rightCensored]
      upper <- rep(Inf, length(csidx$c))
      upper[leftCensored] <- lcl[csidx$c][leftCensored]
      mt <- mtmvnorm(cdist$mean,cdist$var,lower,upper,doComputeVariance=FALSE)
      Yh_jc <- mt$tmean
      Ys_jc <- Yh_jc - M_jc%*%y[tmpIdxU]
      Sxx <- Sxx + t(Xs_jc)%*%iSigma_jc%*%Xs_jc
      Sxy <- Sxy + t(Xs_jc)%*%iSigma_jc%*%Ys_jc
    }
    prmtrX <<- solve(Sxx,Sxy)
    updateTrend()
  }

  updateSigma <- function()
  {
    rss <- 0
    for(i in U)
      rss <- rss + (eps[i] - prmtrAR%*%eps[(i-1):(i-p)])^2
    en <- length(U)

    for(csidx in cenStrIdx)
    {

      #uncensored obs
      if(csidx$nu>0)
      {
        tmpIdx <- c(csidx$z,csidx$u)
        varMat <- covAR2(prmtrAR,sigma,tmpIdx-tmpIdx[1]+1)
        meanVec <- trend[tmpIdx]
        if(!is.null(csidx$z))
        {
          cdist <- conditionalDistMvnorm(y[csidx$z], 1:length(csidx$z), meanVec, varMat)
          Sigma_ju <- cdist$var
          delta <- y[csidx$u] - cdist$mean
        }
        else
        { #no preceding Z series
          delta <- eps[csidx$u]
          Sigma_ju <- varMat
        }
        rss <- rss + quad.form.inv(Sigma_ju,delta)*sigma^2
      }

      #censored obs
      tmpIdx <- c(csidx$z,c(csidx$beg:csidx$end))
      varMat <- covAR2(prmtrAR,sigma,tmpIdx-tmpIdx[1]+1)
      meanVec <- trend[tmpIdx]
      tmpIdxU <- c(csidx$z, csidx$u) #all uncensored obs needed
      #if(length(tmpIdxU)==0) #this shouldn't happen
      #  cdist <- list(mean=meanVec,var=varMat)
      #else
        cdist <- conditionalDistMvnorm(y[tmpIdxU],
                                     tmpIdxU-tmpIdx[1]+1,
                                     meanVec,
                                     varMat)
      Sigma_jc <- cdist$var
      iSigma_jc <- solve(Sigma_jc)

      ##compute Yh_jc and V_jc
      leftCensored <- ci[csidx$c] < 0 #data are left censored
      rightCensored <- ci[csidx$c] > 0 #right censored
      lower <- rep(-Inf,length(csidx$c))
      lower[rightCensored] <- ucl[csidx$c][rightCensored]
      upper <- rep(Inf, length(csidx$c))
      upper[leftCensored] <- lcl[csidx$c][leftCensored]
      mt <- mtmvnorm(cdist$mean,cdist$var,lower,upper,doComputeVariance=TRUE)
      delta <- mt$tmean - cdist$mean
      rss <- rss + quad.form(iSigma_jc,delta)*sigma^2
      en <- en  + csidx$nc + csidx$nu - matrix.trace(mt$tvar%*%iSigma_jc)
    }
    sigma <<- as.numeric(sqrt(rss/en))
  }

  logLklhd <- function(phi)
  {
    gamma0 <- sigma^2/(1-phi^2)

    ret <- 0
    for(i in U)
      ret <- ret + dnorm(eps[i],mean=phi*eps[i-1],sd=sigma,log=T)
	  
 
    if(ci[1]!=0)
    {
      #browser()
      csidx <- cenStrIdx[[1]]
      v_jc <- csidx$nc

      #uncensored obs
      ret <- ret + dnorm(eps[csidx$end],mean=0,sd=sqrt(gamma0),log=T)
      #censored obs
      M1 <- matrix(nrow=v_jc,ncol=v_jc)
      for(i in 1:v_jc)
        for(j in 1:v_jc)
          M1[i,j] <- phi^abs(i-j)

      M2 <- matrix(nrow=v_jc,ncol=1)
      for(i in 1:v_jc)
        M2[i,1] <- phi^(v_jc+1-i)
      M <- M1-M2%*%t(M2)
      Sigma_jc <- gamma0*M
      
      tmpEps <- eps[csidx$end]
      eta_jc <- as.vector(trend[csidx$c]+tmpEps*M2)
      
      leftCensored <- ci[csidx$c] < 0 #data are left censored
      rightCensored <- ci[csidx$c] > 0 #right censored
      lower <- rep(-Inf,v_jc)
      lower[rightCensored] <- ucl[csidx$c][rightCensored]
      upper <- rep(Inf, length(csidx$c))
      upper[leftCensored] <- lcl[csidx$c][leftCensored]
      prob <- pmvnorm(lower,upper,mean=as.vector(eta_jc),sigma=Sigma_jc)
      ret <- ret + log(prob)
    }

    #ignore <- T
		#if(!ignore){
		#	message("keeping censored string")
	   
    for(csidx in cenStrIdx)
    {
      #message("computing score ar1: ",csidx$beg,"--",csidx$end)
      if(csidx$beg==1 | (csidx$end==nObs & csidx$nu==0)) #special case at the beginning or end
        next
      #part 1: uncensored obs
      v_jc <- csidx$nc
      #nu==1 here
      nt <- csidx$nc + csidx$nu #number of total elements in the cenStr

      eta_ju <- trend[csidx$u] + phi^nt*eps[csidx$z]
      Sigma_ju <- (1-phi^(2*nt))*gamma0

      ret <- ret + dnorm(y[csidx$u],mean=eta_ju,sd=sqrt(Sigma_ju),log=T)

      #part2 censored obs
      ##compute Sigma_jc


      M1 <- matrix(nrow=v_jc,ncol=v_jc)
      #M1 is in the form of
      #[ 1, phi,...,phi^{v_jc-1}]
      #[ phi,1,,...,phi^{v_jc-2}]
      #...
      for(i in 1:v_jc)
        for(j in 1:v_jc)
          M1[i,j] <- phi^abs(i-j)

      M2 <- matrix(nrow=v_jc,ncol=2)
      for(i in 1:v_jc)
        M2[i,] <- c(phi^i, phi^(nt-i))
      M3 <- matrix(nrow=2,ncol=2)
      M3[1,1] <- M3[2,2] <- 1
      M3[1,2] <- M3[2,1] <- -phi^nt #M3 is the inverse of covariance matrix
      M3 <- M3/(1-phi^(2*nt))
      M <- M1-M2%*%M3%*%t(M2)
      Sigma_jc <- gamma0*M
      tmpEps <- eps[c(csidx$z,csidx$u)]
      eta_jc <- trend[csidx$c] +M2%*%M3%*%tmpEps #M3 is already the inverse

      ##compute Yh_jc and V_jc
      leftCensored <- ci[csidx$c] < 0 #data are left censored
      rightCensored <- ci[csidx$c] > 0 #right censored
      lower <- rep(-Inf,v_jc)
      lower[rightCensored] <- ucl[csidx$c][rightCensored]
      upper <- rep(Inf, length(csidx$c))
      upper[leftCensored] <- lcl[csidx$c][leftCensored]
      prob <- pmvnorm(lower,upper,mean=as.vector(eta_jc),sigma=Sigma_jc)
      ret <- ret + log(prob)
    }
    
    #check the last cenStr
    if(length(cenStrIdx)>0)
    {
      csidx <- cenStrIdx[[length(cenStrIdx)]]
      if( csidx$nu==0 )
      {
        #last one is censored, i.e., all values in cenStr are censored
         #the follow computes the contribution of the censored obs
        v_jc <- csidx$nc
        M1 <- matrix(nrow=v_jc,ncol=v_jc)
        for(i in 1:v_jc)
          for(j in 1:v_jc)
            M1[i,j] <- phi^abs(i-j)

        M2 <- matrix(nrow=v_jc,ncol=1)
        for(i in 1:v_jc)
          M2[i,1] <- phi^i
        M <- M1-M2%*%t(M2)
        Sigma_jc <- gamma0*M

        tmpEps <- eps[csidx$z]
        eta_jc <- as.vector(trend[csidx$c]+tmpEps*M2)

        ##compute Yh_jc and V_jc
        leftCensored <- ci[csidx$c] < 0 #data are left censored
        rightCensored <- ci[csidx$c] > 0 #right censored
        lower <- rep(-Inf,v_jc)
        lower[rightCensored] <- ucl[csidx$c][rightCensored]
        upper <- rep(Inf, length(csidx$c))
        upper[leftCensored] <- lcl[csidx$c][leftCensored]
        prob <- pmvnorm(lower,upper,mean=eta_jc,sigma=Sigma_jc)
        ret <- ret + log(prob)
      }
    }
		#}
		ret
  }

  scoreAR1 <- function(phi)
  {
    gamma0 <- sigma^2/(1-phi^2)
    dgamma0 <- 2*phi/(1-phi^2)*gamma0
    score <- 0
    for(i in U)
      score <- score + (eps[i]-phi*eps[i-1])*eps[i-1]/sigma^2
    #check if the first obs is censored
    if(ci[1]!=0) #the first obs is censored
    {
      csidx <- cenStrIdx[[1]]
      v_jc <- csidx$nc

      #the first obs is censored
      #implying the first cenStrIdx begins at the beginning of the sequence
      score <- score - phi/(1-phi^2) + phi/sigma^2*(eps[csidx$end])^2 #contribution of Y_j^u

      #the follow computes the contribution of the censored obs
      M1 <- matrix(nrow=v_jc,ncol=v_jc)
      dM1 <- matrix(nrow=v_jc,ncol=v_jc)
      for(i in 1:v_jc)
        for(j in 1:i)
        {
          M1[i,j] <- phi^abs(i-j)
          if(i == j)
            dM1[i,j] <- 0
          else
            dM1[i,j] <- abs(i-j)*phi^(abs(i-j)-1)

          if(i>j)
          {
            M1[j,i] <- M1[i,j]
            dM1[j,i] <- dM1[i,j]
          }
        }

      M2 <- matrix(nrow=v_jc,ncol=1)
      dM2 <- matrix(nrow=v_jc,ncol=1)
      for(i in 1:v_jc)
      {
        M2[i,1] <- phi^(v_jc+1-i)
        dM2[i,1] <- (v_jc+1-i)*phi^(v_jc-i)
      }
      M <- M1-M2%*%t(M2)
      dM <- dM1 - dM2%*%t(M2) - M2%*%t(dM2)
      Sigma_jc <- gamma0*M
      dSigma_jc <- dgamma0*M +gamma0*dM
      iSigma_jc <- solve(Sigma_jc)
      diSigma_jc <- -iSigma_jc%*%dSigma_jc%*%iSigma_jc

      tmpEps <- eps[csidx$end]
      eta_jc <- as.vector(trend[csidx$c]+tmpEps*M2)
      deta_jc <- tmpEps*dM2

      ##compute Yh_jc and V_jc
      leftCensored <- ci[csidx$c] < 0 #data are left censored
      rightCensored <- ci[csidx$c] > 0 #right censored
      lower <- rep(-Inf,v_jc)
      lower[rightCensored] <- ucl[csidx$c][rightCensored]
      upper <- rep(Inf, length(csidx$c))
      upper[leftCensored] <- lcl[csidx$c][leftCensored]
      mt <- mtmvnorm(as.vector(eta_jc),Sigma_jc,lower,upper,doComputeVariance=TRUE)
      delta <- mt$tmean - eta_jc
      score <- score - 0.5*matrix.trace(iSigma_jc%*%dSigma_jc)
                     - 0.5*quad.form(diSigma_jc,delta)
                     + t(deta_jc)%*%iSigma_jc%*%delta
                     - 0.5*matrix.trace(diSigma_jc%*%mt$tvar)

    }
    for(csidx in cenStrIdx)
    {
      #message("computing score ar1: ",csidx$beg,"--",csidx$end)
      if(csidx$beg==1 | (csidx$end==nObs & csidx$nu==0)) #special case at the beginning or end
        next
      #part 1: uncensored obs
      v_jc <- csidx$nc
      #nu==1 here
      nt <- csidx$nc + csidx$nu #number of total elements in the cenStr

      eta_ju <- trend[csidx$u] + phi^nt*eps[csidx$z]
      deta_ju <- nt*phi^(nt-1)*eps[csidx$z]
      Sigma_ju <- (1-phi^(2*nt))*gamma0
			dSigma_ju <- -2*(nt*phi^(2*nt-1)*(1-phi^2)-phi*(1-phi^(2*nt)))/(1-phi^2)^2*sigma^2
      diSigma_ju <- -dSigma_ju/Sigma_ju^2
      #print(diSigma_ju)
      #diSigma_ju <- -2*(phi*(1-phi^(2*nt))-(1-phi^2)*nt*phi^(2*nt-1))/(1-phi^(2*nt))^2/sigma^2
      #print(diSigma_ju)
      #stop("")

      delta <- y[csidx$u] - eta_ju
      score <- score - 0.5*dSigma_ju/Sigma_ju  #ZB paper, 3.4b, part1
      score <- score - 0.5*delta^2*diSigma_ju #., 3.4b part2
      score <- score + deta_ju/Sigma_ju*delta #., 3.4b part3

      #part2 censored obs
      ##compute Sigma_jc
      M1 <- matrix(nrow=v_jc,ncol=v_jc)
      dM1 <- matrix(nrow=v_jc,ncol=v_jc)
      for(i in 1:v_jc)
        for(j in 1:v_jc)
        {
          M1[i,j] <- phi^abs(i-j)
          if(i == j)
            dM1[i,j] <- 0
          else
            dM1[i,j] <- abs(i-j)*phi^(abs(i-j)-1)
        }

      M2 <- matrix(nrow=v_jc,ncol=2)
      dM2 <- matrix(nrow=v_jc,ncol=2)
      for(i in 1:v_jc)
      {
        M2[i,] <- c(phi^i, phi^(nt-i))
        dM2[i,] <- c(i*phi^(i-1), (nt-i)*phi^(nt-i-1))
      }
      M3 <- matrix(nrow=2,ncol=2)
      dM3 <- matrix(nrow=2,ncol=2)
      M3[1,1] <- M3[2,2] <- 1
      M3[1,2] <- M3[2,1] <- -phi^nt #M3 is the inverse of covariance matrix
			#checkhere
      M3 <- M3/(1-phi^(2*nt))
      dM3 <- 2*nt*phi^(2*nt-1)/(1-phi^(2*nt))*M3
      dM3[1,2] <- dM3[1,2] - nt*phi^(nt-1)/(1-phi^(2*nt))
      dM3[2,1] <- dM3[1,2]
      M <- M1-M2%*%M3%*%t(M2)
      dM <- dM1 - dM2%*%M3%*%t(M2) - M2%*%dM3%*%t(M2) - M2%*%M3%*%t(dM2)
      Sigma_jc <- gamma0*M
      dSigma_jc <- dgamma0*M +gamma0*dM
      iSigma_jc <- solve(Sigma_jc)
      diSigma_jc <- -iSigma_jc%*%dSigma_jc%*%iSigma_jc

      tmpEps <- eps[c(csidx$z,csidx$u)]
      eta_jc <- trend[csidx$c] +M2%*%M3%*%tmpEps
      deta_jc <- (dM2%*%M3 + M2%*%dM3)%*%tmpEps

      ##compute Yh_jc and V_jc
      leftCensored <- ci[csidx$c] < 0 #data are left censored
      rightCensored <- ci[csidx$c] > 0 #right censored
      lower <- rep(-Inf,v_jc)
      lower[rightCensored] <- ucl[csidx$c][rightCensored]
      upper <- rep(Inf, length(csidx$c))
      upper[leftCensored] <- lcl[csidx$c][leftCensored]
      mt <- mtmvnorm(as.vector(eta_jc),Sigma_jc,lower,upper,doComputeVariance=TRUE)
      delta <- mt$tmean - eta_jc
      score <- score - 0.5*matrix.trace(iSigma_jc%*%dSigma_jc)
                     - 0.5*quad.form(diSigma_jc,delta)
                     + t(deta_jc)%*%iSigma_jc%*%delta
                     - 0.5*matrix.trace(diSigma_jc%*%mt$tvar)
    }

    #check the last cenStr
    if(length(cenStrIdx)>0)
    {
      csidx <- cenStrIdx[[length(cenStrIdx)]]
      if( csidx$nu==0 )
      {
        #last one is censored, i.e., all values in cenStr are censored
         #the follow computes the contribution of the censored obs
        v_jc <- csidx$nc
        M1 <- matrix(nrow=v_jc,ncol=v_jc)
        dM1 <- matrix(nrow=v_jc,ncol=v_jc)
        for(i in 1:v_jc)
          for(j in 1:v_jc)
          {
            M1[i,j] <- phi^abs(i-j)
            if(i == j)
              dM1[i,j] <- 0
            else
              dM1[i,j] <- abs(i-j)*phi^(abs(i-j)-1)
          }

        M2 <- matrix(nrow=v_jc,ncol=1)
        dM2 <- matrix(nrow=v_jc,ncol=1)
        for(i in 1:v_jc)
        {
          M2[i,1] <- phi^i
          dM2[i,1] <- i*phi^(i-1)
        }
        M <- M1-M2%*%t(M2)
        dM <- dM1 - dM2%*%t(M2) - M2%*%t(dM2)
        Sigma_jc <- gamma0*M
        dSigma_jc <- dgamma0*M +gamma0*dM
        iSigma_jc <- solve(Sigma_jc)
        diSigma_jc <- -iSigma_jc%*%dSigma_jc%*%iSigma_jc

        tmpEps <- eps[csidx$z]
        eta_jc <- as.vector(trend[csidx$c]+tmpEps*M2)
        deta_jc <- dM2*tmpEps

        ##compute Yh_jc and V_jc
        leftCensored <- ci[csidx$c] < 0 #data are left censored
        rightCensored <- ci[csidx$c] > 0 #right censored
        lower <- rep(-Inf,v_jc)
        lower[rightCensored] <- ucl[csidx$c][rightCensored]
        upper <- rep(Inf, length(csidx$c))
        upper[leftCensored] <- lcl[csidx$c][leftCensored]
        mt <- mtmvnorm(eta_jc,Sigma_jc,lower,upper,doComputeVariance=TRUE)
        delta <- mt$tmean - eta_jc
        score <- score - 0.5*matrix.trace(iSigma_jc%*%dSigma_jc)
                       - 0.5*quad.form(diSigma_jc,delta)
                       + t(deta_jc)%*%iSigma_jc%*%delta
                       - 0.5*matrix.trace(diSigma_jc%*%mt$tvar)
      }
    }
	#message("score: ",score, " @prmtrAR: ",phi)
   score
  }


  updatePrmtrAR <- function()
  {
    if(p == 1)
    {
      currAR <- prmtrAR
      #updatePrmtrAR1(); message("updated ar by score:", prmtrAR)
      ret <- optim(currAR,logLklhd,method="Brent",lower=-0.99,upper=0.99,control=list(fnscale=-1))
      #message("updated ar by llk  :", ret$par)
      prmtrAR <<- ret$par

      debug <- F
      if(debug){
        idx <- seq(0.4,0.95,0.01)
        disScore <- numeric(length(idx))
        llk <- numeric(length(idx))
        for(i in 1:length(idx)) 
        {
          disScore[i] <- scoreAR1(idx[i])
          llk[i] <- logLklhd(idx[i])
        }
        par(new=F)
        par(mar=c(5,4,4,5)+.1)
        plot(idx,disScore,type="l",ylim = range(c(disScore,llk)))
        grid()
        par(new=T)
        plot(idx,llk,type="l",col="blue",xaxt="n",yaxt="n")
        grid()
        axis(4,col="blue",col.axis="blue")
        #readline("Press <return to continue\n")
      }
    }
    else
      stop("function not implemented for AR order  > 1")
  }

  updatePrmtrAR1 <- function()
  {
    #specialization code for AR1
    message("solving score equation root for phi")
    nr <- uniroot(scoreAR1,c(-0.99,0.99) )
    message("solving score equation root for phi, done.")
    prmtrAR <<- nr$root
    #message("updated AR prmtr:", prmtrAR)
  }

  estimate <- function()
  {
    delta <- 1
    iter <- 1
    while(delta>tol && iter < max.iter)
    {
      prevPrmtr <- c(prmtrX,prmtrAR,sigma)
      #updateTrend()
      updatePrmtrX() #updateTrend()
      updatePrmtrAR()
      updateSigma()
      currPrmtr <- c(prmtrX,prmtrAR,sigma)
      delta <- sqrt(sum((currPrmtr-prevPrmtr)^2))/sqrt(sum(prevPrmtr^2))
      if(verbose)
        message("Iter:",iter,", delta: ", delta, " prmtrX:", paste0(prmtrX,del=", "),"prmtrAR:",prmtrAR,", sigma:",sigma)
      iter <- iter + 1
    }
  }

  estimate()

  list(U = U,
       cenStrIdx = cenStrIdx,
       prmtrX = prmtrX,
       prmtrAR = prmtrAR,
       sigma = sigma
       )
}
