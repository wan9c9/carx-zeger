source("./simulation_setting.R",echo=T)
system(paste("cp ./simulation_setting.R ",ind_result_dir)) #copy setting to dir
setwd(ind_result_dir)

library(matrixStats)
truePrmtr <- c(prmtrX,prmtrAR,sigma)

cat("\nnObs:", nObs)
cat("\nnRep:", nRep)
cat("\nlarge sample cr:",cr)
cat("\nabs cl:",abscl)
cat("\ntrue prmtr:\n")
print(truePrmtr)

simTable0 <- matrix(nrow=nRep,ncol=length(truePrmtr))
colnames(simTable0) <- c(paste0("X",1:length(prmtrX)),paste0("AR",1:length(prmtrAR)),"sigma")

simTable1 <- matrix(nrow=nRep,ncol=length(truePrmtr))
colnames(simTable1) <- colnames(simTable0)

timeUsed <- matrix(nrow=nRep,ncol=2)
colnames(timeUsed) <- c("wang","zeger")
censorRate<- numeric(nRep)
cenStrLen <- numeric(nObs)

for(i in 1:nRep)
{
  load(file=sprintf("rep_%05d.RData",i-1))
  #load(file=paste0(ind_result_dir,sprintf("rep_%05d.RData",i-1)))
  r <- ret
  prmtr <- c(r$zeger$prmtrX,r$zeger$prmtrAR,r$zeger$sigma)
  for(cs in r$zeger$cenStrIdx)
    cenStrLen[cs$nc] <- cenStrLen[cs$nc] + 1
  simTable1[i,] <- prmtr
  simTable0[i,] <- r$wang$prmtrEstd
  censorRate <- r$wang$censorRate
  timeUsed[i,1] <- r$twang[3]
  timeUsed[i,2] <- r$tzeger[3]
  i <- i+1
}

bias0 <- colMeans(simTable0) - truePrmtr
sd0 <- colSds(simTable0)
mse0 <- bias0^2 + sd0^2
summary0 <- rbind(t(bias0),t(sd0),t(mse0))
rownames(summary0) <- c("bias","sd","mse")

bias1 <- colMeans(simTable1) - truePrmtr
sd1 <- colSds(simTable1)
mse1 <- bias1^2 + sd1^2
summary1 <- rbind(t(bias1),t(sd1),t(mse1))
rownames(summary1) <- c("bias","sd","mse")

timeStat <- cbind(colMeans(timeUsed),colSds(timeUsed))
rownames(timeStat) <- c("Wang", "Zeger")
colnames(timeStat) <- c("mean", "sd")


write.table(simTable0,file="zegerResult.txt")
write.table(summary0,file="summary0.txt")

write.table(simTable1,file="wangResult.txt")
write.table(summary1,file="summary1.txt")

write.table(timeUsed,file="timeUsed.txt")
write.table(timeStat,file="timeStat.txt")

#```

print("The average censoring rate is ")
#```{r}
mcr <- mean(censorRate)
print(mcr)
write(mcr,"averageCensorRate.txt")
#```

print("The bias, s.d. and mse of our method:")
#```{r}
print(summary0)
write.table(summary0,"summary_wc.txt")
#```

print("Similar statistics by zeger:")
#```{r}
print(summary1)
#```

print("The MSE ratio (wang/zeger):")
#```{r}
mseRatio <- mse0/mse1
print(mseRatio)

save(mse0,file="mse0.RData")
write.table(mse0,file="mse0.txt")
save(mse1,file="mse1.RData")
write.table(mse1,file="mse1.txt")
save(mseRatio,file="mseRatio.RData")
write.table(mseRatio,file="mseRatio.txt")
#```

print("Averge time used by both methods:")
#```{r}
timeStat <- cbind(colMeans(timeUsed),colSds(timeUsed))
rownames(timeStat) <- c("Wang", "Zeger")
colnames(timeStat) <- c("mean", "sd")
print(timeStat)
save(timeStat,file="timeStat.RData")
#```

#```{r hist-cenStrLen}
idx <- nObs
for(i in nObs:1)
  if(!all(cenStrLen[i:nObs]==0))
  {
    idx <-i
    break
  }
svg("hist_cenStrLen.svg")
plot(cenStrLen[1:idx]/nRep/nObs,xlab="Length of censored string", ylab="Empirical proability")
save(cenStrLen,file="cenStrLen.RData")
dev.off()
#```
