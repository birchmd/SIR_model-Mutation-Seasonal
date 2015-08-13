## This program implements the SI model with mutation of the virus ##

source("../SI_mut_consts.R")
source("./SI_mut_fncs.R")

calcParams <- calcParamsList(alphaMin=1e-15, alphaMax=100, N=2000, Tmax=1, Tstep=0.05)
initList <- list(I_0=0.05, S_0=0.95, mean=3, var=1, initDist="norm")
solutionList <- odeSolve(calcParams, initList, paramList)

#--- ouput time evolution of (normalized) i(alpha,t) --- #
#plotFrames(solutionList,filedir="./img/",by=1)

with(solutionList, {
  plot(times[1:nrow(i)],S,type="l",ylim=c(0,max(c(S,I))),ylab=" ",xlab="t",col="red")
  lines(times[1:nrow(i)], I,col="blue")
  
  plot(times[1:nrow(i)],varalpha,type="l",xlab="t")
  plot(times[1:nrow(i)],alphabar,type="l",xlab="t")
})

dyn.unload(paste(sep="/", rootDir, "numericImplimentation", compDrvFile))
#save.image(file="SI_mut.RData")
