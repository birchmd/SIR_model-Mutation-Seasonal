i_dist <- function(alpha, initList, paramList) { #initial distribution of infectious individuals

  if(initList$initDist == "lognorm"){
    sigma <- sqrt(log((initList$var + initList$mean^2)/(initList$mean^2)))
    mu <- log(initList$mean) - 0.5*initList$var^2
    result <- dlnorm(alpha,mu,sigma) #lognormal
  }else if(initList$initDist == "norm"){
    mu <- initList$mean
    sigma <- sqrt(initList$var)
    result <- dnorm(alpha,mu,sigma) #normal
  }
  
  return(result * initList$I_0 / trapz(alpha, result))
}

# --- solve the PDE as set of coupled ODEs --- #
odeSolve <- function(calcParams, initList, paramList, returnDist){
  alpha <- seq(calcParams$alphaMin,calcParams$alphaMax,len=calcParams$N) #virulence range to use
  times <- seq(0,calcParams$Tmax,by=calcParams$Tstep) #times vector for ODE solver
  
  i <- i_dist(alpha, initList, paramList) #inital infectious individual density
  
  y <- c(c(initList$S_0),i) #state vector for ODE solver
  
  if(!(is.loaded("derivs", PACKAGE="compiledDerivatives"))){
    dyn.load(paste(sep="/", rootDir, "numericImplimentation", compDrvFile))
  }
  
  parms <- c(paramList$nu, paramList$mu, paramList$betaMult,
             paramList$gamma, calcParams$alphaMin,calcParams$alphaMax, paramList$D,
             paramList$seasAmp, paramList$seasPer)
  
  solution <- ode(y,times,"derivs",parms=parms,
                  dllname="compiledDerivatives",
                  maxsteps=1000,method="lsoda") #solve system of ODEs
  S <- solution[,2]
  i <- as.matrix(solution[,seq(3,calcParams$N+2)])
  if(min(i) < 0) i <- 0.5*(abs(i)+i) #eliminate negative infectious densities
  
  I <- rep(0,nrow(i))
  for (j in 1:nrow(i)){ #integrate i(alpha,t) to get I(t)
    I[j] = trapz(alpha, as.vector(i[j,]))
  }
  
  alphabar <- rep(0,nrow(i))
  varalpha <- rep(0,nrow(i))
  skewalpha <- rep(0,nrow(i))
  alphafourth <- rep(0,nrow(i))
  for (j in 1:nrow(i)){ #compute the mean, variance and skewness of the virulence
    alphabar[j] <- trapz(alpha, alpha*as.vector(i[j,]))/I[j]
    varalpha[j] <- trapz(alpha, (alpha - alphabar[j])^2*as.vector(i[j,]))/I[j]
    skewalpha[j] <- trapz(alpha, ((alpha - alphabar[j])/varalpha[j])^3*as.vector(i[j,]))/I[j]
    alphafourth[j] <- trapz(alpha, (alpha - alphabar[j])^4*as.vector(i[j,]))/I[j]
  }
  
  return(list(times=times, alpha=alpha, S=S, I=I, i=i, alphabar=alphabar, varalpha=varalpha, 
              skewalpha=skewalpha, alphafourth=alphafourth, title="ode"))
}

momDiffs1 <- function(t,y,parms, ...){
  S <- y[1]
  I <- y[2]
  alphabar <- y[3]
  
  A <- S*I*(beta(alphabar) + 0.5*parms$sigma*betaPrimePrime(alphabar))*(1 + parms$seasAmp*sin(2*pi*(t/parms$seasPer)))
  dS <- parms$nu - A - parms$mu*S
  dI <- A - alphabar*I - parms$mu*I
  dalphabar <- parms$sigma*(S*betaPrime(alphabar) - 1) + parms$D*parms$i0
  
  return(list(c(dS=dS, dI=dI, dalphabar=dalphabar)))
}

momDiffs2 <- function(t,y,parms, ...){
  S <- y[1]
  I <- y[2]
  alphabar <- y[3]
  varalpha <- y[4]
  
  A <- S*I*(beta(alphabar) + 0.5*varalpha*betaPrimePrime(alphabar))*(1 + parms$seasAmp*sin(2*pi*(t/parms$seasPer)))
  dS <- parms$nu - A - parms$mu*S
  dI <- A - alphabar*I - parms$mu*I
  dalphabar <- varalpha*(S*betaPrime(alphabar) - 1) + parms$D*parms$i0
  dvaralpha <- varalpha^2*S*betaPrimePrime(alphabar) + 2*parms$D*(1 - alphabar*parms$i0)
  
  return(list(c(dS=dS, dI=dI, dalphabar=dalphabar, dvaralpha=dvaralpha)))
}

# /* using first order moment equations with fixed variance,
#    and fixed i(0,t)/I(t), assuming they both take their
#    equilibrium values */
momSolve1 <- function(calcParams, initList, paramList, sigma){
  i0 <- 0
  parms <- c(list(sigma=sigma, i0=i0), paramList)
  times <- seq(0,calcParams$Tmax,by=calcParams$Tstep) #times vector for ODE solver
  y <- c(S=initList$S_0, I=initList$I_0, alphabar=initList$mean)
  solution <- lsoda(y, times, momDiffs1, parms)
  S <- solution[,2]
  I <- solution[,3]
  alphabar <- solution[,4]
  return(list(times=times, S=S, I=I, alphabar=alphabar, varalpha=rep(sigma, times=length(times)), title="mom1"))
}

# /* using second order moment equations with
#    and fixed i(0,t)/I(t), assuming it is at the
#    equilibrium value */
momSolve2 <- function(calcParams, initList, paramList){
  i0 <- 0
  parms <- c(list(i0=i0), paramList)
  times <- seq(0,calcParams$Tmax,by=calcParams$Tstep) #times vector for ODE solver
  y <- c(S=initList$S_0, I=initList$I_0, alphabar=initList$mean, varalpha=initList$var)
  solution <- lsoda(y, times, momDiffs2, parms)
  S <- solution[,2]
  I <- solution[,3]
  alphabar <- solution[,4]
  varalpha <- solution[,5]
  return(list(times=times, S=S, I=I, alphabar=alphabar, varalpha=varalpha, title="mom2"))
}

# --- ouput time evolution of (normalized) i(alpha,t) --- #
plotFrames <- function(solutionList, filedir = "./img/", filename="i_alpha-t_", by=5, tmax=Inf, alphamax=Inf, newPlotDevice=T){
  count <- 1
  with(solutionList, {
    ymax <- max(diag(x=1/I,nrow=length(I),ncol=length(I)) %*% i)
    tmax.ind <- tail(which(times < tmax, arr.ind=T), n=1)
    alphamax.ind <- tail(which(alpha < alphamax, arr.ind=T), n=1)
    for(j in seq(1,tmax.ind,by=by)){
      if(newPlotDevice){
        png(filename=paste(filedir,filename,count,".png",sep=""))
      }
      y <- i[j,]/I[j]
      plot(alpha[1:alphamax.ind],y[1:alphamax.ind],typ="l",ylim=c(0,ymax),ylab="i(alpha,t)",
           main=paste("t=",times[j]), xlab="alpha")
      if(newPlotDevice){
        dev.off()
      }
      count <- count + 1
    }
  })
}

plot_S_I_alphabar_varalpha <- function(solutionList, tmax=Inf){
  maxPlotTime.ind <- tail(which(solutionList$times <= tmax, arr.ind=T),n=1)
  timeseq <- 2:maxPlotTime.ind
  with(solutionList, {
    plot(times[timeseq],S[timeseq],type="l",ylim=c(0,max(c(S,I))),ylab=" ",xlab="t",col="red", main="S: red; I: blue", log="x")
    lines(times[timeseq], I[timeseq],col="blue")
    
    plot(times[timeseq],varalpha[timeseq],type="l",xlab="t", ylab="varalpha", log="x")
    plot(times[timeseq],alphabar[timeseq],type="l",xlab="t", ylab="alphabar", log="x")
  })
}

plot_scaled <- function(solutionList, tmax=Inf){
  maxPlotTime.ind <- tail(which(solutionList$times <= tmax, arr.ind=T),n=1)
  timeseq <- 2:maxPlotTime.ind
  with(solutionList, {
    plot(times[timeseq],S[timeseq]/max(S),type="l", ylim=c(0,1),ylab=" ",xlab="t",col="red", 
         main=paste(solutionList$title, "S: red; I: blue; ab: black; va: green", sep="\n", log="x"))
    lines(times[timeseq], I[timeseq]/max(I),col="blue")
    
    lines(times[timeseq],varalpha[timeseq]/max(varalpha),col="green")
    lines(times[timeseq],alphabar[timeseq]/max(alphabar))
  })
}

plot_phaseAlphabar <- function(solutionList){
  with(solutionList, {
    plot(I, alphabar,type="l",ylab="alphabar",xlab="I", main=solutionList$title)
  })
}

plot_phaseVaralpha <- function(solutionList){
  with(solutionList, {
    plot(I, varalpha,type="l",ylab="varalpha",xlab="I", main=solutionList$title)
  })
}


compareTwoSolutions <- function(solutionList1, solutionList2, tMin=1e-15, tMax=Inf){
  Smax <- max(solutionList1$S, solutionList2$S)
  Imax <- max(solutionList1$I, solutionList2$I)
  alphabarmax <- max(solutionList1$alphabar, solutionList2$alphabar)
  varalphamax <- max(solutionList1$varalpha, solutionList2$varalpha)
  
  Smin <- min(solutionList1$S, solutionList2$S)
  Imin <- min(solutionList1$I, solutionList2$I)
  alphabarmin <- min(solutionList1$alphabar, solutionList2$alphabar)
  varalphamin <- min(solutionList1$varalpha, solutionList2$varalpha)

  times.ind <- which(solutionList1$times > tMin & solutionList1$times < tMax)
  maintitle <- paste(paste(solutionList1$title, ": black", sep=""), paste(solutionList2$title, " red", sep=""), sep="; ")
  
  with(solutionList1, {
    plot(times[times.ind], S[times.ind], ylim=c(Smin,Smax), typ="l", xlab="t", ylab="S", log="x")
  })
  with(solutionList2, {
    lines(times[times.ind], S[times.ind], col="red")
  })
  title(main=maintitle)
  
  with(solutionList1, {
    plot(times[times.ind], I[times.ind], ylim=c(Imin,Imax), typ="l", xlab="t", ylab="I", log="x")
  })
  with(solutionList2, {
    lines(times[times.ind], I[times.ind], col="red")
  })
  title(main=maintitle)
  
  with(solutionList1, {
    plot(times[times.ind], alphabar[times.ind], ylim=c(alphabarmin,alphabarmax), typ="l", xlab="t", 
         ylab="alphabar", log="x")
  })
  with(solutionList2, {
    lines(times[times.ind], alphabar[times.ind], col="red")
  })
  title(main=maintitle)
  
  with(solutionList1, {
    plot(times[times.ind], varalpha[times.ind], ylim=c(varalphamin,varalphamax), typ="l", xlab="t", 
         ylab="varalpha", log="x")
  })
  with(solutionList2, {
    lines(times[times.ind], varalpha[times.ind], col="red")
  })
  title(main=maintitle)
}

compareThreeSolutions <- function(solutionList1, solutionList2, solutionList3, tMin=1e-15, tMax=Inf){
  Smax <- max(solutionList1$S, solutionList2$S, solutionList3$S)
  Imax <- max(solutionList1$I, solutionList2$I, solutionList3$I)
  alphabarmax <- max(solutionList1$alphabar, solutionList2$alphabar, solutionList3$alphabar)
  varalphamax <- max(solutionList1$varalpha, solutionList2$varalpha, solutionList3$varalpha)
  
  Smin <- min(solutionList1$S, solutionList2$S, solutionList3$S)
  Imin <- min(solutionList1$I, solutionList2$I, solutionList3$I)
  alphabarmin <- min(solutionList1$alphabar, solutionList2$alphabar, solutionList3$alphabar)
  varalphamin <- min(solutionList1$varalpha, solutionList2$varalpha, solutionList3$varalpha)
  
  times.ind <- which(solutionList1$times > tMin & solutionList1$times < tMax)
  maintitle <- paste(paste(solutionList1$title, ": black", sep=""), paste(solutionList2$title, ": red", sep=""), 
                 paste(solutionList3$title, ": blue", sep=""), sep="; ")
  
  with(solutionList1, {
    plot(times[times.ind], S[times.ind], ylim=c(Smin,Smax), typ="l", xlab="t", ylab="S", log="x")
  })
  with(solutionList2, {
    lines(times[times.ind], S[times.ind], col="red")
  })
  with(solutionList3, {
    lines(times[times.ind], S[times.ind], col="blue")
  })
  title(main=maintitle)
  
  with(solutionList1, {
    plot(times[times.ind], I[times.ind], ylim=c(Imin,Imax), typ="l", xlab="t", ylab="I", log="x")
  })
  with(solutionList2, {
    lines(times[times.ind], I[times.ind], col="red")
  })
  with(solutionList3, {
    lines(times[times.ind], I[times.ind], col="blue")
  })
  title(main=maintitle)
  
  with(solutionList1, {
    plot(times[times.ind], alphabar[times.ind], ylim=c(alphabarmin,alphabarmax), typ="l", xlab="t", 
         ylab="alphabar", log="x")
  })
  with(solutionList2, {
    lines(times[times.ind], alphabar[times.ind], col="red")
  })
  with(solutionList3, {
    lines(times[times.ind], alphabar[times.ind], col="blue")
  })
  title(main=maintitle)
  
  with(solutionList1, {
    plot(times[times.ind], varalpha[times.ind], ylim=c(varalphamin,varalphamax), typ="l", xlab="t",
         ylab="varalpha", log="x")
  })
  with(solutionList2, {
    lines(times[times.ind], varalpha[times.ind], col="red")
  })
  with(solutionList3, {
    lines(times[times.ind], varalpha[times.ind], col="blue")
  })
  title(main=maintitle)
}