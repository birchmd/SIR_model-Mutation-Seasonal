mom_alphabarStar_RootFunc <- function(alphabar, parms){
  a <- beta(alphabar)/betaPrime(alphabar)
  b <- -0.5*parms$D*betaPrimePrime(alphabar)/betaPrime(alphabar)
  c <- alphabar + parms$mu
  return(abs(a-sqrt(b)-c))
}

# --- alphabar_* based on moment equations --- #
mom_alphabarStar <- function(parms){
  alphabarStar <- DEoptim(mom_alphabarStar_RootFunc, lower=0, upper=25, DEoptim.control(VTR=0, itermax=300),
                          parms=parms)$optim$bestmem[1]
  return(alphabarStar)
}

# --- varalpha_* based on moment equations --- #
mom_varalphaStar <- function(alphabarStar, parms){
  return(sqrt(-2*parms$D*betaPrime(alphabarStar)/betaPrimePrime(alphabarStar)))
}

mom_SStar <- function(alphabar){
  return(1/betaPrime(alphabar))
}

mom_IStar <- function(alphabar, parms){
  return((parms$nu - parms$mu/betaPrime(alphabar))/(alphabar + parms$mu))
}

# --- ODE system for limiting infectious distribution --- #
diffs <- function(alpha, y, parms, ...){
  dy1 <- y[2]
  dy2 <- (parms$mu + alpha - parms$SStar*beta(alpha))*y[1]/parms$D
  return(list(c(dy1,dy2)))
}


# --- Function to compute S_* --- #
SStar_func <- function(alpha, lambda, paramList){
  alphaBar <- trapz(alpha, alpha*lambda)
  betaBar <- trapz(alpha, beta(alpha)*lambda)
  
  return((paramList$mu + alphaBar)/betaBar)
}

# --- Function to compute I_* --- #
IStar_func <- function(alpha, lambda, paramList){
  alphaBar <- trapz(alpha, alpha*lambda)
  betaBar <- trapz(alpha, beta(alpha)*lambda)
  
  return(paramList$nu/(paramList$mu + alphaBar) - paramList$mu/betaBar)
}

# --- Function to compute equilibrium distribution of infectious individuals (lambda) --- #
compLambda <- function(alpha, SStar, paramList, useInterval=T, Nalpha=10000){
  parms <- list(mu=paramList$mu, D=paramList$D, SStar=SStar)
  if (useInterval){
    y <- c(dnorm(-5),0)
    solution <- lsoda(y, alpha, diffs, parms)
    lambda <- solution[,2]
  }else{
    alphabarStar <- mom_alphabarStar(paramList)
    varalphaStar <- mom_varalphaStar(alphabarStar, paramList)
    alpha2 <- seq(max(alphabarStar-5*sqrt(varalphaStar), 0), alphabarStar+5*sqrt(varalphaStar), len=Nalpha)
    y <- c(dnorm(-5),0)
    solution <- lsoda(y, alpha2, diffs, parms)
    lambda <- solution[,2]
    lambda <- approx(x=alpha2, y=lambda, xout=alpha)$y
    lambda[which(is.na(lambda))] <- 0
  }
  return(lambda)
}

approx2 <- function(x, y, xout){
  
}

# --- returns true if vector has negative elements --- #
hasNeg <- function(x){
  return(length(which(x < 0)) > 0)
}

# --- Compute S_* by the "wag the dog method" (Griffiths QM, p. 54-55) --- #
findSStar <- function(paramList, tol=1e-5, minOrder=1e-15, Nalpha=10000){
  #minOrder is smallest order of magnitude to adjust
  order <- 1 #order of magnitude currently adjusting
  alphabarStar <- mom_alphabarStar(paramList)
  varalphaStar <- mom_varalphaStar(alphabarStar, paramList)
  alpha <- seq(max(alphabarStar-5*sqrt(varalphaStar), 0), alphabarStar+5*sqrt(varalphaStar), len=Nalpha)
  SStar <- 1/betaPrime(alphabarStar)
  n <- length(alpha)/10 #length of tail
  
  
  while(T){
    lambda <- compLambda(alpha, SStar, paramList)
    if(hasNeg(lambda)){ #oscillating solutions, SStar is too big
      SStar <- SStar - order
    }else{
      lambda <- lambda/trapz(alpha,lambda) #normalize
      if(max(tail(lambda,n=n)) < tol){ #converged!
        return(SStar)
      }else{ #too big, reduce adjustment order order
        SStar <- SStar + order
        order <- 0.1*order
        if(order < minOrder){
          #print(SStar,digits=15)
          warning("Adjustment order too small!")
          return(SStar - 10*order)
        }
        SStar <- SStar - order
        #print(SStar,digits=15)
      }
    }
  }
}