library("caTools")
library("DEoptim")
library("deSolve")

rootDir <- "C:/Users/Michael/Documents/School/NSERC_USRA_2014/SI_model-Mutation-Seasonal"
if (!(file.exists(paste(sep="/", rootDir, "SI_mut_consts.R")))){
  stop("Root directory not set correctly! Please edit the value in SI_mut_consts.R")
}

# --- Population constants --- #
## Note: pick 1/mu to be one time unit. Pick N_0 := nu/mu
## (i.e. the number of people born in one time unit) to
## be one people unit. In these units nu=mu=1.
nu <- 1 #birth rate
mu <- 1 #per captia death rate
# ---------------------------- #

# --- Disease constants --- #
betaMult <- 3 #beta(alpha) = betaMult*alpha^(1/gamma)
gamma <- 2
diffConst <- 0.05 #diffusion constant (for mutation)
seasAmp <- 0.1 #amplitude of seasonal forcing
seasPer <- 0.1 #period of seasonal forcing
# ------------------------- #

ext <- if (.Platform$OS.type=="unix") "so" else "dll" #extension for compiled shared libraries
compDrvFile <- paste("compiledDerivatives",ext,sep=".")
dyn.load(paste(sep="/", rootDir, "numericImplimentation", compDrvFile))

# --- transmission function --- #
beta <- function(alpha){
  return(betaMult*alpha^(1/gamma))
}

betaPrime <- function(alpha){
  return(betaMult*alpha^(1/gamma - 1)/gamma)
}

betaPrimePrime <- function(alpha){
  return((1-gamma)*betaMult*alpha^(1/gamma - 2)/(gamma^2))
}

paramList <- as.list(c(nu=nu, mu=mu, betaMult=betaMult, gamma=gamma, D=diffConst, seasAmp=seasAmp, seasPer=seasPer))

# --- Set up calculation constants --- #
calcParamsList <- function(alphaMin=1e-15, alphaMax=15, N=2000, Tmax=50, Tstep=0.2){
  return(list(alphaMin=alphaMin, alphaMax=alphaMax, N=N, Tmax=Tmax, Tstep=Tstep))
}

# --- equilibrium values --- #
source(paste(sep="/",rootDir,"endemicEquilibrium/SI_mut_eqDist_fcns.R"))
calcEq <- function(alpha, paramList, exact=F){
  if(exact){
    SStar <- findSStar(paramList,minOrder=1e-16)
    eqDist <- compLambda(alpha, SStar, paramList, useInterval=F)
    eqDist <- eqDist/trapz(alpha,eqDist) #normalize
    IStar <- IStar_func(alpha,eqDist, paramList)
    alphabar_star <- trapz(alpha,alpha*eqDist)
    varalpha_star <- trapz(alpha,(alpha-alphabar_star)^2*eqDist)
  }else{
    alphabar_star <- mom_alphabarStar(paramList)
    IStar <- mom_IStar(alphabar_star, paramList)
    SStar <- mom_SStar(alphabar_star)
    varalpha_star <- mom_varalphaStar(alphabar_star, paramList)
    eqDist <- dnorm(alpha, mean=alphabar_star, sd=sqrt(varalpha_star))
  }
  return(list(S_star=SStar, I_star=IStar, eqDist=eqDist, 
              alphabar_star=alphabar_star, varalpha_star=varalpha_star, alpha=alpha))
}
