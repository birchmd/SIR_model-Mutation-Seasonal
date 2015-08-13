## This program computes the equilibrium distrubution of infectious
## individuals in the the SI model with mutation of the virus.

library("caTools")
library("deSolve")
source("../SI_mut_consts.R")
source("./SI_mut_eqDist_fcns.R")


# --- Constants for calculation --- #
N <- 10000 #number of points to solve ODE at
alphaMin <- 0 #bounds on interval to solve ODE on
alphaMax <- 40
# --------------------------------- #

alpha <- seq(alphaMin, alphaMax, length=N)
SStar <- findSStar(paramList,minOrder=1e-16)
lambda <- compLambda(alpha, SStar, paramList, useInterval=F)
lambda <- lambda/trapz(alpha,lambda) #normalize
alphabar <- trapz(alpha,alpha*lambda)
varalpha <- trapz(alpha,(alpha-alphabar)^2*lambda)

plot(alpha,lambda,typ="l")
abline(0,0)

sink("./eq_output.txt",append=F,split=T)
print(paste("S_*:",SStar))
print(paste("I_*:",IStar_func(alpha,lambda, paramList)))
print(paste("alphabar_*:",alphabar))
print(paste("varalpha_*:",varalpha))
plot(alpha, (SStar*beta(alpha) - alpha - paramList$mu)/paramList$D,typ="l")
abline(0,0)
unlink("./eq_output.txt")
