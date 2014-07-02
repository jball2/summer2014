# dushoff2004.R
# JRCP 01 July 2014

rm(list=ls())
require(deSolve)

## sir() -- Function for numerical analysis of model 1
sir <- function(t,y,parms){
  # The with() function gives access to the named values of parms within the
  # local environment created by the function
  with(c(as.list(y),parms),{
    
    transRate <- beta.t(t,beta0,beta1)*S*I/N
    
    dSdt <- (N-S-I)/LL-transRate
    dIdt <- transRate - I/DD
    # Note: Population size is constant, so don't need to specify dRdt
    list(c(dSdt,dIdt))
  })
}

beta.t <- function(time,beta0,beta1){
  beta0 * (1 + beta1*cos(2*pi*time))
}

R0.t <- function(time,DD,beta0,beta1){
  beta.t(time,beta0,beta1)*DD
}

endemicEq <- function(params){
  with(as.list(params),{
    S.star <- N/(beta0*DD)
    I.star <- (N-S.star)/(1+LL/DD)
    return(c(S=S.star,I=I.star))
  })
}

intrinsicPeriod <- function(params){
  with(as.list(params),
       2*pi*sqrt(DD*LL/(DD*beta0-1))
  )}

# NOTE: methods say value of beta1 used is 0.02, but Figure 1 appears to have been produced using
# a value of 0.04

fig1A <- c(N = 500000,
           LL = 4,
           DD = 0.02,
           beta0 = 500,
           beta1 = 0.04)
init1A <- endemicEq(fig1A)

fig1B <- c(N = 500000,
           LL = 8,
           DD = 0.025,
           beta0 = 400,
           beta1 = 0.04)
init1B <- endemicEq(fig1B)


runSIR <- function(params,time.out = seq(0,20,.05),plot=T,...){
  ts <- lsoda(
    y = endemicEq(params),        # Initial conditions for population
    times = time.out,             # Timepoints for evaluation
    func = sir,                   # Function to evaluate
    parms = params                # Vector of parameters
  )
  out <- list(ts = data.frame(ts),params = params)
  if(plot){
    plotSIR(out,compare=T,...)
  }
  return(out)
}

plotSIR <- function(run,compare=F,ymax=4000,...){
  with(run,{
    plot(ts$time,ts$I,
         xlim=c(10,20),ylim=c(0,ymax),
         xlab="Time (years)",ylab="Number infected",         
         type="l",bty="n",...)
    axis(2,seq(0,ymax,ymax/8),labels=NA,tcl=-.2)
    if(compare){
      with(as.list(params),{
        S.star.t <- N/R0.t(ts$time,DD,beta0,beta1)
        I.star.t <- (N-S.star.t)/(1+LL/DD)
        lines(ts$time,I.star.t,lwd=3)
      })}
  })
}

par(mfcol=c(2,1),mar=c(2,2,1,0))
ts1A <- runSIR(fig1A,col="red",lwd=3)
ts1B <- runSIR(fig1B,col="red",lwd=3)

intrinsicPeriod(fig1A)
intrinsicPeriod(fig1B)

