# Replication of deterministic model from Dushoff et al (2004)
# JL 06.28.14

library(deSolve)

# Notes from the paper:
# Seasonal oscillations in influenza incidence presumably due to some mechanism
# that causes seasonal changes in the effective transmission rate from person to person
# We show that it may be impossible to establish the underlying cause of seasonality
# in epidemics because the large observed oscillations in incidence can be generated
# by seasonal changes in the transmission rate that are too small to measure.
# if the period of endogenous SIRS oscillations (T) is near the period of seasonal forcing
# (1 year) then these two factors may resonate to produce greatly amplified oscillations in 
# incidence

# Deterministic model
# det.mod is a function with arguments t,y and params.
# the input to argument y is a vector of the state variables-these are variables we are
# interested in and that change over time, so the input of y is a vector containing
# the initial number of: S and I hosts
det.mod<-function(t,y,params){ 
    S <- y[1] # susceptible hosts
    I <- y[2] # infectious hosts
    
# the input to the params argument is also a vector. parameters are the variables that
# influence the state variables, including number of hosts, the transmission
# rate- beta etc

# parameters
    N<-params[1] # total population size
    beta0<-params[2]
		beta1<-params[3]
    dur.inf<-params[4] # d
    dur.imm<-params[5] # l

# time-dependent beta
beta.t.func<-function(b0,b1,time){
	beta.t<-b0*(1 + b1*cos(2*pi*time))
	return(beta.t)
}

beta.t<-beta.t.func(b0=beta0,b1=beta1,time=t)

# ODE's
    dS.dt <- (N-S-I)/dur.imm - beta.t*I*S/N
    dI.dt <- beta.t*I*S/N - I/dur.inf
    # list containing the derivatives
    xdot <- c(dS.dt,dI.dt)
    return(list(xdot))
}

flu.sim<-function(N,beta0,beta1,dur.imm,dur.inf,i0=1/N,t.max=20,ts=0.001){
  I.0<-i0*N
  S.0<-N-i0*N
  # timesteps
  times<-seq(0,t.max,ts) 
  # initial conditions
  init<-c(sus=S.0,inf=I.0)
  return(data.frame(lsoda(init,times,det.mod,c(N,beta0,beta1,dur.inf,dur.imm))))
}

# a weak resonance
a<-flu.sim(N=500000,dur.imm=4,dur.inf=0.02,beta0=500,beta1=0.02)
plot(a$time, a$inf,ylim=c(0,4000),xlim=c(10,20),type="l",ylab="Number infected",xlab="Time (years)")

# b strong resonance
b<-flu.sim(N=500000,dur.imm=8,dur.inf=0.025,beta0=400,beta1=0.02)
plot(b$time, b$inf,ylim=c(0,4000),xlim=c(10,20),type="l",ylab="Number infected",xlab="Time (years)")

# results qualtitatively similar but not as big a difference between a and b

#d<- 6 7,8,9,10 days mean infectious period
#l<- 4 5,6,7,8 years average duration of immunity
#R0<- 4-16 
# In a SIRS model, the intrinsic period of oscillation is approximately: 
int.per.func<-function(d,l,R0){
	t= 2*pi*sqrt((d*l)/(R0-1))
	return(t)
}
# if the period of endogenous SIRS oscillations (T) is near the period of seasonal forcing
# (1 year) then these two factors may resonate to produce greatly amplified oscillations in 
# incidence

# check calcs of T for a and b
# R0 is dur.inf*beta0
0.02*500
a.t<-int.per.func(0.02,4,10)
0.025*400
b.t<-int.per.func(0.025,8,10)


