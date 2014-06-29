# Functions deterministic model in Dushoff et al. 2004
# June 2014

## SIRS model for seasonal infection model

library(deSolve)

# Parameter values taken from Fig 1a

parms1a <- c(N=500000,          # total population size
					 L=4,								# average duration of immunity (in years)
					 D=0.02, 						# mean infectious period (in years)
					 beta.0=500,				# R0/D (individuals/year)
					 beta.1=0.02				# scaling factor for contact function
)

# Seasonal contact function 
Beta <- function(x,B0,B1){
	B0*(1+B1*cos(2*pi*x))
}

MAXTIME <- 100
TIMESTEP <- 0.005

sirs <- function(t,y,parms){
	with(c(as.list(y),parms),{ 
			
		dSdt <- (N-S-I)/L - Beta(t,beta.0,beta.1)*I*S/N 
		dIdt <- Beta(t,beta.0,beta.1)*I*S/N - I/D
		list(c(dSdt,dIdt))
	})
}


# time series
ts <- data.frame(lsoda(c(S=499999,I=1),seq(0,MAXTIME,TIMESTEP),sirs,parms1a)) 
plot(ts$time,ts$I,type="l",xlim=c(10,20),ylim=c(0,4000), main="Number infected", #bty="n",
		 xlab="Time (years)",  ylab="", cex.main=2, cex.lab=1.5, cex.axis=1.25, 
		 lwd =3, col="red")




