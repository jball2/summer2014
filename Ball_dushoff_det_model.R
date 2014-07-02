#Figure 1a from Dynamical Resonance paper
library ("deSolve")

#Parameters for figure 1a
D<- 0.02 #years
L<-4 #duration of immunity in years
beta.0 <- 500 #R0/D, individuals per year
beta.1 <- 0.02 #per year
N<- 500000

#Create vector with parameters above
f1aparams <- c(N,L,D,beta.0,beta.1)


##Create function for beta
      #NOTE: beta naught CANNOT be beta.0 since it's been assigned a value already--> use b0 instead

beta<- function(b0,b1,t) {
  b0 * (1+b1*cos(2*pi*t))
}

##Create the SIRS Model and embed the dsdt and didt functions within it
        #I got sorta stuck here.  My first instinct was to separately define dsdt and didt and then
        #try to throw together a function that used them both.  This is the approach Becky took and it
        #is much more elegant and makes sense.
sirs <- function(N,beta.0,beta.1){
  dsdt <- (N-S-I)/L - beta(beta.0,beta.1,t)*I*S/N  #define dsdt within sirs
  
  didt <- beta(beta.0,beta.1,t)*I*S/N-I/D   #define didt within sirs
  
  list(c(dsdt,didt))      }


#Maxtime and timestep for x axis; I knew to do this part
maxtime<- 50
timestep<- 0.01



##Create time series that solves for number infected
series <-data.frame(lsoda(c(S=499999,I=1),
                          seq(0,maxtime,timestep),sirs,f1aparams))


#####I am getting an error message: "Error in del/by : non-numeric argument to binary operator"

#can't figure out what it means...help?



#next step would be to actually plot the time series across seq(0,maxtime,timestep)















#Becky did MAXTIME<-100 and TIMESTEP<- 0.05
#This will generate a sequence of t from 0 to 100 by 0.05 steps with:
#seq(0,100,0.05)
maxtime<- 100
timestep <- 0.05

#Create SIRS model

sirs<- function(t,y,parms){
    with(c(as.list(y),parms),{
      dsdt <- (N-S-I)/L - beta(beta.0,beta.1,t)*S*I/N #equation 1
      
      didt <- beta(beta.0,beta.1,t)*S*I/N - I/D #equation 2
  
      list(c(dsdt,didt)))}}


###BECKY'S SIRS
sirs <- function(t,y,parms){
  
   with(c(as.list(y),parms),{ 
    
    		
      
      	dsdt <- (N-S-I)/L - beta(t,beta.0,beta.1)*I*S/N 
    
    	didt <- beta(t,beta.0,beta.1)*I*S/N - I/D
    
  		list(c(dsdt,didt))
    
  	})
  
  }

##Time Series
series<- data.frame(lsoda(c(S=499999,I=1),
              seq(0,maxtime,timestep),sirs,f1aparams))
series <-data.frame(lsoda(c(S=499999,I=1),
                          seq(0,maxtime,timestep),sirs,f1aparams))

plot(series$time,series$I,type="l", xlim=c(10,20),ylim=c(0,4000), main="Number Infected",
        xlab="Time(years)", ylab="Number Infected", cex.main=2, cex.lab=1.5, lwd=2, col="green")


##REMOVE EVERYTHING AND JUST USE BECKY'S CODE
remove(list=ls())

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
series <-data.frame(lsoda(c(S=499999,I=1),seq(0,maxtime,timestep),sirs,f1aparams))

plot(ts$time,ts$I,type="l",xlim=c(10,20),ylim=c(0,4000), main="Number infected", #bty="n",
      
      	 xlab="Time (years)",  ylab="", cex.main=2, cex.lab=1.5, cex.axis=1.25, 
      
      	 lwd =3, col="red")


