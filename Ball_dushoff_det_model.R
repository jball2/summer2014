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








