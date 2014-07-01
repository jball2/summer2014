require(deSolve)

f <- function(t, y, parms) with( as.list(c(y, parms)), {
  xfer <- b0 * (1 + b1*cos(2*pi*t)) * S * I / N
  dS <- (N - S - I)/L - xfer
  dI <- xfer - I/D
  return(list(c(dS, dI)))
})

y0 <- c(S=500000-2250, I=2250)
parmsA <- c(N=500000, L=4, D= 0.02, b0=500, b1=0.02)
parmsB <- c(N=500000, L=8, D= 0.025, b0=400, b1=0.02)


plotter<-function(parms) {
  out <- ode(y0, seq(0,20,0.1), f, parms)
  ymax <- 2*mean(out[which(out[,'time'] > 10),'I'])
  plot(out[,'time'], out[,'I'], type = "l", xlim=c(10,20), ylim=c(0,ymax), ylab="I", xlab="year")
}

plotter(parmsA)
plotter(parmsB)