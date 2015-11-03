



library("deSolve")
parameters <- c(a = 0.25,
                b = 4.0 ,
                F1=8.0,
                G1=1.0
)

px=runif(1,0,1)
py=runif(1,0,1)
pz=runif(1,0,1)

#px=0
#py=0
#pz=0


state <- c(X = 1.478+px,
           Y = 1.0+py,
           Z = 0.0+pz)
  

Lorenz<-function(t, state, parameters) {
  with(as.list(c(state, parameters)),{
    # rate of change
    dX <- -Y^2 - Z^2-a*X + a*F1
    dY <- X*Y - b*X*Z - Y + G1
    dZ <- b*X*Y + X*Z - Z     
    # return the rate of change
    list(c(dX, dY, dZ))
  })   # end with(as.list ...
}


times <- seq(0, 10000, by = 1)

out <- ode(y = state, times = times, func = Lorenz, parms = parameters)
#head(out)

layout(matrix(c(1,2,3),nrow=3))
#par(oma = c(0, 0, 3, 0))
#plot(out, xlab = "time", ylab = "-")
plot(out[,"X"],xlab = "time", ylab = "X")
plot(out[,"Y"],xlab = "time", ylab = "Y")
plot(out[,"Z"],xlab = "time", ylab = "Z")
mtext(outer = TRUE, side = 3, "Lorenz model G=1.0", cex = 1.5)

mean(out[,"X"])
mean(out[,"Y"])
mean(out[,"Z"])













