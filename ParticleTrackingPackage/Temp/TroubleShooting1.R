


voltage <- 1
DACMAX <- 3.3
A <- 0.8            # Amplitude
offset <- 0 #double A plus a small amount
tf <- 300         ## of sec
pi <- 3.14159
fi <- 1           #frequency initial
ff <- .01
N <- ((log(ff/fi))/(log(2)))
R <- (N/tf)
P <- 0.0
X <- 0.0
myT <- seq(0,300000,length.out=30000)

val <- c()
for(t in myT)
{
     if (t < 1000*tf)
     {
          #t <- (millis() - start)
          P <- 2^(R*((t)/1000))
          X <- (sin(2*pi*((fi*(-1+(P)))/(R*log(2)))))
          voltage<- ((-A)*(asin(X)))+offset
          outputValue <- round((4095.0/DACMAX)*voltage+2048.0)
          val <- c(val, outputValue)#analogWrite(DAC1, outputValue)
          #//delay(2)
     }
     else
     {
          t <- 1000*tf
          P <- 2^(R*((t)/1000))
          X <- (sin(2*pi*((fi*(-1+(P)))/(R*log(2)))))
          voltage<- ((-A)*(asin(X)))+offset
          outputValue <- round((4095.0/DACMAX)*voltage+2048.0)
          val <- c(val, outputValue)#analogWrite(DAC1, outputValue)
     }
}



bf <- butter(3, 0.1)                        # 10 Hz low-pass filter
t <- trackList$meta$allTimes
x <- colMeans(trackList$getMatrix(), na.rm=T)
y <- filtfilt(bf, x)

duh <- getSweep(sin=T, amplitude=2048, fi=1, ff=0.01, ti=-2, t=trackList$meta$allTimes)
duh2 <- getSweep(sin=F, amplitude=2048, fi=1, ff=0.01, ti=-2, t=trackList$meta$allTimes)
plot(duh$t, -duh$v, col="red", lwd=1, type='l')#, xlim=c(0,5))
lines(duh2$t, -duh2$v, col="blue", lwd=1)
abline(h=0)

plot(myT/1000, getDerivative(val-2048, myT/1000), type='l')
lines(t, y*50, col="blue", lwd=1)
lines(duh$t, -duh$v, col="red", lwd=1)
abline(h=0)



plot(myT/1000, getDerivative(val-2048, myT/1000), type='l', xlim=c(0,5))
lines(t, y*50, col="blue", lwd=1)
lines(duh$t, -duh$v, col="red", lwd=1)
abline(h=0)

plot(myT/1000, getDerivative(val-2048, myT/1000), type='l', xlim=c(150,300))
lines(t, y*50, col="blue", lwd=1)
lines(duh$t, -duh$v, col="red", lwd=1)
abline(h=0)

