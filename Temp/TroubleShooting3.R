A <- 1
fi <- 1
ff <- 0.1
t <- seq(0,300, 0.1)
sweepDuration <- 300
N <- log(ff/fi)/log(2)
R <- N / (sweepDuration) # Usually tf-ti but time adjusting t by ti
phi <- 0
b <- 0

ti <- 0
predicted <- A*sin(2*pi*((fi*(-1+2^(R*(t-ti))))/(R*log(2))) - phi) + b

plot(t, predicted, type='l')

A <- 1
fi <- 1
ff <- 0.1
t <- seq(0,300, 0.1)
sweepDuration <- 300
N <- log(ff/fi)/log(2)
R <- N / (sweepDuration) # Usually tf-ti but time adjusting t by ti
phi <- 0
b <- 0

ti <- -1
predicted <- A*sin(2*pi*((fi*(-1+2^(R*(t-ti))))/(R*log(2))) - phi) + b
lines(t, predicted, col='red')
