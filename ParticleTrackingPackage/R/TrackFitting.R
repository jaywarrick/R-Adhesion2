#' @import NMOF
NULL

#' sseBulk
#'
#' @param trackList A TrackList object
#' @param trackMatrix The trackMatrix of the TrackList object (so it doesn't need to be calculated for every call to this function)
#' @param amplitude A numeric value of the amplitude of the sweep function.
#' @param phaseShift A numeric value of the phase shift of the sweep function.
#' @param timeScalingFactor A numeric value of the time scaling factor
#' @param ti A numeric value of the initial time of the first frame
#'
#' @return the sum square error between the TrackList data and a sweep with the provided amplitude and phaseShift
sseBulk <- function(trackList, trackMatrix, amplitude, phaseShift, timeScalingFactor, ti)
{
     # Cols are frames, rows are track ids
     # we pass the 'trackMatrix' so we don't have to obtain it each iteration

     # get the sweep (for this we need the 'trackList')
     sweep <- getSweep(amplitude=amplitude, phaseShift=phaseShift, offset=pi, sin=trackList$meta$sin, ti=ti, fi=trackList$meta$fi, ff=trackList$meta$ff, sweepDuration=timeScalingFactor*trackList$meta$sweepDuration, t=trackList$meta$tAll, guess=NULL)
     # For each index in tAll (i.e., for each frame)
     sse <- sum((t(trackMatrix)-sweep$v)^2, na.rm=TRUE) # Do the transpose because the subtract function typically makes the subtracted vector vertical
     cat("(", amplitude, ",", timeScalingFactor, ",", ti, ") = ", sse, "\n", sep="")
     return(sse)
}

#' seBulkGS
#'
#' @param x default
#' @param trackList default
#' @param trackMatrix default
#' @param amplitude default
#' @param timeScalingFactor default
#'
#' @return the sum square error between the TrackList data and a sweep with the provided amplitude and phaseShift
sseBulkGS <- function(x, trackList, trackMatrix, amplitude, timeScalingFactor)
{


     # Cols are frames, rows are track ids
     # we pass the 'trackMatrix' so we don't have to obtain it each iteration

     # get the sweep (for this we need the 'trackList')
     sweep <- getSweep(amplitude=amplitude, phaseShift=x[2], offset=pi, sin=trackList$meta$sin, ti=x[1], fi=trackList$meta$fi, ff=trackList$meta$ff, sweepDuration=timeScalingFactor*trackList$meta$sweepDuration, t=trackList$meta$tAll, guess=NULL)
     # For each index in tAll (i.e., for each frame)

     #sse <- sum((t(trackMatrix)-sweep$v)^2, na.rm=TRUE) # Do the transpose because the subtract function typically makes the subtracted vector vertical
     sse <- mad((t(trackMatrix)-sweep$v), na.rm=TRUE) # Do the transpose because the subtract function typically makes the subtracted vector vertical

     cat("(", x[1], ",", x[2], ") = ", sse, "\n", sep="")
     return(sse)
}

#' sseBulk2
#'
#' @param trackList default
#' @param trackMatrix default
#' @param amplitude default
#' @param phaseShift default
#' @param ti default
#' @param fi default
#' @param ff default
#'
#' @return the sum square error between the TrackList data and a sweep with the provided amplitude and phaseShift
sseBulk2 <- function(trackList, trackMatrix, amplitude, phaseShift, ti, fi, ff)
{
     # Cols are frames, rows are track ids
     # we pass the 'trackMatrix' so we don't have to obtain it each iteration
     timeScalingFactor = 1
     # get the sweep (for this we need the 'trackList')
     sweep <- getSweep(amplitude=amplitude, phaseShift=phaseShift, offset=pi, sin=trackList$meta$sin, ti, fi=fi, ff=ff, sweepDuration=timeScalingFactor*trackList$meta$sweepDuration, t=trackList$meta$tAll, guess=NULL)
     # For each index in tAll (i.e., for each frame)
     sse <- sum((t(trackMatrix)-sweep$v)^2, na.rm=TRUE) # Do the transpose because the subtract function typically makes the subtracted vector vertical
     cat("(", amplitude, ",", 1.00, ",", ti, ",", fi, ",", ff, ") = ", sse, "\n", sep="")
     return(sse)
}

#' Fit the data to determine the phase shift
#'
#' Guess the other parameters from the average of the data.
#'
#' @param trackList A TrackList object
#' @param tiGuess A numeric value guessing the frame at which time actually equals zero (i.e., the fluid frequency sweep is started).
getBulkPhaseShift <- function(trackList, tiGuess=0)
{
     phaseShift <- pi
     # require(stats)
     trackMatrix <- trackList$getMatrix()
     amplitude <- max(as.numeric(trackList$getProp(fun=function(x){r <- x$range('x', rel=TRUE); r <- (r[2]-r[1])/5; return(r)})))
     guess <- c(timeScalingFactor=1, ti=tiGuess)
     #guess <- c(amplitude=amplitude, phaseShift=1.2*pi)
     #amplitudeLimits=c(0.1*amplitude, amplitude)
     #phaseShiftLimits=c(1*pi,1.4*pi)
     timeScalingFactorLimits=c(0.99, 1.01)
     tiLimits=c(-0.5/trackList$meta$fi, 0.5/trackList$meta$fi)
     #sweepDurationLimits=c(0.98*trackList$meta$sweepDuration, 1.02*trackList$meta$sweepDuration)
     #fiLimits=c(0.98*trackList$meta$fi, 1.02*trackList$meta$fi)
     #ffLimits=c(0.98*trackList$meta$ff, 1.02*trackList$meta$ff)
     bestFit <- optim(par=guess,
                      function(par, trackList, trackMatrix){sseBulk(trackList=trackList, trackMatrix=trackMatrix, amplitude=amplitude, phaseShift=phaseShift, timeScalingFactor=par['timeScalingFactor'], ti=par['ti'])},
                      method='L-BFGS-B',
                      #method='SANN',
                      lower=c(min(timeScalingFactorLimits), min(tiLimits)),
                      upper=c(max(timeScalingFactorLimits), max(tiLimits)),
                      control=list(trace=0),
                      trackList=trackList,
                      trackMatrix=trackMatrix)
     return(list(par=c(phaseShift=phaseShift, amplitude=amplitude, timeScalingFactor=as.numeric(bestFit$par['timeScalingFactor']), ti=as.numeric(bestFit$par['ti']), offset=as.numeric(bestFit$par['offset'])), fit=bestFit))
}

#' Fit the data to determine the phase shift (version 2)
#'
#' Guess the other parameters from the average of the data.
#'
#' @param trackList A TrackList object
#' @param tiGuess A numeric value guessing the frame at which time actually equals zero (i.e., the fluid frequency sweep is started).
getBulkPhaseShift2 <- function(trackList, tiGuess=0)
{
     phaseShift <- pi
     #require(stats)
     trackMatrix <- trackList$getMatrix()
     amplitude <- max(as.numeric(trackList$getProp(fun=function(x){r <- x$range('x', rel=TRUE); r <- (r[2]-r[1])/5; return(r)})))
     guess <- c(ti=tiGuess, fi=trackList$meta$fi, ff=trackList$meta$ff)
     #guess <- c(amplitude=amplitude, phaseShift=1.2*pi)
     #amplitudeLimits=c(0.1*amplitude, amplitude)
     #phaseShiftLimits=c(1*pi,1.4*pi)
     #timeScalingFactorLimits=c(0.99, 1.01)
     tiLimits=c(-0.55/trackList$meta$fi, 0.55/trackList$meta$fi)
     #sweepDurationLimits=c(0.98*trackList$meta$sweepDuration, 1.02*trackList$meta$sweepDuration)
     fiLimits=c(0.98*trackList$meta$fi, 1.05*trackList$meta$fi)
     ffLimits=c(0.98*trackList$meta$ff, 1.05*trackList$meta$ff)
     bestFit <- optim(par=guess,
                      function(par, trackList, trackMatrix){sseBulk2(trackList=trackList, trackMatrix=trackMatrix, amplitude=amplitude, phaseShift=phaseShift, ti=par['ti'], fi=par['fi'], ff=par['ff'])},
                      method='L-BFGS-B',
                      #method='SANN',
                      lower=c(min(tiLimits), min(fiLimits), min(ffLimits)),
                      upper=c(max(tiLimits), max(fiLimits), max(ffLimits)),
                      control=list(trace=0),
                      trackList=trackList,
                      trackMatrix=trackMatrix)
     return(list(par=c(phaseShift=phaseShift, amplitude=amplitude, timeScalingFactor=1.00, ti=as.numeric(bestFit$par['ti']), fi=as.numeric(bestFit$par['fi']), ff=as.numeric(bestFit$par['ff']), offset=as.numeric(bestFit$par['offset'])), fit=bestFit))
}


#' Fit the data to determine the phase shift (version GridSearch, GS)
#'
#' Guess the other parameters from the average of the data.
#'
#' @param ti default
#' @param phaseShift default
#' @param cores default
#' @param trackList A TrackList object
getBulkPhaseShiftGS <- function(trackList, ti=seq(-1,1,1/30), phaseShift=seq(-pi,pi,pi/30), cores=1)
{
     trackMatrix <- trackList$getMatrix()

     #Choose decent amplitude
     amplitude <- max(as.numeric(trackList$getProp(fun=function(x){r <- x$range('x', rel=TRUE); r <- (r[2]-r[1])/5; return(r)})))

     #Time scaling factor is basically useless, we can be about 0.015 seconds off by guessing 0.035s frame rate after 300s (i.e, a half a frame).
     #Vary phaseShift, ti

     if(cores > 1)
     {
          mc.control <- list(mc.cores=cores)
          res <- gridSearch(sseBulkGS, levels=list(ti=ti, phaseShift=phaseShift), trackList=trackList, trackMatrix=trackMatrix, amplitude=amplitude, timeScalingFactor=1, method='multicore', mc.control=mc.control)
     }
     else
     {
          res <- gridSearch(sseBulkGS, levels=list(ti=ti, phaseShift=phaseShift), trackList=trackList, trackMatrix=trackMatrix, amplitude=amplitude, timeScalingFactor=1, method='loop')
     }

     finalPhaseShift <- res$minlevels[2]
     finalTi=res$minlevels[1]
     limitFlag <- 0
     if(finalPhaseShift %in% range(ti) | finalTi %in% range(phaseShift))
     {
     	limitFlag <- 1
     }

     return(list(par=c(phaseShift=finalPhaseShift, amplitude=amplitude, timeScalingFactor=1, ti=finalTi, offset=0, limitFlag=limitFlag), fit=NULL, errorResults=res))
}

#' sseLogNormGS1
#'
#' @param x default
#' @param dataX default
#' @param dataY default
#'
#' @return the sum square error between the data and a logNorm cumulative curve
sseLogNormGS1 <- function(x, dataX, dataY)
{
     temp <- logNorm(x=dataX, mu=x[1], sigma=x[2])
     # For each index in tAll (i.e., for each frame)
     sse <- sum((dataY-temp)^2, na.rm=TRUE) # Do the transpose because the subtract function typically makes the subtracted vector vertical
     cat("(", x[1], ",", x[2], ") = ", sse, "\n", sep="")
     return(sse)
}

#' sseLogNormGS2
#'
#' @param x default
#' @param dataX default
#' @param dataY default
#'
#' @return the sum square error between the data and a logNorm cumulative curve
sseLogNormGS2 <- function(x, dataX, dataY)
{
     temp <- logNorm2(dataX, x[1], x[2], x[3], x[4], x[5])
     # For each index in tAll (i.e., for each frame)
     sse <- sum((dataY-temp)^2, na.rm=TRUE) # Do the transpose because the subtract function typically makes the subtracted vector vertical
     cat("(", x[1], ",", x[2], ",", x[3], ",", x[4], ",", x[5], ") = ", sse, "\n", sep="")
     return(sse)
}

#' sseLogNorm1
#'
#' @param x default
#' @param y default
#' @param mu default
#' @param sigma default
sseLogNorm1 <- function(x, y, mu, sigma)
{
     temp <- logNorm(x=x, mu=mu, sigma=sigma)
     # For each index in tAll (i.e., for each frame)
     sse <- sum((y-temp)^2, na.rm=TRUE) # Do the transpose because the subtract function typically makes the subtracted vector vertical
     cat("(", mu, ",", sigma, ") = ", sse, "\n", sep="")
     return(sse)
}

#' sseLogNorm
#'
#' @param x default
#' @param mu default
#' @param sigma default
logNorm <- function(x=seq(0,10,0.1), mu=1, sigma=1)
{
     return(1-(1/2+(1/2)*erf((log(x)-log(mu))/(sqrt(2*sigma)))))
}

#' sseLogNorm2GS
#'
#' @param x default
#' @param dataX default
#' @param dataY default
sseLogNorm2GS <- function(x, dataX, dataY)
{
     # If gridSearch, then variables are all passed in via x
     temp <- logNorm2(x=dataX, alpha=x[1], mu1=x[2], sigma1=x[3], mu2=x[4], sigma2=x[5])
     # For each index in tAll (i.e., for each frame)
     sse <- sum((dataY-temp)^2, na.rm=TRUE) # Do the transpose because the subtract function typically makes the subtracted vector vertical
     cat("(", x[1], ",", x[2], ",", x[3], "," , x[4], ",", x[5], ") = ", sse, "\n", sep="")
     return(sse)
}

#' sseLogNorm2
#'
#' @param x default
#' @param y default
#' @param alpha default
#' @param mu1 default
#' @param sigma1 default
#' @param mu2 default
#' @param sigma2 default
sseLogNorm2 <- function(x, y, alpha, mu1, sigma1, mu2, sigma2)
{
     # If gridSearch, then variables are all passed in via x
     temp <- logNorm2(x=x, alpha=alpha, mu1=mu1, sigma1=sigma1, mu2=mu2, sigma2=sigma2)
     # For each index in tAll (i.e., for each frame)
     sse <- sum((y-temp)^2, na.rm=TRUE) # Do the transpose because the subtract function typically makes the subtracted vector vertical
     cat("(", alpha, ",", mu1, ",", sigma1, "," , mu2, ",", sigma2, ") = ", sse, "\n", sep="")
     return(sse)
}

#' logNorm2
#'
#' @param x default
#' @param alpha default
#' @param mu1 default
#' @param sigma1 default
#' @param mu2 default
#' @param sigma2 default
#'
#' @export
logNorm2 <- function(x=seq(0,10,0.1), alpha=0.5, mu1=1, sigma1=1, mu2=1.5, sigma2=1.5)
{
     return(alpha*(1-(1/2+(1/2)*erf((log(x)-log(mu1))/(sqrt(2*sigma1))))) + (1-alpha)*(1-(1/2+(1/2)*erf((log(x)-log(mu2))/(sqrt(2*sigma2))))))
}

#' getGuessGS
#'
#' @param x default
#' @param y default
#' @param mu default
#' @param sigma default
#' @param cores default
#'
#' @export
getGuessGS <- function(x, y, mu=lseq(0.0001, 0.5, 50), sigma=seq(0.1,5,0.1), cores=1)
{
     mc.control <- list(mc.cores=cores)
     levels <- list(mu=mu, sigma=sigma)
     res <- gridSearch(sseLogNormGS1, levels, dataX=x, dataY=y, method='multicore', mc.control=mc.control)
     muRet <- res$minlevels[1]
     sigmaRet <- res$minlevels[2]
     limitFlag <- FALSE
     if(muRet == min(mu) || muRet == max(mu))
     {
          limitFlag = TRUE
     }
     if(sigmaRet == min(sigma) || sigmaRet == max(sigma))
     {
          limitFlag = TRUE
     }
     return(list(mu=muRet, sigma=sigmaRet, limitFlag=limitFlag))
}

#' getGuess
#'
#' @param x default
#' @param y default
#' @param thresh default
#' @param max default
#' @param fraction default
#' @param cores default
#'
#' @export
getGuess <- function(x, y, thresh, max=1, fraction=0.2, cores=1)
{
     mu <- lseq(0.0001,0.5,100)
     sigma <- seq(0.1,5,length.out=50)
     levels <- list(mu=mu, sigma=sigma)
     mc.control <- list(mc.cores=cores)
     res <- gridSearch(sseLogNormGS, levels, dataX=x, dataY=y, method='multicore', mc.control=mc.control)

     # get a guess for tau based upon hints
     temp <- smooth.spline(x=tau, y=y, spar=0.5)
     temp2 <- spline(x=temp$x, y=temp$y)
     tau = which(temp2$y <= thresh)
     if(isempty(tau))
     {
          tau = which.max(temp2$y)
     }
     else
     {
          tau <- tau[1]
     }
     if(temp2$y[tau] > max)
     {
          tau = which(temp2$y <= max*frac)
     }
     tau <- temp2$x[tau]

     # Now get sigma
     theRange <- range(temp2$y)
     rangeMiddle <- mean(theRange)
     hi <- mean(c(rangeMiddle, theRange[2]))
     lo <- mean(c(rangeMiddle, theRange[1]))
     hiN <- which(temp2$y <= hi)[1]
     loN <- which(temp2$y <= lo)[1]

     # This is the fraction of sigma=1 expected for this difference in cumulative distribution
     idealDiff <- (qnorm(temp2$y[hiN], mean=0, sd=1) - qnorm(temp2$y[loN], mean=0, sd=1))
     if(idealDiff == 0)
     {
          stop("Couldn't calculate a valid 'ideal' sigma change because likely the chosen datapoints are the same.")
     }

     # This is the actual difference observed for this difference in cumulative distribution
     actualDiff <- abs(log(temp2$x[hiN]) - log(temp2$x[loN]))

     # Therefore, if idealDiff corresponds to a sigma of 1, the sigma for actualDiff is just actualDiff/idealDiff
     sigma <- actualDiff/idealDiff

     return(list(mu=tau, sigma=sigma, minLevels=res$minlevels))
}

#' fitLogNorm1
#'
#' @param x default
#' @param y default
#' @param guess default
#' @param method default
#'
#' @export
fitLogNorm1 <- function(x, y, guess, method='L-BFGS-B')
{
     bestFit <- optim(par=guess,
                      function(par, x, y){sseLogNorm1(x=x, y=y, mu=par['mu'], sigma=par['sigma'])},
                      method=method,
                      lower=c(1e-8,1e-2),
                      #  upper=c(max(muLimits), max(sigmaLimits)),
                      control=list(trace=0),
                      x=x,
                      y=y)
     sst <- sum((y-mean(y))^2)
     sse <- sseLogNorm1(x=x, y=y, mu=bestFit$par['mu'], sigma=bestFit$par['sigma'])
     r2 <- 1-(sse/sst)
     #      print(paste0('guess = ', guess))
     #      print(paste0('muLimits = ', muLimits))
     #      print(paste0('sigmaLimits = ', sigmaLimits))
     # print(bestFit$par)
     return(list(par=c(mu=bestFit$par[['mu']], sigma=bestFit$par[['sigma']]), fit=bestFit, sse=sse, sst=sst, r2=r2))
}

#' fitLogNorm2
#'
#' @param x default
#' @param y default
#' @param guess default
#' @param alpha default
#' @param mu2Guess default
#' @param sigma2Guess default
#' @param method default
#'
#' @export
fitLogNorm2 <- function(x, y, guess, alpha, mu2Guess, sigma2Guess, method='L-BFGS-B')
{
     sse_single <- sseLogNorm1(x, y, mu=guess[['mu']], sigma=guess[['sigma']])

     #      alphaLimits <- c(0.01,0.99)
     #      sigma1Limits=as.numeric(c(guess[['sigma']]*(1-alpha), guess[['sigma']]*(1/(1-alpha))))
     #      sigma2Limits=as.numeric(c(newSigma*alpha, newSigma*(1/alpha)))
     #      mu1Limits=as.numeric(c(guess[['mu']]*(1-alpha), guess[['mu']]*(1/(1-alpha))))
     #      mu2Limits=as.numeric(c(newMu*alpha, newMu*(1/alpha)))

     guess <- c(alpha=alpha, mu1=guess[['mu']], sigma1=guess[['sigma']], mu2=mu2Guess, sigma2=sigma2Guess)
     bestFit <- optim(par=guess,
                      function(par, x, y){sseLogNorm2(x=x, y=y, alpha=par['alpha'], mu1=par['mu1'], sigma1=par['sigma1'], mu2=par['mu2'], sigma2=par['sigma2'])},
                      method=method,
                      lower=c(0,1e-8,1e-2,1e-8,1e-2),
                      #upper=c(max(alphaLimits), max(mu1Limits), max(sigma1Limits), max(mu2Limits), max(sigma2Limits)),
                      control=list(trace=0),
                      x=x,
                      y=y)
     sst <- sum((y-mean(y))^2)
     sse <- sseLogNorm2(x=x, y=y, alpha=bestFit$par[['alpha']], mu1=bestFit$par['mu1'], sigma1=bestFit$par['sigma1'], mu2=bestFit$par['mu2'], sigma2=bestFit$par['sigma2'])
     r2_single <- 1-(sse_single/sst)
     r2_double <- 1-(sse/sst)
     #      print(paste0('mu1Limits = ', mu1Limits))
     #      print(paste0('sigma1Limits = ', sigma1Limits))
     #      print(paste0('mu2Limits = ', mu2Limits))
     #      print(paste0('sigma2Limits = ', sigma2Limits))
     return(list(par=c(alpha=bestFit$par[['alpha']], mu1=bestFit$par[['mu1']], sigma1=bestFit$par[['sigma1']], mu2=bestFit$par[['mu2']], sigma2=bestFit$par[['sigma2']]), fit=bestFit, sse=sse, sst=sst, r2_single=r2_single, r2_double=r2_double))
}

#' fitLogNormGS
#'
#' @param x default
#' @param y default
#' @param mu default
#' @param sigma default
#' @param alphaGuess default
#' @param mu2Guess default
#' @param sigma2Guess default
#' @param method default
#' @param cores default
#'
#' @export
fitLogNormGS <- function(x, y, mu=lseq(0.0001, 0.5, 50), sigma=seq(0.01,1.51,0.01), alphaGuess=0.9, mu2Guess=-1, sigma2Guess=-1, method='L-BFGS-B', cores=1)
{
     guessGS <- getGuessGS(x, y, mu=mu, sigma=sigma, cores=cores)

     limitFlag <- guessGS[['limitFlag']]
     guessGS$limitFlag <- NULL # remove from list so it doesn't mess up optim function.

     bestFit1 <- fitLogNorm1(x, y, guess=guessGS, method=method)

     guess <- bestFit1$par

     if(mu2Guess < 0)
     {
          mu2Guess = (1+guess[['sigma']])*guess[['mu']]
     }
     if(sigma2Guess < 0)
     {
          sigma2Guess = guess[['sigma']]
     }

     bestFit2 <- fitLogNorm2(x, y, guess=guess, alpha=alphaGuess, mu2Guess=mu2Guess, sigma2Guess=sigma2Guess, method=method)


     return(list(par1=guess, par2=bestFit2$par, r2_single=bestFit2$r2_single, r2_double=bestFit2$r2_double, logNormFit1=bestFit1, logNormFit2=bestFit2, limitFlag=limitFlag))
}

# temp45 <- fitLogNormGS(x, y, cores=4)
# temp45$par1
# temp45$par2
# temp45$r2_single
# temp45$r2_double
# plot(x, y, log='x')
# lines(x, do.call(logNorm, c(list(x=x), temp45$par1)), lwd=3, col='red')
# lines(x, do.call(logNorm2, c(list(x=x), temp45$par2)), lwd=3, col='blue')

# duh <- trackList$getPercentAdhered(velocityThreshold=1)
# tau <- getShearStress(f=getFrequencies(duh$time)$f, pixelAmplitude = 180, mu=0.0052)
#
# plot(tau, duh$percentAdhered/100, log='x')
# bestFit1 <- fitLogNorm(tau, duh$percentAdhered/100, thresh=0.1, SANN=TRUE)
# lines(tau, do.call(logNorm, c(list(x=tau), bestFit1$par)), lwd=3, col='red')
# bestFit2 <- fitLogNorm2(tau, duh$percentAdhered/100, guess=guess, alpha=0.1, newMu=0.01, newSigma=1, SANN=FALSE)
# lines(tau, do.call(logNorm2, c(list(x=tau), bestFit2$par)), lwd=3, col='blue')
#
# # duh2 <- smooth.spline(x=tau, y=duh$percentAdhered/100, spar=0.5)
# #
# # lines(seq(0.002,0.1, 0.002), logNorm(x=seq(0.002,0.1, 0.002), mu=0.01/2, sigma=0.3))
# bestFit
#
# mc.control <- list(mc.cores=4)
# alpha <- seq(0,1,length.out=11)
# mu <- lseq(0.0001,0.5,50)
# sigma <- seq(0.1,5,0.1)
# levels <- list(mu=mu, sigma=sigma)
# levels2 <- list(alpha=alpha, mu1=mu, sigma1=sigma, mu2=mu, sigma2=sigma)
# res <- gridSearch(sseLogNormGS, levels, dataX=tau, dataY=duh$percentAdhered/100, method='multicore', mc.control=mc.control)
# lines(tau, logNorm(x=tau, mu=res$minlevels[1], sigma=res$minlevels[2]), lwd=3, col='red')
# res2 <- gridSearch(sseLogNormGS, levels, dataX=tau, dataY=duh$percentAdhered/100, method='multicore', mc.control=mc.control)
# guess <- c(mu=res$minlevels[1], sigma=res$minlevels[2])
# res <- gridSearch(sseLogNorm, levels, mc.control = mc.control)
#
# temp <- function(x, mu, sigma){mu*x + sigma*x^2}
