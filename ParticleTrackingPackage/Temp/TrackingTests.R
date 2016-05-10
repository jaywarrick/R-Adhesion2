library(pracma)
library(GraphAlignment)

#' @importFrom pracma isempty
#' @importFrom GraphAlignment LinearAssignment
NULL

#' plotAssignments
#'
#' @param pointSet0 PointSet object respresenting the points to link from
#' @param pointSet1 PointSet object representing the point to link to
#' @param digits numeric value
#' @param maxDist numeric value
#' @param direction numeric vector
#' @param perpendicularPenaltyFactor numeric value
#' @param uniformityDistThresh numeric value
#' @param silent logical value
#'
#' @export
plotAssignments <- function(pointSet0, pointSet1, digits=1, maxDist=0.5, direction=c(1,0,0), perpendicularPenaltyFactor=1, uniformityDistThresh=-1, silent=FALSE)
{
     t0 <- pointSet0$pts
     t1 <- pointSet1$pts
     results <- directionalLinearAssignment(pointSet0=pointSet0, pointSet1=pointSet1, maxDist=maxDist, direction=direction, perpendicularPenaltyFactor=perpendicularPenaltyFactor, uniformityDistThresh=uniformityDistThresh, digits=digits)
     px <- results$links$id1
     if(!silent)
     {
          plot(t0$x, t0$y, type='p', pch=20, col='red', xlab='X', ylab='Y')
          points(t1$x, t1$y, type='p', pch=20, col='blue')
          for(i in 1:length(px))
          {
               lines(c(t0$x[px[i]], t1$x[i]), c(t0$y[px[i]], t1$y[i]))
          }
     }

     # Calculate accuracy
     wrongs <- 0
     t0i <- row.names(t0)
     t1i <- row.names(t1)
     pxExpected <- c()
     for(i in 1:length(t1i))
     {
          t1match <- which(t0i == t1i[i])
          if(!isempty(t1match))
          {
               pxExpected <- c(pxExpected, t1match[1])
          }
          else
          {
               pxExpected <- c(pxExpected, NaN)
          }
     }
     total <- numel(pxExpected)
     iToCompare <- which(!is.nan(pxExpected))
     rights <- sum(px[iToCompare] == pxExpected[iToCompare])
     wrongs <- total-rights
     return(list(px=px, errorRate=100*wrongs/total))
}

#' plotAssignments2
#'
#' @param pointSet0 PointSet object respresenting the points to link from
#' @param pointSet1 PointSet object representing the point to link to
#' @param assignmentResults list with a 'links' item that is a data.frame with columns 'id0', and 'id1' indicating which id0's link to which id1's from pointSet0 to pointSet1
#'
#' @export
plotAssignments2 <- function(pointSet0, pointSet1, assignmentResults)
{
     pointSet0$plotPointSet(cex=1, col='red')
     pointSet1$plotPointSet(col='blue', cex=1, add=TRUE)
     if(nrow(assignmentResults$links) > 0)
     {
          for(i in 1:nrow(assignmentResults$links))
          {
               id0 <- assignmentResults$links$id0[i]
               id1 <- assignmentResults$links$id1[i]
               # print(paste("linking ", id0, " to ", id1))
               arrows(x0=pointSet0$getPoint(id=id0)$x, y0=pointSet0$getPoint(id=id0)$y, x1=pointSet1$getPoint(id=id1)$x, y1=pointSet1$getPoint(id=id1)$y, length = 0.1)
          }
     }
}

#' testSimpleTracking
#'
#' Simple test script that generates two point sets and links them
#' forward and backward. One of them has one more point than the other.
#'
#' @export
testSimpleTracking <- function()
{

     # Test forward and backward tracking with randomly ordered points with different numbers of points
     set.seed(1234)
     pts <- getPoints(n=10, movement=c(0,0,0), direction=c(1,0,0), parallelNoise=0.1/4, perpendicularNoise=((0.1/4)/10))
     pts$t1 <- rbind(pts$t1, data.frame(x=0,y=0,z=0))
     pts$t1 <- pts$t1[sample(nrow(pts$t1)),]
     pointSet0 <- new('PointSet')
     pointSet0$initializeWithDataFrame(frame=1, data=pts$t0[,c('x','y')])
     pointSet1 <- new('PointSet')
     pointSet1$initializeWithDataFrame(frame=2, data=pts$t1[,c('x','y')])
     resultsForward <- directionalLinearAssignment(pointSet0=pointSet0, pointSet1=pointSet1, maxDist=0.1, direction=c(1,0), perpendicularPenaltyFactor=10, uniformityDistThresh=0, digits=4)
     plotAssignments2(pointSet0, pointSet1, resultsForward)
     print(paste0("Checking Forward Test. Result = ", min(resultsForward$px == c(2, 10, 13, 9, 1, 4, 8, 7, 6, 5, 3, 19, 14, 20, 16, 15, 12, 11, 21, 18, 17))))
     resultsBackward <- directionalLinearAssignment(pointSet0=pointSet1, pointSet1=pointSet0, maxDist=0.1, direction=c(1,0), perpendicularPenaltyFactor=10, uniformityDistThresh=0, digits=4)
     plotAssignments2(pointSet1, pointSet0, resultsBackward)
     print(paste0("Checking Backward Test. Result = ", min(resultsBackward$px == c(5, 1, 11, 6, 10, 9, 8, 7, 4, 2, 13, 21, 3, 18, 19, 15, 12, 16, 17, 14, 20))))

}

#' testLAPAccuracy_directionality
#'
#' @param min numeric value
#' @param max numeric value
#' @param i numeric value
#' @param perpendicularPenaltyFactor numeric value
#'
#' @export
testLAPAccuracy_directionality <- function(min=10, max=1000, i=5, perpendicularPenaltyFactor=10)
{
     ns <- seq(min, max, length.out=i)
     accuracy1 <- c()
     accuracy2 <- c()
     for(n in ns)
     {
          pts <- getPoints(n=n, movement=c(0,0,0), direction=c(1,0,0), parallelNoise=0.1/4, perpendicularNoise=((0.1/4)/perpendicularPenaltyFactor))
          pointSet0 <- new('PointSet')
          pointSet0$initializeWithDataFrame(frame=1, data=pts$t0[,c('x','y')])
          pointSet1 <- new('PointSet')
          pointSet1$initializeWithDataFrame(frame=2, data=pts$t1[sample(nrow(pts$t1)),][,c('x','y')])
          temp1 <- plotAssignments(pointSet0=pointSet0, pointSet1=pointSet1, digits=1, maxDist=0.4, direction=c(1,0), perpendicularPenaltyFactor=1)
          temp2 <- plotAssignments(pointSet0=pointSet0, pointSet1=pointSet1, digits=1, maxDist=0.4, direction=c(1,0), perpendicularPenaltyFactor=perpendicularPenaltyFactor)
          ptaccuracy1 <- c(accuracy1, temp1$errorRate)
          accuracy2 <- c(accuracy2, temp2$errorRate)
          print(accuracy1)
          print(accuracy2)
     }
     duh2 <- data.frame(n=ns, ran=accuracy1, dir=accuracy2)
     plot(duh2$n, duh2$ran, type='l', col='black', main='Tracking Errors', xlab='Cell Density [cells/area]', ylab='Tracking Error Rate [%]')
     lines(duh2$n, duh2$dir, type='l', col='blue')
     legend('topleft', c('isotropic','anisotropic'), lty=c(1,1), col=c('black','blue'), inset = .05, bty='n')
     return(duh2)
}

# ccc # ccc testLAPAccuracy_outOfFrame
# ccc # ccc
# ccc # ccc @param n numeric value
# ccc # ccc @param mags numeric vector
# ccc # ccc @param perpendicularPenaltyFactor numeric value
# ccc # ccc
# ccc # ccc @export
# ccc testLAPAccuracy_outOfFrame <- function(n=10, mags=seq(0,0.1,length.out=3), perpendicularPenaltyFactor=10)
# ccc {
# ccc      accuracy1 <- c()
# ccc      accuracy2 <- c()
# ccc      oof <- c() # Out Of Frame (oof)
# ccc      for(m in mags)
# ccc      {
# ccc           # Create pts to track
# ccc           pts <- getPoints(n=n, movement=c(m,0,0), direction=c(1,0,0), parallelNoise=0.1/4, perpendicularNoise=((0.1/4)/perpendicularPenaltyFactor))
# ccc           t1OutOfFrame <- subset(pts$t1, get('x') > 1)
# ccc           t1InFrame <- subset(pts$t1, get('x') <= 1)
# ccc           newPoints <- pts
# ccc           newPoints$t1 <- t1InFrame
# ccc
# ccc           # Track them
# ccc           temp1 <- plotAssignments(pointSet0=pts$t0, pointSet1=pts$t1, digits=1, maxDist=m+2*0.1, direction=c(1,0,0), perpendicularPenaltyFactor=perpendicularPenaltyFactor)
# ccc           temp2 <- plotAssignments(pointSet0=pts$t0, pointSet1=pts$t1, digits=1, maxDist=m+2*0.1, direction=c(1,0,0), perpendicularPenaltyFactor=perpendicularPenaltyFactor)
# ccc
# ccc           # Store the results
# ccc           accuracy1 <- c(accuracy1, temp1$errorRate)
# ccc           accuracy2 <- c(accuracy2, temp2$errorRate)
# ccc           oof <- c(oof, nrow(t1OutOfFrame))
# ccc           print(accuracy1)
# ccc           print(accuracy2)
# ccc      }
# ccc
# ccc      poof <- 100*(oof/n) # Percent Out Of Frame (poof)
# ccc      duh2 <- data.frame(poof=poof, allIn=accuracy1, someOut=accuracy2)
# ccc      plot(poof, duh2$someOut, ylim=c(0,max(duh2$someOut)), type='l', col='black', main='Tracking Errors', xlab='Percent Out of Frame [%]', ylab='Tracking Error Rate [%]')
# ccc      lines(poof, duh2$allIn, type='l', col='blue')
# ccc      legend('topleft', c('some out of frame','all in frame'), lty=c(1,1), col=c('black','blue'), inset = .05, bty='n')
# ccc      return(duh2)
# ccc }
# ccc
# ccc # ccc testLAPAccuracy_outOfFrame2
# ccc # ccc
# ccc # ccc @param n numeric value
# ccc # ccc @param mags numeric vector
# ccc # ccc @param perpendicularPenaltyFactor numeric value
# ccc # ccc
# ccc # ccc @export
# ccc testLAPAccuracy_outOfFrame2 <- function(n=10, mags=seq(0,0.1,length.out=6), perpendicularPenaltyFactor=10)
# ccc {
# ccc      accuracy1 <- c()
# ccc      accuracy2 <- c()
# ccc      oof <- c() # Out Of Frame (oof)
# ccc      for(m in mags)
# ccc      {
# ccc           # Create pts to track
# ccc           pts <- getPoints(n=n, movement=c(m,0,0), direction=c(1,0,0), parallelNoise=0.1/4, perpendicularNoise=((0.1/4)/perpendicularPenaltyFactor))
# ccc           t0OutOfFrame <- subset(pts$t0, get('x') < 0.1)
# ccc           t0InFrame <- subset(pts$t0, get('x') >= 0.1)
# ccc           t1OutOfFrame <- subset(pts$t1, get('x') > 1)
# ccc           t1InFrame <- subset(pts$t1, get('x') <= 1)
# ccc           newPoints <- pts
# ccc           newPoints$t0 <- t0InFrame
# ccc           newPoints$t1 <- t1InFrame
# ccc           # Track them
# ccc           temp1 <- plotAssignments(pts=pts, digits=1, maxDist=m+2*0.1, direction=c(1,0,0), perpendicularPenaltyFactor=10)
# ccc           temp2 <- plotAssignments(pts=newPoints, digits=1, maxDist=m+2*0.1, direction=c(1,0,0), perpendicularPenaltyFactor=10)
# ccc
# ccc           # Store the results
# ccc           accuracy1 <- c(accuracy1, temp1$errorRate)
# ccc           accuracy2 <- c(accuracy2, temp2$errorRate)
# ccc           oof <- c(oof, nrow(t1OutOfFrame))
# ccc           print(accuracy1)
# ccc           print(accuracy2)
# ccc      }
# ccc
# ccc      poof <- 100*(oof/n) # Percent Out Of Frame (poof)
# ccc      duh2 <- data.frame(poof=poof, allIn=accuracy1, someOut=accuracy2)
# ccc      plot(poof, duh2$someOut, ylim=c(0,max(duh2$someOut)), type='l', col='black', main='Tracking Errors', xlab='Percent Out of Frame [%]', ylab='Tracking Error Rate [%]')
# ccc      lines(poof, duh2$allIn, type='l', col='blue')
# ccc      legend('topleft', c('some out of frame','all in frame'), lty=c(1,1), col=c('black','blue'), inset = .05, bty='n')
# ccc      return(duh2)
# ccc }
# ccc
# ccc # ccc testLAPAccuracy_outOfFrame3
# ccc # ccc
# ccc # ccc @param n numeric value
# ccc # ccc @param mags numeric vector
# ccc # ccc @param perpendicularPenaltyFactor numeric value
# ccc # ccc @param silent logical value
# ccc # ccc
# ccc # ccc @export
# ccc testLAPAccuracy_outOfFrame3 <- function(n=500, mags=seq(0.1,0.1,length.out=30), perpendicularPenaltyFactor=10, silent=TRUE)
# ccc {
# ccc      accuracy1 <- c()
# ccc      accuracy2 <- c()
# ccc      accuracy3 <- c()
# ccc      accuracy4 <- c()
# ccc      accuracy5 <- c()
# ccc      accuracy6 <- c()
# ccc      poof <- c() # Percent Out Of Frame (poof)
# ccc      for(m in mags)
# ccc      {
# ccc           # Create pts to track
# ccc           #           perpendicularPenaltyFactor <- 10
# ccc           #           m <- 0.1
# ccc           #           n <- 10
# ccc           pts <- getPoints(n=n, xlim=c(-0.5,1.5), ylim=c(0.5,1.5), zlim=c(0,0), movement=c(m,0,0), direction=c(1,0,0), parallelNoise=0.1/4, perpendicularNoise=((0.1/4)/perpendicularPenaltyFactor))
# ccc           t1InFrameRows <- which(pts$t1$x >= 0 & pts$t1$x <= 1, pts$t1$y >= 0 & pts$t1$y <= 1)
# ccc           t0InFrameRows <- which(pts$t0$x >= 0 & pts$t0$x <= 1, pts$t0$y >= 0 & pts$t0$y <= 1)
# ccc           completeData <- intersect(t1InFrameRows, t0InFrameRows)
# ccc           allData <- union(t1InFrameRows, t0InFrameRows)
# ccc
# ccc           fractionOutOfFrame <- (length(allData) - length(completeData))/(length(allData))
# ccc           completePoints <- list(t0=pts$t0[t0InFrameRows,], t1=pts$t1[t0InFrameRows,])
# ccc           incompletePoints <- list(t0=pts$t0[t0InFrameRows,], t1=pts$t1[t1InFrameRows,])
# ccc           # Track them
# ccc
# ccc           temp1 <- plotAssignments(pts=completePoints, digits=1, maxDist=m+2*0.1, direction=c(1,0,0), perpendicularPenaltyFactor=1, uniformityDistThresh=-1, silent=silent)
# ccc           temp2 <- plotAssignments(pts=completePoints, digits=1, maxDist=m+2*0.1, direction=c(1,0,0), perpendicularPenaltyFactor=10, uniformityDistThresh=-1, silent=silent)
# ccc           temp3 <- plotAssignments(pts=completePoints, digits=1, maxDist=m+2*0.1, direction=c(1,0,0), perpendicularPenaltyFactor=10, uniformityDistThresh=3*0.1/4, silent=silent)
# ccc           temp4 <- plotAssignments(pts=incompletePoints, digits=1, maxDist=m+2*0.1, direction=c(1,0,0), perpendicularPenaltyFactor=1, uniformityDistThresh=0, silent=silent)
# ccc           temp5 <- plotAssignments(pts=incompletePoints, digits=1, maxDist=m+2*0.1, direction=c(1,0,0), perpendicularPenaltyFactor=10, uniformityDistThresh=-1, silent=silent)
# ccc           temp6 <- plotAssignments(pts=incompletePoints, digits=1, maxDist=m+2*0.1, direction=c(1,0,0), perpendicularPenaltyFactor=10, uniformityDistThresh=3*0.1/4, silent=silent)
# ccc
# ccc           # Store the results
# ccc           accuracy1 <- c(accuracy1, temp1$errorRate)
# ccc           accuracy2 <- c(accuracy2, temp2$errorRate)
# ccc           accuracy3 <- c(accuracy3, temp3$errorRate)
# ccc           accuracy4 <- c(accuracy4, temp4$errorRate)
# ccc           accuracy5 <- c(accuracy5, temp5$errorRate)
# ccc           accuracy6 <- c(accuracy6, temp6$errorRate)
# ccc           poof <- c(poof, 100*fractionOutOfFrame)
# ccc           print(accuracy1)
# ccc           print(accuracy2)
# ccc           print(accuracy3)
# ccc           print(accuracy4)
# ccc           print(accuracy5)
# ccc           print(accuracy6)
# ccc      }
# ccc      daMax <- max(accuracy1,accuracy2,accuracy3,accuracy4,accuracy5,accuracy6)
# ccc      duh2 <- data.frame(poof=poof, error1=accuracy1, error2=accuracy2, error3=accuracy3, error4=accuracy4, error5=accuracy5, error6=accuracy6)
# ccc      duh2 <- duh2[order(duh2$poof),]
# ccc      plot(poof, duh2$error1, ylim=c(0,daMax), type='l', lty=1, col='blue', main='Tracking Errors', xlab='Percent Out of Frame [%]', ylab='Tracking Error Rate [%]')
# ccc      lines(poof, duh2$error2, lty=3, col='blue')
# ccc      lines(poof, duh2$error3, lty=5, col='blue')
# ccc      lines(poof, duh2$error4, lty=1, col='black')
# ccc      lines(poof, duh2$error5, lty=3, col='black')
# ccc      lines(poof, duh2$error6, lty=5, col='black')
# ccc      legend('topleft', c('Complete, Isotropic, Unconstrained', 'Complete, Anisotropic, Unconstrained', 'Complete, Anisotropic, Constrained', 'Incomplete, Isotropic, Unconstrained', 'Incomplete, Anisotropic, Unconstrained', 'Incomplete, Anisotropic, Constrained'), lty=c(1,3,5,1,3,5), col=c('blue','blue','blue','black','black','black'), inset = .05, bty='n', cex=0.5)
# ccc      return(duh2)
# ccc }
# ccc
# ccc # ccc comparison
# ccc # ccc
# ccc # ccc @param mags numeric vector
# ccc # ccc @param cells numeric vector
# ccc # ccc @param reps numeric value
# ccc # ccc
# ccc # ccc @export
# ccc comparison <- function(mags=seq(0,0.1,0.025), cells=seq(50,350,50), reps=50)
# ccc {
# ccc      results <- data.frame(mag=c(), cells=c(), poofMean=c(), error1Mean=c(), error2Mean=c(), error3Mean=c(), error4Mean=c(), error5Mean=c(), error6Mean=c(), poofMad=c(), error1Mad=c(), error2Mad=c(), error3Mad=c(), error4Mad=c(), error5Mad=c(), error6Mad=c())
# ccc      for(m in mags)
# ccc      {
# ccc           for(c in cells)
# ccc           {
# ccc                duh <- testLAPAccuracy_outOfFrame3(n=c, mags <- seq(m,m,length.out=reps))
# ccc                means <- colMedians(as.matrix(duh))
# ccc                names(means) <- c('poofMean', 'error1Mean', 'error2Mean', 'error3Mean', 'error4Mean', 'error5Mean', 'error6Mean')
# ccc                mads <- colMads(as.matrix(duh))
# ccc                names(mads) <- c('poofMad', 'error1Mad', 'error2Mad', 'error3Mad', 'error4Mad', 'error5Mad', 'error6Mad')
# ccc                results <- rbind(results, data.frame(mag=m, cells=c, as.list(means), as.list(mads)))
# ccc           }
# ccc      }
# ccc      return(results)
# ccc }
# ccc
# ccc # ccc plotMethod
# ccc # ccc
# ccc # ccc @param results data.frame
# ccc # ccc @param name character value
# ccc # ccc @param normalize logical value
# ccc # ccc
# ccc # ccc @export
# ccc plotMethod <- function(results, name, normalize=FALSE)
# ccc {
# ccc      xlim <- c(0, max(results$poofMean))
# ccc      ylim <- c(0, max(results[,name]))
# ccc      if(normalize)
# ccc      {
# ccc           ylim <- c(0, max(results[,c('error1Mean','error2Mean','error3Mean','error4Mean','error5Mean','error6Mean')]))
# ccc      }
# ccc
# ccc      mags <- unique(results$mag)
# ccc      cells <- unique(results$cells)
# ccc      plot(x=c(), y=c(), xlim=xlim, ylim=ylim, type='l', xlab='Points Out Of Frame [%]', ylab='Error Rate [%]')
# ccc      myCount <- 1
# ccc      for(c in cells)
# ccc      {
# ccc           data <- subset(results, cells==c)
# ccc           lines(data$poofMean, data[,name], lty=myCount)
# ccc           myCount <- myCount + 1
# ccc      }
# ccc }

#' getPoints
#'
#' @param n numeric value
#' @param xlim numeric vector
#' @param ylim numeric vector
#' @param zlim numeric vector
#' @param movement numeric vector
#' @param direction numeric vector
#' @param parallelNoise numeric value
#' @param perpendicularNoise numeric value
#'
#' @export
getPoints <- function(n=100, xlim=c(0,1), ylim=c(0,1), zlim=c(0,1), movement=c(0.1,0,0), direction=c(1,0,0), parallelNoise=0.1/4, perpendicularNoise=parallelNoise/10)
{
     print('Generating pts')
     # Define motions
     vectors <- createMovements(n=n, movement=movement, direction=direction, parallelNoise=parallelNoise, perpendicularNoise=perpendicularNoise)

     # Define locations
     a0 <- runif(n, min=xlim[1], max=xlim[2])
     b0 <- runif(n, min=ylim[1], max=ylim[2])
     c0 <- runif(n, min=zlim[1], max=zlim[2])

     t0 <- data.frame(x=a0, y=b0, z=c0)
     t1 <- t0 + vectors

     # Eliminate z for now by setting all z's to 0
     t0$z <- 0
     t1$z <- 0

     # Return
     return(list(t0=t0, t1=t1))
}

#' getRotationMatrix
#'
#' @param v1 numeric vector
#' @param v2 numeric vector
#'
#' @export
getRotationMatrix <- function(v1=c(1,0,0), v2=c(1,1,1))
{
     mag1 <- sqrt(sum(v1^2))
     u1 <- v1/mag1
     mag2 <- sqrt(sum(v2^2))
     u2 <- v2/mag2

     if(sum(u1==u2)==3)
     {
          return(diag(1,3,3))
     }

     v <- cross(u1, u2)
     sAngle <- sqrt(sum(v*v))
     cAngle <- sum(u1*u2)
     vx <- matrix(c(0, v[3], -v[2], -v[3], 0, v[1], v[2], -v[1], 0), 3, 3)

     R <- diag(1,3,3) + vx + (vx%*%vx)*((1-cAngle)/(sAngle^2))

     return(R)
}

#' createMovements
#'
#' @param n numeric value
#' @param movement numeric value
#' @param direction numeric vector
#' @param parallelNoise numeric value
#' @param perpendicularNoise numeric value
#'
#' @export
createMovements <- function(n, movement, direction=movement, parallelNoise, perpendicularNoise)
{
     # We use this method to create "movement" along the x-direction
     # We will eventually rotate them to be along the direction of the mo
     r <- sqrt(sum((movement^2)))
     xNoise <- rnorm(n, 0, parallelNoise)
     yNoise <- rnorm(n, 0, perpendicularNoise)
     zNoise <- rnorm(n, 0, perpendicularNoise)
     x <- xNoise + r
     y <- yNoise
     z <- zNoise

     temp <- rbind(x,y,z)
     R <- getRotationMatrix(v1=c(1,0,0), v2=direction)

     # Rotate all the movement to be along the direction supplied in movement
     results <- numeric(0)
     for(i in 1:length(x))
     {
          results <- cbind(results, R%*%temp[,i])
     }

     return(data.frame(x=as.numeric(results[1,]), y=as.numeric(results[2,]), z=as.numeric(results[3,])))
}

#' getPoints
#'
#' @param n numeric value
#' @param xlim numeric vector
#' @param ylim numeric vector
#' @param zlim numeric vector
#' @param movement numeric vector
#' @param direction numeric vector
#' @param parallelNoise numeric value
#' @param perpendicularNoise numeric value
#'
#' @export
getPointSets <- function(n=100, xlim=c(0,1), ylim=c(0,1), zlim=c(0,1), movement=c(0.1,0,0), direction=c(1,0,0), parallelNoise=0.1/4, perpendicularNoise=parallelNoise/10)
{
     print('Generating pts')
     # Define motions
     vectors <- createMovements(n=n, movement=movement, direction=direction, parallelNoise=parallelNoise, perpendicularNoise=perpendicularNoise)

     # Define locations
     a0 <- runif(n, min=xlim[1], max=xlim[2])
     b0 <- runif(n, min=ylim[1], max=ylim[2])
     c0 <- runif(n, min=zlim[1], max=zlim[2])

     t0 <- data.frame(x=a0, y=b0, z=c0)
     t1 <- t0 + vectors

     # Eliminate z for now by setting all z's to 0
     t0$z <- 0
     t1$z <- 0

     # Make the pointSets
     ps0 <- new('PointSet')
     ps0$initializeWithDataFrame(frame=1, data=t0)
     ps1 <- new('PointSet')
     ps1$initializeWithDataFrame(frame=2, data=t1)

     # Return
     return(list(t0=ps0, t1=ps1))
}
