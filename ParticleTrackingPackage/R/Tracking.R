
# Installe the GraphAlignment package from bioconductor if necessary
if ( !require("GraphAlignment") )
{
     source("http://www.bioconductor.org/biocLite.R")
     biocLite("GraphAlignment", suppressUpdates = T)
}

# duh <- getUL(t0=t0, t1=t1, linkCostFunction=directionalCost, theta=0, penaltyFactors=c(0,0,1))

##### Cost Functions #####

getBlockingCost <- function(maxDist)
{
     return(10*(maxDist*maxDist))
}

# REMEMBER TO CHECK AND SEE WHEN THIS IS GETTING CALLED
# IF AFTER WE SET THRESHOLD, THEN THIS CHANGES WHAT WE THINK
# ABOUT THE CUTOFF VALUE AND WE CAN POTENTIALLY CHANGE THE FUNCTION
# TO USE MAXDIST
getCutoff <- function(UL, blockingCost, factorOverMax=10)
{
     potentialMatches <- UL[UL < blockingCost]
     if(isempty(potentialMatches))
     {
          return(blockingCost/10)
     }
     else
     {
          return(factorOverMax*max(UL[UL < blockingCost]))
     }
}

getLinkMatrix <- function(points, maxDist=15, direction=c(1,0,0), directionality=10)
{
     t0 <- points$t0
     t1 <- points$t1

     # Check inputs
     if(length(direction) < 3 || ncol(t0) < 3 || ncol(t1) < 3)
     {
          print("Invalid arguments to getLinkMatrix function. Returning NULL.")
          return(NULL)
     }

     # make sure direction is a unit vector by dividing all components by its magnitude
     direction <- direction/(sqrt(sum(direction^2)))

     # need to calculate distance for all combinations of points
     # To do this we replicate the t0 and t1 matrices appropriately and do matrix math
     # Note that t0 and t1 can have different lengths
     nt0 <- nrow(t0)
     nt1 <- nrow(t1)
     t0Rep <- matrix(rep(t(t0), nt1), ncol=ncol(t0), byrow=TRUE)
     t1Seq <- rep(1:nrow(t1), each=nt0)
     t1Rep <- t1[t1Seq,]

     # Calculate movement vectors for each combination
     movement <- as.matrix(t1Rep-t0Rep)
     parallelMovement <- movement%*%direction # dot product of direction with each row of the movement matrix
     movementSign <- sign(parallelMovement)
     perpendicularMovement2 <- rowSums(movement*movement) - parallelMovement*parallelMovement

     # Calculate cost
     cost <- parallelMovement*parallelMovement + ((directionality*directionality)*perpendicularMovement2) # Penalize perpendicular movement in proportion to the expected ratio of motion expected parallel vs perpendicular to the specified direction
     cost[cost >= (maxDist*maxDist)] <- getBlockingCost(maxDist)
     signedCost <- cost*movementSign
     cost <- matrix(cost, nt0, nt1)
     signedCost <- matrix(signedCost, nt0, nt1)

     return(list(cost=cost, signedCost=signedCost))
}

DirectionalLinearAssignment <- function(points, maxDist=10, direction=c(1,0,0), directionality=10, uniformityDistThresh=-1, digits=4)
{
     # Do first round of tracking as first guess
     #cat("Getting initial cost matrix\n")
     linkMatrix <- getLinkMatrix(points, maxDist=maxDist, direction=direction, directionality=directionality)
     cost <- getCostMatrix(UL=linkMatrix$cost, digits=digits, maxDist=maxDist)

     #cat("Linking (round 1)\n")
     px <- LinearAssignment(cost)

     # A negative number indicates that we should not apply a uniformity constraint
     if(uniformityDistThresh >= 0)
     {
          uniformityCostThresh <- (directionality*uniformityDistThresh)^2 + uniformityDistThresh^2

          #cat("Determining directionality\n")
          diagonal <- 1:length(px)
          # First find valid rows and cols
          nrows <- nrow(points$t0)
          ncols <- nrow(points$t1)
          # Only look at pxs that reference links between t0 and t1 and not broken links
          validPx <- which(px <= nrows)
          validPx <- validPx[validPx <= ncols] # only care about linked points (this selects valid columns)
          linkedCosts <- numeric(length(validPx))
          for(i in 1:length(validPx))
          {
               linkedCosts[i] <- linkMatrix$signedCost[px[validPx[i]],validPx[i]]
          }

          linkedSign <- -1*sign(sum(sign(linkedCosts))) # determine which sign to penalize/remove (i.e., replace with blocking cost)
          if(!is.na(linkedSign)) # this catches condition where there are essentially all blocking costs (i.e., no valid links)
          {
               if(linkedSign > 0)
               {
                    linkMatrix$cost[linkMatrix$signedCost > uniformityCostThresh] <- getBlockingCost(maxDist)
               }
               else
               {
                    linkMatrix$cost[linkMatrix$signedCost < -1*uniformityCostThresh] <- getBlockingCost(maxDist)
               }
          }

          cost <- getCostMatrix(UL=linkMatrix$cost, digits=digits, maxDist=maxDist)

          #cat("Linking (round 2)\n\n")
          px <- LinearAssignment(cost)
     }

     return(list(px=px, cost=cost))
}

##### Matrix Assembly #####

getURorLL <- function(n, blockingCost, cutoff)
{
     UR <- matrix(blockingCost, n, n);

     # Set the cutoff along the diagonal (top left to bottom right)
     for (i in 1:n)
     {
          UR[i, i] <- cutoff
     }

     return(UR);
}

getLR <- function(UL, blockingCost)
{
     # minVal <- min(UL)
     cutoff <- getCutoff(UL=UL, blockingCost=blockingCost)
     LR = t(UL)
     r <- nrow(LR)
     c <- ncol(LR)
     for (i in 1:r)
     {
          for (j in 1:c)
          {
               if(numel(LR[i, j] < blockingCost) == 0)
               {
                    print('whoa!')
               }
               if (LR[i, j] < blockingCost)
               {
                    LR[i, j] <- cutoff
               }
          }
     }
     return(LR)
}

getCostMatrix <- function(UL, digits=4, maxDist)
{
     if(maxDist == 0)
     {
          print('Warning. Likely error because maxDist is set to 0!')
     }
     blockingCost <- getBlockingCost(maxDist)
     cutoff <- getCutoff(UL=UL, blockingCost=blockingCost)
     UR <- getURorLL(n=nrow(UL), blockingCost=blockingCost, cutoff=cutoff)
     LL <- getURorLL(n=ncol(UL), blockingCost=blockingCost, cutoff=cutoff)
     LR <- getLR(UL=UL, blockingCost=blockingCost)
     ret <- rbind(UL,LL)
     temp <- rbind(UR,LR)
     ret <- cbind(ret,temp)
     ret <- round((10^digits)*ret) # turn it into an integer problem for the LinearAssignment function of the GraphAlignment package
     return(ret)
}

##### Visualization #####

plotAssignments <- function(points, digits=4, maxDist=0.5, direction=c(1,0,0), directionality=1, uniformityDistThresh=-1, silent=FALSE)
{
     t0 <- points$t0
     t1 <- points$t1
     results <- DirectionalLinearAssignment(points=points, maxDist=maxDist, direction=direction, directionality=directionality, uniformityDistThresh=uniformityDistThresh, digits=digits)
     px <- results$px
     cost <- results$cost
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
     return(list(px=px, errorRate=100*wrongs/total, costMat=cost))
}


##### Point Generation #####

getPoints <- function(n=100, xlim=c(0,1), ylim=c(0,1), zlim=c(0,1), movement=c(0.1,0,0), direction=c(1,0,0), parallelNoise=0.1/4, perpendicularNoise=parallelNoise/10)
{
     print('Generating points')
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
##### Commented Code #####
# getLinkingCost <- function(p1, p2, maxDist, penaltyFactors)
# {
#      # Do the cost calculation
#      cost <- getCost(p1, p2)
#      cutoff <- maxDist*maxDist
#
#      # Apply the threshold linking distance
#      if (cost > cutoff)
#      {
#           return(getBlockingCost(maxDist))
#      }
#
#      # Apply penalties
#      penalty <- 1;
#      for (i in 1:length(p1))
#      {
#           ndiff <- getNormDiff(p1[i], p2[i])
#           if (is.nan(ndiff))
#           {
#                next
#           }
#           factor <- penaltyFactors[i]
#           penalty <- penalty + factor * 1.5 * ndiff
#      }
#
#      # calculate final total cost (total cost = cost * pentalty^2 = (distance + penalty)^2)
#      cost <- (cost * penalty * penalty)
#
#      # Apply the threshold linking distance
#      if (cost > cutoff)
#      {
#           return(getBlockingCost(maxDist))
#      }
#
#      # Return the total cost
#      return(cost);
# }
#
# getNormDiff <- function(a, b)
# {
#      if ( a == -b )
#      {
#           return(0)
#      }
#      else
#      {
#           return(abs( a - b ) / ( ( a + b ) / 2 ))
#      }
# }


##### Tests #####
testLAPAccuracy_directionality <- function(min=10, max=1000, i=5, directionality=10)
{
     ns <- seq(min, max, length.out=i)
     accuracy1 <- c()
     accuracy2 <- c()
     for(n in ns)
     {
          points <- getPoints(n=n, movement=c(0,0,0), direction=c(1,0,0), parallelNoise=0.1/4, perpendicularNoise=((0.1/4)/directionality))
          temp1 <- plotAssignments(points=points, digits=4, maxDist=0.4, direction=c(1,0,0), directionality=1)
          temp2 <- plotAssignments(points=points, digits=4, maxDist=0.4, direction=c(1,0,0), directionality=directionality)
          accuracy1 <- c(accuracy1, temp1$errorRate)
          accuracy2 <- c(accuracy2, temp2$errorRate)
          print(accuracy1)
          print(accuracy2)
     }
     duh2 <- data.frame(n=ns, ran=accuracy1, dir=accuracy2)
     plot(duh2$n, duh2$ran, ylim=c(0,max(duh2$someOut)), type='l', col='black', main='Tracking Errors', xlab='Cell Density [cells/area]', ylab='Tracking Error Rate [%]')
     lines(duh2$n, duh2$dir, type='l', col='blue')
     legend('topleft', c('isotropic','anisotropic'), lty=c(1,1), col=c('black','blue'), inset = .05, bty='n')
     return(duh2)
}

testLAPAccuracy_outOfFrame <- function(n=10, mags=seq(0,0.1,length.out=3))
{
     accuracy1 <- c()
     accuracy2 <- c()
     oof <- c() # Out Of Frame (oof)
     for(m in mags)
     {
          # Create points to track
          points <- getPoints(n=n, movement=c(m,0,0), direction=c(1,0,0), parallelNoise=0.1/4, perpendicularNoise=((0.1/4)/directionality))
          t1OutOfFrame <- subset(points$t1, x > 1)
          t1InFrame <- subset(points$t1, x <= 1)
          newPoints <- points
          newPoints$t1 <- t1InFrame

          # Track them
          temp1 <- plotAssignments(points=points, digits=4, maxDist=m+2*0.1, direction=c(1,0,0), directionality=10)
          temp2 <- plotAssignments(points=newPoints, digits=4, maxDist=m+2*0.1, direction=c(1,0,0), directionality=10)

          # Store the results
          accuracy1 <- c(accuracy1, temp1$errorRate)
          accuracy2 <- c(accuracy2, temp2$errorRate)
          oof <- c(oof, nrow(t1OutOfFrame))
          print(accuracy1)
          print(accuracy2)
     }

     poof <- 100*(oof/n) # Percent Out Of Frame (poof)
     duh2 <- data.frame(poof=poof, allIn=accuracy1, someOut=accuracy2)
     plot(poof, duh2$someOut, ylim=c(0,max(duh2$someOut)), type='l', col='black', main='Tracking Errors', xlab='Percent Out of Frame [%]', ylab='Tracking Error Rate [%]')
     lines(poof, duh2$allIn, type='l', col='blue')
     legend('topleft', c('some out of frame','all in frame'), lty=c(1,1), col=c('black','blue'), inset = .05, bty='n')
     return(duh2)
}

testLAPAccuracy_outOfFrame2 <- function(n=10, mags=seq(0,0.1,length.out=6), directionality=10)
{
     accuracy1 <- c()
     accuracy2 <- c()
     oof <- c() # Out Of Frame (oof)
     for(m in mags)
     {
          # Create points to track
          points <- getPoints(n=n, movement=c(m,0,0), direction=c(1,0,0), parallelNoise=0.1/4, perpendicularNoise=((0.1/4)/directionality))
          t0OutOfFrame <- subset(points$t0, x < 0.1)
          t0InFrame <- subset(points$t0, x >= 0.1)
          t1OutOfFrame <- subset(points$t1, x > 1)
          t1InFrame <- subset(points$t1, x <= 1)
          newPoints <- points
          newPoints$t0 <- t0InFrame
          newPoints$t1 <- t1InFrame
          # Track them
          temp1 <- plotAssignments(points=points, digits=4, maxDist=m+2*0.1, direction=c(1,0,0), directionality=10)
          temp2 <- plotAssignments(points=newPoints, digits=4, maxDist=m+2*0.1, direction=c(1,0,0), directionality=10)

          # Store the results
          accuracy1 <- c(accuracy1, temp1$errorRate)
          accuracy2 <- c(accuracy2, temp2$errorRate)
          oof <- c(oof, nrow(t1OutOfFrame))
          print(accuracy1)
          print(accuracy2)
     }

     poof <- 100*(oof/n) # Percent Out Of Frame (poof)
     duh2 <- data.frame(poof=poof, allIn=accuracy1, someOut=accuracy2)
     plot(poof, duh2$someOut, ylim=c(0,max(duh2$someOut)), type='l', col='black', main='Tracking Errors', xlab='Percent Out of Frame [%]', ylab='Tracking Error Rate [%]')
     lines(poof, duh2$allIn, type='l', col='blue')
     legend('topleft', c('some out of frame','all in frame'), lty=c(1,1), col=c('black','blue'), inset = .05, bty='n')
     return(duh2)
}

testLAPAccuracy_outOfFrame3 <- function(n=500, mags=seq(0.1,0.1,length.out=30), directionality=10, silent=TRUE)
{
     accuracy1 <- c()
     accuracy2 <- c()
     accuracy3 <- c()
     accuracy4 <- c()
     accuracy5 <- c()
     accuracy6 <- c()
     poof <- c() # Percent Out Of Frame (poof)
     for(m in mags)
     {
          # Create points to track
          #           directionality <- 10
          #           m <- 0.1
          #           n <- 10
          points <- getPoints(n=n, xlim=c(-0.5,1.5), ylim=c(0.5,1.5), zlim=c(0,0), movement=c(m,0,0), direction=c(1,0,0), parallelNoise=0.1/4, perpendicularNoise=((0.1/4)/directionality))
          t1InFrameRows <- which(points$t1$x >= 0 & points$t1$x <= 1, points$t1$y >= 0 & points$t1$y <= 1)
          t0InFrameRows <- which(points$t0$x >= 0 & points$t0$x <= 1, points$t0$y >= 0 & points$t0$y <= 1)
          completeData <- intersect(t1InFrameRows, t0InFrameRows)
          allData <- union(t1InFrameRows, t0InFrameRows)

          fractionOutOfFrame <- (length(allData) - length(completeData))/(length(allData))
          completePoints <- list(t0=points$t0[t0InFrameRows,], t1=points$t1[t0InFrameRows,])
          incompletePoints <- list(t0=points$t0[t0InFrameRows,], t1=points$t1[t1InFrameRows,])
          # Track them

          temp1 <- plotAssignments(points=completePoints, digits=4, maxDist=m+2*0.1, direction=c(1,0,0), directionality=1, uniformityDistThresh=-1, silent=silent)
          temp2 <- plotAssignments(points=completePoints, digits=4, maxDist=m+2*0.1, direction=c(1,0,0), directionality=10, uniformityDistThresh=-1, silent=silent)
          temp3 <- plotAssignments(points=completePoints, digits=4, maxDist=m+2*0.1, direction=c(1,0,0), directionality=10, uniformityDistThresh=3*0.1/4, silent=silent)
          temp4 <- plotAssignments(points=incompletePoints, digits=4, maxDist=m+2*0.1, direction=c(1,0,0), directionality=1, uniformityDistThresh=0, silent=silent)
          temp5 <- plotAssignments(points=incompletePoints, digits=4, maxDist=m+2*0.1, direction=c(1,0,0), directionality=10, uniformityDistThresh=-1, silent=silent)
          temp6 <- plotAssignments(points=incompletePoints, digits=4, maxDist=m+2*0.1, direction=c(1,0,0), directionality=10, uniformityDistThresh=3*0.1/4, silent=silent)

          # Store the results
          accuracy1 <- c(accuracy1, temp1$errorRate)
          accuracy2 <- c(accuracy2, temp2$errorRate)
          accuracy3 <- c(accuracy3, temp3$errorRate)
          accuracy4 <- c(accuracy4, temp4$errorRate)
          accuracy5 <- c(accuracy5, temp5$errorRate)
          accuracy6 <- c(accuracy6, temp6$errorRate)
          poof <- c(poof, 100*fractionOutOfFrame)
          print(accuracy1)
          print(accuracy2)
          print(accuracy3)
          print(accuracy4)
          print(accuracy5)
          print(accuracy6)
     }
     daMax <- max(accuracy1,accuracy2,accuracy3,accuracy4,accuracy5,accuracy6)
     duh2 <- data.frame(poof=poof, error1=accuracy1, error2=accuracy2, error3=accuracy3, error4=accuracy4, error5=accuracy5, error6=accuracy6)
     duh2 <- duh2[order(duh2$poof),]
     plot(poof, duh2$error1, ylim=c(0,daMax), type='l', lty=1, col='blue', main='Tracking Errors', xlab='Percent Out of Frame [%]', ylab='Tracking Error Rate [%]')
     lines(poof, duh2$error2, lty=3, col='blue')
     lines(poof, duh2$error3, lty=5, col='blue')
     lines(poof, duh2$error4, lty=1, col='black')
     lines(poof, duh2$error5, lty=3, col='black')
     lines(poof, duh2$error6, lty=5, col='black')
     legend('topleft', c('Complete, Isotropic, Unconstrained', 'Complete, Anisotropic, Unconstrained', 'Complete, Anisotropic, Constrained', 'Incomplete, Isotropic, Unconstrained', 'Incomplete, Anisotropic, Unconstrained', 'Incomplete, Anisotropic, Constrained'), lty=c(1,3,5,1,3,5), col=c('blue','blue','blue','black','black','black'), inset = .05, bty='n', cex=0.5)
     return(duh2)
}

comparison <- function(mags=seq(0,0.1,0.025), cells=seq(50,350,50), reps=50)
{
     results <- data.frame(mag=c(), cells=c(), poofMean=c(), error1Mean=c(), error2Mean=c(), error3Mean=c(), error4Mean=c(), error5Mean=c(), error6Mean=c(), poofMad=c(), error1Mad=c(), error2Mad=c(), error3Mad=c(), error4Mad=c(), error5Mad=c(), error6Mad=c())
     for(m in mags)
     {
          for(c in cells)
          {
               duh <- testLAPAccuracy_outOfFrame3(n=c, mags <- seq(m,m,length.out=reps))
               means <- colMedians(as.matrix(duh))
               names(means) <- c('poofMean', 'error1Mean', 'error2Mean', 'error3Mean', 'error4Mean', 'error5Mean', 'error6Mean')
               mads <- colMads(as.matrix(duh))
               names(mads) <- c('poofMad', 'error1Mad', 'error2Mad', 'error3Mad', 'error4Mad', 'error5Mad', 'error6Mad')
               results <- rbind(results, data.frame(mag=m, cells=c, as.list(means), as.list(mads)))
          }
     }
     return(results)
}

plotMethod <- function(results, name, normalize=FALSE)
{
     xlim <- c(0, max(results$poofMean))
     ylim <- c(0, max(results[,name]))
     if(normalize)
     {
          ylim <- c(0, max(results[,c('error1Mean','error2Mean','error3Mean','error4Mean','error5Mean','error6Mean')]))
     }

     mags <- unique(results$mag)
     cells <- unique(results$cells)
     plot(x=c(), y=c(), xlim=xlim, ylim=ylim, type='l', xlab='Points Out Of Frame [%]', ylab='Error Rate [%]')
     count <- 1
     for(c in cells)
     {
          data <- subset(results, cells==c)
          lines(data$poofMean, data[,name], lty=count)
          count = count + 1
     }
}

##### Track Cells #####

