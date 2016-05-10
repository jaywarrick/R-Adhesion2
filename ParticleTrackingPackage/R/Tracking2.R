library(Matrix)
library(matrixStats)
library(pracma)

#' @import Matrix
#' @import matrixStats
#' @importFrom pracma cross isempty numel
NULL

#' getStandardCostMatrix
#'
#' This function calls 'getDirectionalCostMatrix' with a default direction vector of along the first dimension
#' of the data and a perpendicularityPenaltyFactor of 1 (i.e., perpendicular motion is not penalazed). The absolute
#' value of the matrix is returned to remove the sign on each of the values returned by the function.
#'
#' @param data0 matrix where rows represent individual 'before points' to be linked (or not) with points from data1
#' @param data1 matrix where rows represent individual 'after points' to be linked (or not) with points from data0
#' @param maxDist The maximum allowable euclidean distance allowed for linking
#' @param digits numeric value indicating the number of places after the decimal point to maintain for precision
#'
#' @export
getStandardCostMatrix <- function(data0, data1, digits, maxDist)
{
     return(abs(getDirectionalCostMatrix( data0 = data0, data1 = data1, digits, maxDist = maxDist, direction = c(1,numeric(ncol(data0)-1)), perpendicularPenaltyFactor = 1)))
}

#' getDirectionalCostMatrix
#'
#' This function produces a cost matrix of the format described for point linking in the following paper.
#'
#' #' Khuloud Jaqaman, Dinah Loerke, Marcel Mettlen, Hirotaka Kuwata, Sergio Grinstein, Sandra L Schmid & Gaudenz Danuser
#' "Robust single-particle tracking in live-cell time-lapse sequences" Nature Methods - 5, 695 - 702 (2008).
#'
#' The max linkable cost is the square of the max distance. The cost used to block links is 10 times the
#' max linkable cost. The cost of not linking is set equal to the max linkable cost.
#'
#' This function retains the sign on the cost result so that it can be used in a 'directional' (i.e., is motion
#' in the same direction as the direction vector, positive, or in the opposite direction, negative)
#'
#' @param data0 matrix where rows represent individual 'before' points to be linked (or not) with points from data1
#' @param data1 matrix where rows represent individual 'after' points to be linked (or not) with points from data0
#' @param digits numeric value indicating how many places after the decimal to retain in terms of accuracy from the point information (applied to all dimensions of the data, scale data appropriately)
#' @param maxDist The maximum allowable euclidean distance allowed for linking
#' @param direction numeric vector with the same number of elements as the columns in x (i.e., spatial dimensions x, y, ...)
#' @param perpendicularPenaltyFactor numeric value to multiply the delta distance between two points that is perpendicular to the direction vector
#'
#' @export
getDirectionalCostMatrix <- function(data0, data1, digits, maxDist, direction, perpendicularPenaltyFactor)
{
     if(!is.matrix(data0) | !is.matrix(data1))
     {
          stop("data0 and data1 need to be matrices.")
     }
     if(ncol(data0) != ncol(data1))
     {
          stop("data0 and data1 must have the same number of columns")
     }
     if(ncol(data0) != base::length(direction))
     {
          stop("The direction vector must define the same number of dimensions as the data. Number of data colums.")
     }

     # Establish important costs
     maxLinkingCost <- maxDist*maxDist
     blockingCost <- 10*maxLinkingCost

     # Reformat the data
     # get the pair-wise combinations of the data by replicating the rows in each data
     # object appropriately
     dataReps <- getDataReps(data0, data1)

     # allow direct access to dataReps$x0 and dataReps$x1
     x0 <- dataReps$x0
     x1 <- dataReps$x1

     # Calculate the costs
     # make sure direction is a unit vector by dividing all components by its magnitude
     direction <- direction/(sqrt(sum(direction^2)))

     # Calculate movement vectors for each combination
     movement <- x1-x0
     parallelMovement <- movement%*%direction # dot product of direction with each row of the movement matrix
     movementSign <- sign(parallelMovement)
     perpendicularMovement2 <- rowSums(movement*movement) - parallelMovement*parallelMovement

     # Calculate cost
     cost <- parallelMovement*parallelMovement + ((perpendicularPenaltyFactor*perpendicularPenaltyFactor)*perpendicularMovement2) # Penalize perpendicular movement in proportion to the expected ratio of motion expected parallel vs perpendicular to the specified direction
     cost[cost >= maxLinkingCost] <- blockingCost

     # Just want to keep track of positive and negative. Sign can be 0 for 0's so just make them 1's
     # so that they don't result in a zero cost when multiplying the costs by the signs again.
     movementSign[movementSign == 0] <- 1
     cost <- cost*movementSign
     cost <- matrix(cost, nrow(data0), nrow(data1))

     # Assemble the cost matrix (apply the digits precision to the )
     results <- getCostMatrix(UL=cost, digits=digits, maxLinkingCost=maxLinkingCost, blockingCost=blockingCost)

     return(results)
}

#' standardLinearAssigment
#'
#' This function calls 'directionalLinearAssignment with a uniformityDistThresh value of -1 (i.e., no
#' uniformity constraint is applied.)
#'
#' @param pointSet0 PointSet object at the 'before' time for which to get link assignments using the provided link cost matrix
#' @param pointSet1 PointSet object at the 'after' time for which to get link assignments using the provided link cost matrix
#' @param costMatrix numeric matrix representing the costs associated with linking the points between pointSet0 to pointSet1 (See documentation for getDirectionalCostMatrix)
#'
#' @export
standardLinearAssignment <- function(pointSet0, pointSet1, costMatrix)
{
     px <- LinearAssignment(costMatrix)

     assignments <- getAssignments(pointSet0=pointSet0, pointSet1=pointSet1, px=px)

     return(assignments)
}

#' directionalLinearAssignment
#'
#' Link points in pointSet0 to pointSet1 with digits precision after the decimal place. This function can favor links made parallel
#' to a particular direction, specified by the direction vector. The direction vector is turned into a unit vector by normalizing
#' by its magnitude. Any "movement" perpendicular to the unit vector is multiplied by the "perpendicularPenaltyFactor". Thus, the total
#' cost of a link is cost = A^2 + (perpendicularPenaltyFactor*B)^2 where A is the euclidean distance parallel to the direction vector
#' and B is the euclidean distance perpendicular to the direction vector. This function can also apply a "uniformity" constraint.
#' That is to say, if all the points should be heading in either the positive or the negative direction along the direction vector,
#' then retrograde motion is deemed illegal. This uniformity constraint is applied when the average (penalized) euclidean "distance"
#' traveled by the points is greater than the "uniformityDistThresh".
#'
#' To avoid penalizing motion perpendicular to the direction vector, use perpendicularPenaltyFactor=1 (default value is 10). To avoid
#' applying a uniformity constraint, use uniformityDistThresh < 0 (e.g., -1, which is the default value).
#'
#' @param pointSet0 PointSet object. Links will be made from this PointSet to pointSet1
#' @param pointSet1 PointSet object. Links will be made from pointSet0 to this PointSet
#' @param digits numeric value indicating the number of places after the decimal place to be kept in regards to precision of point data and associated costs
#' @param maxDist numeric value indicating the maximum (penalized) euclidean distance of links that are deemed possible
#' @param direction numeric vector of the same length as the number of data columns representing the point "locations" in the PointSets (e.g., for position dimensions x and y -> direction=c(1,0)). The default is a direction vector with a value of 1 for the first data column and 0's for the remaining.
#' @param perpendicularPenaltyFactor numeric value used to multiply distance traveled perpendicular to the direction vector.
#' @param uniformityDistThresh numeric value represneting the threshold of the average euclidean "distance" traveled by each point. When the average distance is above this value, the uniformity threshold is applied, when below it is not (to allow for situations where there is no bulk motion)
#'
#' @return data.frame(id0=numeric, id1=numeric) where each row indicates which id in pointSet0 (id0) is linked to which id in pointSet1 (id1)
#'
#' @export
directionalLinearAssignment <- function(pointSet0, pointSet1, digits, maxDist, direction, perpendicularPenaltyFactor=10, uniformityDistThresh=-1)
{
     # print(paste('digits', digits))
     # print(paste('maxDist', maxDist))
     # cat('direction', direction))
     # print(paste('perpendicularPenaltyFactor', perpendicularPenaltyFactor))
     # print(paste('uniformityDistThresh', uniformityDistThresh))

     data0 <- pointSet0$getPointsWithoutId()
     data1 <- pointSet1$getPointsWithoutId()
     # Do first round of tracking as first guess
     # print(maxDist)
     results <- getDirectionalCostMatrix(data0=as.matrix(data0), data1=as.matrix(data1), digits=digits, maxDist=maxDist, direction=direction, perpendicularPenaltyFactor=perpendicularPenaltyFactor)

     px <- LinearAssignment(abs(results$costMatrix))

     # A negative number indicates that we should not apply a uniformity constraint
     if(uniformityDistThresh >= 0)
     {
          uniformityCostThresh <- uniformityDistThresh^2
          uniformityCostThresh <- round((10^digits)*uniformityCostThresh) # apply the 'digits' precision

          # First find valid rows and cols
          nrows <- nrow(data0)
          ncols <- nrow(data1)

          # Only look at pxs that reference links between t0 and t1 and not broken links
          validPx <- which(px <= nrows)

          # only care about linked points (this selects valid columns)
          validPx <- validPx[validPx <= ncols]

          # get the costs for each of the valid links
          linkedCosts <- numeric(length(validPx))
          for(i in 1:length(validPx))
          {
               linkedCosts[i] <- results$costMatrix[px[validPx[i]],validPx[i]]
          }

          # determine which sign to penalize/remove (i.e., replace with blocking cost)
          linkedSign <- -1*sign(sum(sign(linkedCosts)))

          # catches condition where there are essentially all blocking costs (i.e., no valid links, so don't apply constraint)
          if(!is.na(linkedSign))
          {
               # block connections in the wrong direction
               if(linkedSign > 0)
               {
                    results$costMatrix[1:nrows, 1:ncols][results$costMatrix[1:nrows, 1:ncols] > uniformityCostThresh] <- results$blockingCost
               }
               else
               {
                    results$costMatrix[1:nrows, 1:ncols][results$costMatrix[1:nrows, 1:ncols] < -1*uniformityCostThresh] <- results$blockingCost
               }

               # redo the linking procedure with the additional blocked links from the uniformity constraint
               px <- LinearAssignment(abs(results$costMatrix))
          }
          else
          {
               # i.e., there are no valid links so there is no point in applying the uniformity constraint
               print('No valid links found, not applying uniformity constraint')
          }
     }
     ret <- getAssignments(pointSet0=pointSet0, pointSet1=pointSet1, px=px)
     return(ret)
}

#' getAssignments
#'
#' Get the ids of the points associated with a linear assigment result. A list of objects is returned.
#'
#' links = data.frame(id0, id1) indicating which ids in pointSet0 are linked with which ids in pointSet1
#'
#' noLinks = numeric vector of ids in pointSet0 which are not linked to any point in pointSet1
#'
#' newTrack = numeric vector of ids in pointSet1 which are not linked to any point in pointSet0, thus
#' representing the beginning of a new track
#'
#' @param pointSet0 PointSet object from the 'before' time (startpoint of link)
#' @param pointSet1 PointSet object from the 'after' time (endpoint of link)
#' @param linearAssigmentResult data.frame with columns id0 and id1 indicating which ids from pointSet0 are linked to pointSet1. Others are not linked.
#'
#' @export
getAssignments <- function(pointSet0, pointSet1, px)
{
     links <- data.frame(id0=numeric(0), id1=numeric(0))
     noLinks <- numeric(0)
     newTracks <- numeric(0)
     n0 <- pointSet0$pointCount()
     n1 <- pointSet1$pointCount()
     for(i in 1:base::length(px))
     {
          # id0 is essentially px while id1 is essentially the index in px
          if(px[i] <= n0 & i <= n1)
          {
               # print("A")
               # Then in the upper-left hand quadrant
               # Then we have a link between x1 and x0 and the points should share ids
               id0 <- pointSet0$getPoint(index=px[i])$id
               id1 <- pointSet1$getPoint(index=i)$id
               links <- rbind(links, data.frame(id0=id0, id1=id1))
          }
          else if(px[i] > n0 & i <= n1)
          {
               # print("B")
               # Then in the lower-left hand quadrant
               # This point in x1 represents the start of a new track
               # LL gets moved up into first n0 rows means new track
               id1 <- pointSet1$getPoint(index=i)$id
               newTracks <- c(newTracks, id1)
          }
          else if(px[i] <= n0 & i > n1)
          {
               # print("C")
               # Then in the upper right-hand quadrant...
               # We failed to link and it isn't a new track so if it is a legitamate id0, then it is a no link
               # This point in x0 represents the end of a track or an auxillary match to enable other types of matches
               id0 <- pointSet0$getPoint(index=px[i])$id
               noLinks <- c(noLinks, id0)
          }
     }
     return(list(links=links, id0_NoLinks=noLinks, id1_newTracks=newTracks, px=px))
}

#' getDataReps
#'
#' A matrix of data0 data is created by replicating the data using data0Seq <- rep(1:nrow(data0), times=nData1).
#' A matrix of data1 data is created by replicating the data using data1Seq <- rep(1:nrow(data1), each=nData0).
#' Now pair-wise calculations of the two data-sets can be performed using the new matrices.
#'
#' @param data0 matrix of starting poitns with rows representing points and columns representing dimensions
#' @param data1 matrix of ending points with rows representing points and columns representing dimensions
#'
#' @return list(x0=newData0, x1=newData1)
#'
#' @export
getDataReps <- function(data0, data1)
{
     if(!is.matrix(data0) | !is.matrix(data1))
     {
          stop("data0 and data1 must be matricies")
     }
     nData0 <- nrow(data0)
     nData1 <- nrow(data1)
     data0Seq <- rep(1:nrow(data0), times=nData1) ##
     data0Rep <- data0[data0Seq,]
     data1Seq <- rep(1:nrow(data1), each=nData0)
     data1Rep <- data1[data1Seq,]
     if(!is.matrix(data0Rep))
     {
          data0Rep <- as.matrix(as.data.frame(as.list(data0Rep)))
     }
     if(!is.matrix(data1Rep))
     {
          data1Rep <- as.matrix(as.data.frame(as.list(data1Rep)))
     }
     return(list(x0=data0Rep, x1=data1Rep))
}

#' getCostMatrix
#'
#' Create a cost matrix with columns that represent id's of points from time t+1 followed id's of points from time t.
#' The columns represent id's of points from time t followed by id's of points from time t+1.
#'
#' @param UL matrix representing the link costs and Upper Left of the final cost matrix
#' @param digits numeric value representing the number of digits of precision after the decimal point (ok to be negative)
#' @param maxLinkingCost numeric value that is max possible linking cost (typically maxDist^2)
#' @param blockingCost numeric value that the cost to use to block links etc (typically 10*maxDist^2)
#'
#' @export
getCostMatrix <- function(UL, digits, maxLinkingCost, blockingCost)
{
     UR <- getURorLL(n=nrow(UL), maxLinkingCost=maxLinkingCost, blockingCost=blockingCost)
     LL <- getURorLL(n=ncol(UL), maxLinkingCost=maxLinkingCost, blockingCost=blockingCost)
     LR <- getLR(UL=UL, maxLinkingCost=maxLinkingCost, blockingCost=blockingCost)
     ret <- rbind(UL,LL)
     temp <- rbind(UR,LR)
     ret <- cbind(ret,temp)
     ret <- round((10^digits)*ret) # turn it into an integer problem for the LinearAssignment function of the GraphAlignment package
     if(max(abs(ret)) > .Machine$integer.max)
     {
          stop("Reduce the value for digits. This much precision is resulting in costs that are too high for the integer representation used by GraphAlignment::LinearAssigment. Negative values are acceptable for the digits parameter.")
     }
     return(list(costMatrix=ret, maxLinkingCost=round((10^digits)*maxLinkingCost), blockingCost=round((10^digits)*blockingCost)))
}

getURorLL <- function(n, maxLinkingCost, blockingCost)
{
     UR <- matrix(blockingCost, n, n);

     # Set the maxLinkingCost along the diagonal (top left to bottom right)
     if(n > 0)
     {
          for (i in 1:n)
          {
               UR[i, i] <- maxLinkingCost
          }
     }
     return(UR);
}

getLR <- function(UL, maxLinkingCost, blockingCost)
{
     LR = t(UL)

     r <- nrow(LR)
     c <- ncol(LR)

     LR[LR < blockingCost] <- maxLinkingCost

     return(LR)
}
