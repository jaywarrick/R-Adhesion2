#' getWindowWidths
#'
#' Calculate the appropriate window widths for smoothing velocity calculations.
#'
#' Underlying Calculations: v represents the velocity calculated over a single period. Thus,
#' the cell should travel 4*A units of distance per period. This is used to estimate the velocity
#' of the cell using the fit. Therefore, the expected distance the cell will travel between frames
#' is v*dt where dt is the time between the frames of interest. Thus, if we would like to
#' quantify velocity over a travel distance of 10 pixels (i.e., dist=10), then we can take
#' dist/v*dt to get the number of frames (i.e., dt's) needed to cover 10 pixels. However, if maxWidth
#' is exceeded, maxWidth is returned. dt could be variable over time depending on frame sampling rates etc.
#' approximate a dt for each index based on diff's. Repeat last value at end to create a dt
#' vector that is the same as the trackList$meta$allTimes vector
#'
#' @param trackList TrackList object for which to calculate appropriate window widths
#' @param fit The fit results returned by 'getBulkPhaseShift'
#' @param dist The number of pixels that we would like to see (if possible) the particle/cell move before estimating velocity
#' @param maxWidth numeric The maximum number of FRAMES that should be averaged together to estimate velocity (i.e., averaging is good up to a point where it just increases compuation time)
#'
#' @export
getWindowWidths <- function(trackList, fit, dist, maxWidth)
{
     dt <- diff(trackList$meta$allTimes)
     dt <- c(dt, last(dt))
     v <- 4*fit$par[['amplitude']]*(trackList$meta$fi*(trackList$meta$ff/trackList$meta$fi)^(trackList$meta$allTimes/last(trackList$meta$allTimes)))
     widths <- ceiling(dist/(v*dt))
     widths[widths > maxWidth] <- maxWidth
     return(data.frame(frame=trackList$meta$allFrames, width=widths))
}

#' calculateValidFrames
#'
#' This function is intended for use with particles that have periodic motion as defined by the
#' 'getSweep' function. The function takes a fit object for the velocity and determines where the cells
#' switch directions 'validStart' and 'validEnd' represent precentages. In other words, after the switch
#' in direction, valid points start at 'validStart' % of the way to the next switch in direction while
#' 'validEnd' occurs 'validEnd' % of the way to the that same next switch in direction.
#'
#' This filter is good for removing inaccurate values of the velocity when using a triangle waveform because
#' the estimate of velocity will be artificially be lower due to sample aliasing of particle motion.
#' Possible values are 0, 1, 2, and 3 for inflectionPtsToFilter. Multiple at once can be provided.
#'
#' 0 Represents the upward zero-crossing of the position and max positive velocity
#' 1 Represents the max position and downward zero-crossing of the velocity
#' 2 Represents the downward zero-crossing of the position and min (i.e most negative) velocity
#' 4 Represents the min (i.e., most negative) position and upward zero-crossing of the velocity
#'
#' 'inflectionPtsToFilter' refers to the first point of the numbered sections (i.e., 1 refers to the max point)
#' Thus if 1 is included in 'inflectionPtsToFilter', points adjacent to this point will be considered for wobble and filtered.
#' We use 'inflectionPtsToFilter' of c(1,3) to mark the points at which the flow switches direction
#'
#' @param trackList TrackList object for which to calculate valid frames
#' @param fit fit object that represents the fit of a dsecending log frequency sweep to the data.
#' @param validStart numeric value (see above)
#' @param validEnd numeric value (see above)
#'
#' @export
calculateValidFrames <- function(trackList, fit, validStart=0.01, validEnd=0.99)
{
     # Helper function: get data points adjacent to specified inflection timepoint
     getNearests <- function(t, inflectionPoint)
     {
          upper <- which.max((t-inflectionPoint) >= 0)
          if(is.na(upper) || upper == 1)
          {
               return(NA)
          }
          lower <- upper - 1
          return(c(lower, upper))
     }

     sweep <- getSweep(amplitude=fit$par['amplitude'], phaseShift=fit$par['phaseShift'], offset=0, sin=trackList$meta$sin, ti=fit$par['ti'], fi=trackList$meta$fi, ff=trackList$meta$ff, sweepDuration=trackList$meta$sweepDuration, t=trackList$meta$allTimes, guess=NULL, flipped=TRUE, calcVelocity=TRUE)
     inflectionsToAddress <- sweep$inflectionNums %in% c(1,3) # These are times at which flow switches directions
     indicesToRemove <- numeric(0)
     for(i in which(inflectionsToAddress))
     {
          nearests <- getNearests(sweep$t, sweep$inflections[i])
          if(!is.na(nearests)[1])
          {
               indicesToRemove <- c(indicesToRemove, nearests)
          }
     }
     selectedFrames0 <- trackList$meta$allFrames
     selectedFrames0 <- selectedFrames0[-indicesToRemove]

     selectedFrames1 <- numeric(0)
     for(tIndex in 1:base::length(trackList$meta$allTimes))
     {
          # Get the nearest inflection at or beyond this time

          infIndex <- which.max((sweep$inflections >= trackList$meta$allTimes[tIndex]) & inflectionsToAddress)
          if(is.na(infIndex)) next

          # Get the bounding inflections that represent changes in fluid direction
          infT2 <- sweep$inflections[infIndex] # take the inflection we found
          if((infIndex-2) < 1)
          {
               infT1 <- 0
          }
          else
          {
               infT1 <- sweep$inflections[infIndex-2] # also take two inflections prior because each inflection represents pi/2 and we want to go back to the last change in direction which is pi ago.
          }
          dInfT <- infT2-infT1 # define the time interval size between these two inflections

          # Within the if statement calculate the fractional location of this time index in the interval between the two inflections.
          if( (trackList$meta$allTimes[tIndex] >= (infT1 + validStart*dInfT)) & (trackList$meta$allTimes[tIndex] <= (infT1 + validEnd*dInfT)) )
          {
               # If it is within the startValid and endValid bounds, add it to the list of the valid frames
               selectedFrames1 <- c(selectedFrames1, trackList$meta$allFrames[tIndex]) # (tIndex-1) = frame because frames are indices that start at 0
          }
     }

     return(sort(intersect(selectedFrames1, selectedFrames0)))
}

#' smoothData
#'
#' Smooth the data in a numeric vector and return a vector of the same length. The sliding
#' window is simply truncated at both ends as it runs into the bounds.
#'
#' @param windowWidth numeric value indicating the width of the window to be used for averaging at each index of the supplied data vector
#' @param x numeric vector of data to smooth
#'
#' @export
smoothData <- function(x, windowWidth=2)
{
     "Smooth the specified data with a running window average of the specified width\n
     @param windowWidths A data.frame(frame=numeric, width=numeric) that at least contains a width corresponding to each frame in this Track, or a single value indicating the width to be used for all frames.\n"

     if(windowWidth < 2)
     {
          stop("The windowWidth is smaller that the minimum width of 2. Aborting.")
     }

     ret <- sapply(seq_along(x),
                              smoothPt,
                              windowWidth=windowWidth,
                              x=x
     )
     return(ret)
}

#' smoothPt
#'
#' Get the adjustable running window average of the data
#'
#' @param i The index within 'frames' at which to calculate an average over a window centered at this location
#' @param windowWidth A vector of window widths appropriate for each frame in the 'frames' of this track
#' @param x The vector of data for which we will calculate the windowed averages
#'
#' @export
smoothPt <- function(i, windowWidth, x)
{
     # Subtract 1 to represent the number of intervals instead of number of pts to average
     width <- windowWidth - 1

     # calculate the index on the left of the interval
     leftIndex <- i - floor(width/2)
     if(leftIndex < 1)
     {
          leftIndex <- 1
     }

     # calculate the width to reach index on the right of the interval
     if((leftIndex+width) > length(x))
     {
          width <- length(x)-leftIndex
     }

     # return the mean of the x over the interval
     return(mean(x[leftIndex:(leftIndex+width)]))
}

#' getPercentAdhered
#'
#' Calculate the percent adhered at each timepoint using the provided velocity threshold
#'
#' @param trackList TrackList object to be analyzed
#' @param slot character value deignating which column of data in the 'points' field of each track to use to determine percent adhered (default 'vx', non-smoothed velocity)
#' @param velocityThreshold numeric The pixels per second below which (exclusive) a cell is considered adhered - default=3 [pixels/second]
#' @param windowWidth numeric value indicating the number of frames of cell motion to perform a sliding window average on the cell count data to account for changes in cell number over time in the calculation of percent adhered (default = 30).
#'
#' @export
getPercentAdhered <- function(trackList, slot='vxs', velocityThreshold=3, windowWidth=30)
{
     trackMatrix <- trackList$getMatrix(slot=slot, selectedOnly=TRUE, rel=FALSE)
     ret <- list()

     if(is.null(trackList$getSelectedFrames()))
     {
          stop("The selected frames of this TrackList must be set first. See 'setSelectedFrames' of the TrackList class.")
     }
     trackList$meta$pointCounts$smoothedVal <- smoothData(trackList$meta$pointCounts$val, windowWidth)

     frames <- colnames(trackMatrix)
     for(frame in frames)
     {
          velocities <- trackMatrix[,frame]
          velocities <- abs(velocities[!is.na(velocities)]) # Need the !is.na clause because we are looking across tracks. Within a track, we are garanteed to have a velocity for every frame.
          if(!isempty(velocities))
          {
               adhered <- sum(velocities < velocityThreshold)/trackList$meta$pointCounts$smoothedVal[trackList$meta$pointCounts$frame == frame]
               ret[[frame]] <- adhered
          }
          else
          {
               ret[[frame]] <- 0
          }
     }
     percents <- 100*as.numeric(ret)
     times <- trackList$meta$allTimes[trackList$meta$allFrames %in% trackList$getSelectedFrames()]
     return(data.frame(time=times, percentAdhered=percents))
}

