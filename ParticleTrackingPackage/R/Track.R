#' @import methods
NULL

#' Track
#'
#' A class for storing a table of linked pts. All Track objects should have column called
#' 'frame' and 't' in their 'pts' data.frame field to store the frame number and time
#' associated with each point.
#'
#' @title A class representing a Track or locations over time for a specific particle
#' @field id numeric ID of the particle being tracked
#' @field pts data.frame table of x, y, t, frame, vx, vy and optionally smoothed versions of vx and vy called vxs and vys (see method smoothVelocities)
#' @field meta list of Track information
Track <- setRefClass('Track',
                     fields = list(id='numeric', pts='data.frame', meta='list'),
                     methods = list(
                          initializeWithTrackROI = function(id, start, pattern, t0_Frame, timePerFrame)
                          {
                               "Using the value in the 'PolygonPts' value and the 'PatternPts'
                               value from a JEX ROI that represents a track, create an R 'Track'
                               object. Assign the track the id given as 'id'."

                               pairs <- strsplit(pattern,';')[[1]]
                               x <- numeric(0)
                               x0 <- numeric(0)
                               y <- numeric(0)
                               y0 <- numeric(0)
                               frames <- numeric(0)
                               first <- TRUE
                               for(pair in pairs)
                               {
                                    if(first)
                                    {
                                         nums <- strsplit(start,',')[[1]]
                                         x0 <- as.numeric(nums[1])
                                         y0 <- as.numeric(nums[2])
                                         i0 <- as.numeric(nums[3])
                                         first <- FALSE
                                         x <- append(x,x0)
                                         y <- append(y,y0)
                                         frames <- append(frames,i0)
                                    }
                                    else
                                    {
                                         nums <- strsplit(pair,',')[[1]]
                                         x <- append(x, x0+as.numeric(nums[1]))
                                         y <- append(y, y0+as.numeric(nums[2]))
                                         frames <- append(frames,as.numeric(nums[3]))
                                         # print(nums)
                                    }
                               }
                               id <<- as.numeric(as.character(id))
                               pts <<- data.frame(frame=frames, x=x, y=y)
                               updateTimes(t0_Frame, timePerFrame)
                          },
                          updateTimes = function(t0_Frame, timePerFrame)
                          {
                               pts$t <<- (pts$frame - t0_Frame) * timePerFrame
                          },
                          frameCount = function()
                          {
                               "Return the number of frames in this track."

                               return(nrow(pts))
                          },
                          getSlot = function(slot, rel=FALSE, selectedFrames=NULL)
                          {
                               "Get a vector of the values indicated by 'slot' for this point.\n
                               @param slot character value indicating the slot for which to obtain data within the 'pts' data.frame field within the Track\n
                               @param rel logical value indicating whether to get the slot information relative to the mean of that slot information (default = FALSE)\n
                               @param selectedFrames numeric vector of frame numbers for which to grab the slot information (default=NULL, which gets all the slot information available for this track)"

                               if(!is.null(selectedFrames))
                               {
                                    validIndices <- which(pts$frame %in% selectedFrames)
                               }
                               else
                               {
                                    if(base::length(pts$frame) > 0)
                                    {
                                         validIndices <- 1:base::length(pts$frame)
                                    }
                                    else
                                    {
                                         validIndices <- numeric(0)
                                    }
                               }

                               if(rel)
                               {
                                    temp <- pts[,slot] - mean(pts[,slot])
                                    temp <- temp[validIndices] # frames start at 0 while R indices start at 1
                               }
                               else
                               {
                                    temp <- pts[,slot]
                                    temp <- temp[validIndices]
                               }
                               temp <- temp[!is.na(temp)]
                               return(temp)
                          },
                          range = function(slot, rel=FALSE)
                          {
                               "For the indicated 'slot' (i.e., column within the pts field),
                               return the range that this value takes.\n
                               @param slot character value indicating the slot for which to obtain data within the 'pts' data.frame field within the Track\n
                               @param rel logical value indicating whether to get the slot information relative to the mean of that slot information (default = FALSE)"

                               ret <- base::range(getSlot(slot, rel=rel))
                               return(ret)
                          },
                          calculateDerivatives = function(slots, withRespectTo='t', prefix='v')
                          {
                               "Internally populate/update the time derivative calculations for each
                               raw data columns listed in the character vector 'slots' using the
                               'getDerivative(x, t)' function. The derivative is performed relative to
                               another data column (e.g., 't' by default) specified in the 'withRespectTo'
                               parameter. The calcualted derivatives are stored in the 'pts' data.frame
                               field. These new columns are the names of the original data columns specified
                               by the slots parameter prepended by the user-specified 'prefix' parameter
                               (default, prefix = 'v')."

                               slots <- slots[slots %in% names(pts)]
                               if(base::length(slots) > 0)
                               {
                                    for(i in 1:base::length(slots))
                                    {
                                         pts[,paste0(prefix, slots[i])] <<- getDerivative(pts[,slots[i]], pts[,withRespectTo])
                                    }
                               }
                               else
                               {
                                    stop("Track::calculateDerivatives --> Couldn't find a slot matching any of the slots provided. Aborting.")
                               }
                          },
                          calculateSmoothedData = function(windowWidths, slots, suffix='s')
                          {
                               "Smooth the specified slots for this track\n
                               @param windowWidths A data.frame(frame=numeric, width=numeric) that at least contains a width corresponding to each frame in this Track, or a single value indicating the width to be used for all frames.\n
                               @param slots character vector of all the columns in the 'pts' data.frame field for which smoothed versions of the data should be calculated."

                               if(is.data.frame(windowWidths))
                               {
                                    if(sum(names(windowWidths) %in% c('frame','width')) == 2)
                                    {
                                         # Restrict the analysis to just the frames in this track
                                         widths <- windowWidths[windowWidths$frame %in% pts$frame,'width']
                                         if(base::length(widths) != base::length(pts$frame))
                                         {
                                              stop("Track::calculateSmoothedData --> the data.frame provided does not specify a width for each frame in this Track. Aborting.")
                                         }
                                    }
                                    else
                                    {
                                         stop("Track::calculateSmoothedData --> The supplied windowWidths data.frame does not contain a column labeled 'frame' and 'width'. Aborting 'calculateSmoothedData'. Aborting.")
                                    }
                               }
                               else if(is.vector(windowWidths) && is.numeric(windowWidths) && base::length(windowWidths) == 1)
                               {
                                    widths <- rep(widths, nrow(pts))
                               }
                               else
                               {
                                    stop("Track::calculateSmoothedData --> windowWidths must be a data.frame or single numeric value. Aborting.")
                               }

                               # Now working with the 'widths' variable instead of the 'windowWidths' variable
                               if(!all(is.finite(widths)))
                               {
                                    stop("Track::calculateSmoothedData --> some window widths were not finite. Aborting.")
                               }
                               if(min(widths) < 2)
                               {
                                    warning("Track::calculateSmoothedData --> There exists a width in the specified widths that is smaller that the minimum width of 2 (i.e., two frames to average). Setting those widths to 2 to proceed.")
                                    widths[widths < 2] <- 2
                               }
                               if(base::length(widths) != nrow(pts))
                               {
                                    stop("The widths provided are not specified for all the points in the track. Aborting.")
                               }

                               # Actually get down to business
                               slots <- slots[slots %in% names(pts)] # Find as many matching slots as possible.
                               for(slot in slots)
                               {
                                    newSlot <- paste0(slot, suffix)
                                    pts[,newSlot] <<- sapply(seq_along(pts[,slot]),
                                                              getAverage,
                                                              widths=widths,
                                                              x=pts[,slot]
                                    )
                               }
                          },
                          plotTrack = function(slotX='frame', slotY='x', funX=NULL, funY=NULL, relX=FALSE, relY=TRUE, selectedFrames=NULL, add=FALSE, withTitle=TRUE, col='black', lwd=1, lty=1, xlab=slotX, ylab=slotY, type='l', ...)
                          {
                               "Plot the track\n
                               @param slotX The x variable string name\n
                               @param slotY The y variable string name\n
                               @param relX A boolean indicating whether to plot x relative to its mean\n
                               @param relY A boolean indicating whether to plot y relative to its mean\n
                               @param funX function A function to transform the slotX data before plotting (e.g., to scale between pixels per second velocity to microns per second velocity)\n
                               @param funY function A function transform the slotY data before plotting (see funX)\n
                               @param ... Additional args are passed to the 'plot' method\n
                               @param selectedFrames numeric vector indicating which frames to plot or NULL (default) to plot all frames in this track"

                               xData <- getSlot(slotX, relX, selectedFrames=selectedFrames)
                               if(!is.null(funX))
                               {
                                    xData <- funX(xData)
                               }
                               yData <- getSlot(slotY, relY, selectedFrames=selectedFrames)
                               if(!is.null(funY))
                               {
                                    yData <- funY(yData)
                               }
                               if(base::length(xData) < 0 || base::length(yData) < 0)
                               {
                                    message(paste('No valid data to plot for track ', id, '.', sep=''))
                               }
                               if(relX)
                               {
                                    xlab <- paste(xlab, ' - mean(', xlab, ')', sep='')
                               }
                               if(relY)
                               {
                                    ylab <- paste(ylab, ' - mean(', ylab, ')', sep='')
                               }

                               if(isempty(xData) || isempty(yData) || is.na(xData) || is.na(yData))
                               {
                                    message(paste('Nothing to plot for track ', id, '.', sep=''))
                               } else
                               {
                                    if(add)
                                    {
                                         graphics::lines(xData, yData, col=col, lwd=lwd, lty=lty, type=type, ...)
                                    }
                                    else
                                    {
                                         if(withTitle)
                                         {
                                              graphics::plot(xData, yData, col=col, lwd=lwd, lty=lty, xlab=xlab, ylab=ylab, main=as.character(id), type=type, ...)
                                         }
                                         else
                                         {
                                              graphics::plot(xData, yData, col=col, lwd=lwd, lty=lty, xlab=xlab, ylab=ylab, type=type, ...)
                                         }
                                    }
                               }
                          },
                          addPoint = function(frame, ...)
                          {
                               "Add a point to the 'pts' field of this track. It is
                               assumed you are adding things in an appropriate order.
                               Typically the pts are listed in order of frame.\n
                               @param frame A single numeric value indicating the frame corresponding to this frame/point\n
                               @param ... A variable list of numeric values that represent the point's data (e.g., x and y position)\n"

                               pts <<- rbind(pts, data.frame(frame=frame, ...))
                          },
                          setMeta = function(meta)
                          {
                               "Basically a convenience function for setting the meta from a parent TrackList object but fine to use otherwise as well."

                               meta <<- meta
                          }
                     )
)

#' getAverage
#'
#' Get the adjustable running window average of the data
#'
#' @param i The index within 'frames' at which to calculate an average over a window centered at this location
#' @param frames The frames in this track
#' @param widths A vector of window widths appropriate for each frame in the 'frames' of this track
#' @param data The vector of data for which we will calculate the windowed averages
getAverage <- function(i, widths, x)
{
     # Subtract 1 to represent the number of intervals instead of number of pts to average
     width <- widths[i] - 1

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
     mean(x[leftIndex:(leftIndex+width)])
}

#' getDerivative
#'
#' Get the derivative of a vector
#'
#' @param x A numeric vector on which to calculate the derivative
#' @param t A numeric vecotor of times with which to determine dt for derivative calculations
#'
#' @export
getDerivative <- function(x, t)
{
     v <- numeric(0)
     for(i in 1:length(x))
     {
          v <- c(v, localDerivative(x, t, i))
     }
     return(v)
}

#' localDerivative
#'
#' Get the local derivative around a point in a vector accounding for boundary scenarios at the start and end of the vector
#'
#' @param x A numeric vector of data
#' @param t A numeric vector of time for calculating dt of the derivative
#' @param i A numeric value indicating the index in the x and t for which to calculate the local derivative
#'
#' @export
localDerivative <- function(x, t, i)
{
     if(i == 1)
     {
          #return((x[i+1]-x[i])/(t[i+1]-t[i]))
          return(interpolateDerivative(x[i], x[i+1], x[i+2], t[i], t[i+1], t[i+2], t[i]))
     }
     else if(i == length(x))
     {
          #return((x[i]-x[i-1])/(t[i]-t[i-1]))
          return(interpolateDerivative(x[i-2], x[i-1], x[i], t[i-2], t[i-1], t[i], t[i]))
     }
     else
     {
          return(interpolateDerivative(x[i-1], x[i], x[i+1], t[i-1], t[i], t[i+1], t[i]))
     }
}

#' interpolateDerivative
#'
#' This is a three point interpolation of the derivative where the interpolated point
#' is the middle of the 3 pts. This simplifies to the the three-point midpoint formula
#' when the time steps are equal but can handle when timesteps are unequal (i.e., the
#' time-step on either side of the 3 pts is not equal)
#'
#' @param f0 numeric left function value
#' @param f1 numeric middle function value
#' @param f2 numeric right function value
#' @param x0 numeric left independent value
#' @param x1 numeric middle independent value
#' @param x2 numeric right independent value
#' @param xj numeric x value for which to evaluate the function
#'
#' @export
interpolateDerivative <- function(f0, f1, f2, x0, x1, x2, xj)
{
     term1 <- f0*((2*xj-x1-x2)/((x0-x1)*(x0-x2)))
     term2 <- f1*((2*xj-x0-x2)/((x1-x0)*(x1-x2)))
     term3 <- f2*((2*xj-x0-x1)/((x2-x0)*(x2-x1)))
     return(term1 + term2 + term3)
}
