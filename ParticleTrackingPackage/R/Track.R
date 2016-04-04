#' @title A class representing a Track or locations over time for a specific particle
#' @field id numeric ID of the particle being tracked
#' @field points data.frame table of x, y, t, frame, vx, vy and optionally smoothed versions of vx and vy called vxs and vys (see method smoothVelocities)
#' @field validFrames numeric vector Vector frame numbers of the track that are set to 'valid' using 'setValidFrames' method
#' @field meta TrackList information from the parent TrackList object of this Track object
Track <- setRefClass('Track',
                     fields = list(id='numeric', points='data.frame', validFrames='numeric', meta='list'),
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
                               #                                if(base::length(frames) > base::length(tAll))
                               #                                {
                               #                                     stop(cat("Length of tAll seems to short for the data that is being used to initialize the track. frames length = ", base::length(frames), " tAll length = ", base::length(tAll)))
                               #                                }
                               points <<- data.frame(frame=frames, x=x, y=y, t=(frames - t0_Frame) * timePerFrame)
                               calculateVelocities()
                          },
                          setMeta = function(meta)
                          {
                               meta <<- meta
                          },
                          length = function()
                          {
                               "Return the number of frames in this track."

                               return(nrow(points))
                          },
                          getSlot = function(slot, rel=FALSE, validOnly=FALSE)
                          {
                               "Get a vector of the values indicated by 'slot' for this point.
                               This accessor method is useful because the data is stored within
                               a data.frame field called points.
                               'validOnly' indicates whether to return values for the validFrames only or all frames
                               typical slots are x, y, t, frame, vx, vy"

                               if(isempty(validFrames) & validOnly)
                               {
                                    cat("Track ID: ", id, " - Can't determine valid frames as this track either doesn't have any frames that are actually valid or they have not been set yet for this track\n")
                                    return(numeric(0))
                               }
                               if(validOnly)
                               {
                                    validIndices <- points$frame %in% validFrames
                               }
                               else
                               {
                                    validIndices <- points$frame
                               }


                               if(rel)
                               {
                                    temp <- points[,slot] - mean(points[,slot])
                                    if(validOnly)
                                    {
                                         temp <- temp[validIndices] # frames start at 0 while R indices start at 1
                                    }
                               }
                               else
                               {
                                    temp <- points[,slot]
                                    if(validOnly)
                                    {
                                         temp <- temp[validIndices]
                                    }
                               }
                               temp <- temp[!is.na(temp)]
                               return(temp)
                          },
                          range = function(slot, rel=FALSE)
                          {
                               "For the indicated 'slot' (i.e., column within the points field),
                               return the range that this value takes."

                               ret <- base::range(getSlot(slot, rel=rel))
                               return(ret)
                          },
                          calculateVelocities = function()
                          {
                               "Internally populate/update the values 'vx' and 'vy' within
                               the points vield using 'getDerivative(x, y)'"

                               points$vx <<- getDerivative(points$x, points$t)
                               points$vy <<- getDerivative(points$y, points$t)
                          },
                          smoothVelocities = function(allWidths, allFrames)
                          {
                               "Smooth the vx and vy velocity calculations for this track\n
                               @param allWidths A numeric vector of all the appropriate widths for each possible frame in the parent TrackList\n
                               @param allFrames A numeric vector of all the possible frames of the parent TrackList"

                               # Restrict the analysis to just the frames in this track
                               widths <- allWidths[points$frame %in% allFrames]

                               # Apply it to the vx data
                               points$vxs <<- sapply(seq_along(points$frame),
                                                     getAverage,
                                                     frames=points$frame,
                                                     widths=widths,
                                                     data=points$vx
                               )

                               # Apply it to the vy data
                               points$vys <<- sapply(seq_along(points$frame),
                                                     getAverage,
                                                     frames=points$frame,
                                                     widths=widths,
                                                     data=points$vy
                               )
                          },
                          plotTrack = function(slotX='t', slotY='x', funX=NULL, funY=NULL, relX=FALSE, relY=TRUE, validOnly=FALSE, add=FALSE, withTitle=TRUE, col='black', lwd=1, lty=1, xlab=slotX, ylab=slotY, type='l', ...)
                          {
                               "Plot the track\n
                               @param slotX The x variable string name\n
                               @param slotY The y variable string name\n
                               @param relX A boolean indicating whether to plot x relative to its mean\n
                               @param relY A boolean indicating whether to plot y relative to its mean\n
                               @param funX function A function to transform the slotX data before plotting (e.g., to scale between pixels per second velocity to microns per second velocity)\n
                               @param funY function A function transform the slotY data before plotting (see funX)\n
                               @param ... Additional args are passed to the 'plot' method\n
                               @param validOnly A boolean indicating whether to plot 'valid frames' only (see setValidFrames) or all frames in this track"

                               xData <- getSlot(slotX, relX, validOnly=validOnly)
                               if(!is.null(funX))
                               {
                                    xData <- funX(xData)
                               }
                               yData <- getSlot(slotY, relY, validOnly=validOnly)
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
                               if(base::length(list(...)$validOnly)>0)
                               {
                                    print(list(...))
                                    stop("what?")
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
                          addPoint = function(x, y, t, frame)
                          {
                               "Add a point to the 'points' field of this track. It is
                               assumed you are adding things in an appropriate order.
                               Typically the points are listed in order of frame.\n
                               @param x A single numeric value indicating the x position of the point\n
                               @param y A single numeric value indicating the y position of the point\n
                               @param t A single numeric value indicating the time corresponding to this frame/point\n
                               @param frame A single numeric value indicating the frame corresponding to this frame/point"

                               points <<- rbind(points, data.frame(x=x, y=y, t=t, frame=frame))
                          },
                          #                           show = function()
                          #                           {
                          #                                "
                          #                                #' Prints the information of this object to the command line
                          #                                #' output. This overrides the basic 'show' method as a track
                          #                                #' carries a reference to its parent trackList (if it has been
                          #                                #' added to a trackList using the 'setTrack' method) which
                          #                                #' results in recursive printing of TrackList information
                          #                                #' because a 'show' for 'TrackList' calls 'show' on each 'Track'.
                          #                                #' This method eliminates this issue.
                          #                                "
                          #                                cat("Reference class object of class 'Track'\n")
                          #                                for(name in names(Track$fields()))
                          #                                {
                          #                                     if(name %in% c('meta'))#c('x','y','t','vx','vy')))
                          #                                     {
                          #                                          cat("Field '", name, "':\n", sep='')
                          #                                          cat("(can't print... recursive)\n")
                          #                                     }
                          #                                     else
                          #                                     {
                          #                                          cat("Field '", name, "':\n", sep='')
                          #                                          methods::show(.self[[name]])
                          #                                     }
                          #                                }
                          #                                cat("\n")
                          #                           },
                          setValidFrames = function(frames)
                          {
                               "Set the 'valid frames' of the Track. Valid frames are those frames
                               which correspond to steady motion of the particle (i.e., not when
                               the particle/cell is switching directions due to changes in flow)\n
                               @param frames The numeric vector indicating the frames of this track that are 'valid'"

                               validFrames <<- points$frame[points$frame %in% frames]
                          },
                          getTrackSweep = function(amplitude=1, offset=mean(getSlot(slot='x')), validFramesOnly=FALSE, guess=NULL)
                          {
                               "This method is provided as a convenience. It calls 'getSweep' using
                               parameters that exist within the 'Track' and parent 'TrackList'
                               when available or the 'getSweep' defaults.\n
                               @param amplitude A numeric value\n
                               @param offset A numeric value\n
                               @param validFrames A boolean indicating whether to only get sweep values for valid frames only\n
                               @param guess A boolean indicating whether to guess appropriate parameters for this track"

                               if(validFramesOnly)
                               {
                                    frames <- validFrames
                               }
                               else
                               {
                                    frames <- points$frame
                               }
                               args <- list(sin=meta$sin, fi=meta$fi, ff=meta$ff, tAll=meta$tAll, phaseShift=meta$phaseShift, amplitude=amplitude, offset=offset, frames=frames, guess=guess)
                               args <- args[!isempty(args)]
                               return(do.call(getSweep, args))
                          },
                          sseTrack = function(amplitude=50, phaseShift=0, offset=0, validFramesOnly=FALSE)
                          {
                               "Calculates the sum square error (i.e., sse) between this track
                               and a sweep function with the given 'amplitude', 'phaseShift',
                               and 'offset'.\n
                               @param validFramesOnly A boolean to limit the calculation to just the 'validFrames' listed in this Track.\n
                               The additional sweep parameters of sin, tAll, fi, and ff are passed\n
                               from the track's parent 'TrackList'."

                               if(validFramesOnly)
                               {
                                    frames <- validFrames
                               }
                               else
                               {
                                    frames <- points$frame
                               }
                               args <- list(sin=meta$sin, fi=meta$fi, ff=meta$ff, tAll=meta$tAll, amplitude=amplitude, offset=offset, frames=frames)
                               args <- args[!isempty(args)]
                               predicted <- do.call(getSweep, args)
                               data <- object$points$vx                ### Explicitly fitting vx ###
                               indicesToGet <- which(points$frame %in% frames)
                               result <- sum((data[indicesToGet]-predicted$v)^2)
                               return(result)
                          }
                     )
)

#' Get the adjustable running window average of the data
#' @param i The index within 'frames' at which to calculate an average over a window centered at this location
#' @param frames The frames in this track
#' @param widths A vector of window widths appropriate for each frame in the 'frames' of this track
#' @param data The vector of data for which we will calculate the windowed averages
getAverage <- function(i, frames, widths, data)
{
     # Subtract 1 to represent the number of intervals instead of number of points to average
     width <- widths[i] - 1

     # calculate the index on the left of the interval
     leftIndex <- i - floor(width/2)
     if(leftIndex < 1)
     {
          leftIndex <- 1
     }

     # calculate the width to reach index on the right of the interval
     if((leftIndex+width) > length(frames))
     {
          width <- length(frames)-leftIndex
     }

     # return the mean of the data over the interval
     mean(data[leftIndex:(leftIndex+width)])
}

#' Get the derivative of a vector
#' @param x A numeric vector on which to calculate the derivative
#' @param t A numeric vecotor of times with which to determine dt for derivative calculations
getDerivative <- function(x, t)
{
     v <- numeric(0)
     for(i in 1:length(x))
     {
          v <- c(v, localDerivative(x, t, i))
     }
     return(v)
}

#' Get the local derivative around a point in a vector accounding for boundary scenarios at the start and end of the vector
#' @param x A numeric vector of data
#' @param t A numeric vector of time for calculating dt of the derivative
#' @param i A numeric value indicating the index in the x and t for which to calculate the local derivative
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

#' @title This is a three point interpolation of the derivative where the interpolated point
#' is the middle of the 3 points.
#'
#' @description This simplifies to the the three-point midpoint formula
#' when the time steps are equal but can handle when timesteps are unequal (i.e., the
#' time-step on either side of the 3 points is not equal)
#'
#' @param f0 numeric left function value
#' @param f1 numeric middle function value
#' @param f2 numeric right function value
#' @param x0 numeric left independent value
#' @param x1 numeric middle independent value
#' @param x2 numeric right independent value
#' @param xj numeric x value for which to evaluate the function
interpolateDerivative <- function(f0, f1, f2, x0, x1, x2, xj)
{
     term1 <- f0*((2*xj-x1-x2)/((x0-x1)*(x0-x2)))
     term2 <- f1*((2*xj-x0-x2)/((x1-x0)*(x1-x2)))
     term3 <- f2*((2*xj-x0-x1)/((x2-x0)*(x2-x1)))
     return(term1 + term2 + term3)
}
