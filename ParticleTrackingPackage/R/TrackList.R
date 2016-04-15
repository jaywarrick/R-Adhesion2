#' @import methods
#' @importFrom pracma numel isempty
NULL


#' Class representing a list of Track objects. This class is used to do track-wise
#' data manipulations (e.g., sweep fitting)
#' @field tracks list of Track objects
#' @field meta A list of metadata
#'
#' Metadata for log frequency sweep data includes:
#' sin boolean Whether the tracks are sinusoidal or triangular
#' fi numeric The initial frequency of the sweep
#' ff numeric The final frequency of the sweep
#' phaseShift numeric The phase shift of the sweep
#' t0_Frame numeric The frame associated with time=0 of the sweep
#' timePerFrame numeric The time per frame of the data in seconds
#' sweepDurtion numeric The time duration of the sweep in seconds
#' tAll numeric vector The vector of times associate with 'allFrames'
#' allFrames numeric vector The vector containing a complete list of frame numbers existing in the data set
#' validFrames numeric vector The subset of 'allFrames' which are set as valid using 'setValidFrames'
TrackList <- setRefClass('TrackList',
					fields = list(tracks='list', meta='list'), #sin='logical', fi='numeric', ff='numeric', phaseShift='numeric', t0_Frame='numeric', timePerFrame='numeric', sweepDuration='numeric', tAll='numeric', allFrames='numeric', validFrames='numeric'
					methods = list(
						initializeWithJEXROIFile = function(file=NULL, t0_Frame, timePerFrame)
						{
							"Initialize the TrackList with a JEXROI file that has been 'tracked' already (i.e the
							id's of points have already been linked through time using a function like LAP Tracker)"

							require(foreign)
							tracksFile <- read.arff(file)
							tracksFile2 <- reorganize(tracksFile, measurementCols='Metadata')
							for(row in 1:nrow(tracksFile2))
							{
								id <- tracksFile2[row,'Track']
								start <- tracksFile2[row,'polygonPts']
								pattern <- tracksFile2[row,'patternPts']
								newTrack <- new('Track')
								newTrack$initializeWithTrackROI(id=id, start=start, pattern=pattern, t0_Frame=t0_Frame, timePerFrame=timePerFrame)
								setTrack(newTrack)
								frameLimits <- range()
							}
						},
						setStandardMeta = function(t0_Frame=t0_Frame, timePerFrame=timePerFrame)
						{
						     "Set the 'standard' metadata for a TrackList object which typically consists
						     of a frame for which the time is equal to 0 (t0_Frame) and the amount
						     of time between each frame (timePerFrame)."
							meta <<- list()
							meta$t0_Frame <<- t0_Frame
							meta$timePerFrame <<- timePerFrame
							calculateAllFrames()
							calculateTAll()
							calculateVelocities()
							callTrackFun('setMeta', meta)
						},
						setOscillatoryMeta = function(sin, fi, ff, t0_Frame, timePerFrame, sweepDuration)
						{
							"Set the sweep parameters used to predict particle motion during a log frequency sweep\n
							@param sin boolean Whether the tracks are sinusoidal or triangular\n
							@param fi numeric The initial frequency of the sweep\n
							@param ff numeric The final frequency of the sweep\n
							@param t0_Frame numeric The frame associated with time=0 of the sweep\n
							@param timePerFrame numeric The time per frame of the data in seconds\n
							@param sweepDurtion numeric The time duration of the sweep in seconds\n
							The phase shift is not set here as that is determined using the getBulkPhaseShift method"

							meta <<- list()
							meta$sin <<- sin
							meta$fi <<- fi
							meta$ff <<- ff
							meta$t0_Frame <<- t0_Frame
							meta$timePerFrame <<- timePerFrame
							meta$sweepDuration <<- sweepDuration
							calculateAllFrames()
							calculateTAll()
							calculateVelocities()
							callTrackFun('setMeta', meta)
						},
						length = function()
						{
							"Return the number of tracks in this TrackList"

							return(base::length(tracks))
						},
						plotTrackList = function(slot='vx', fun=NULL, rel=FALSE, ...)
						{
							"Plot all the values of the provided 'slot' for all tracks on a single plot.\n
							@param slot string Name of the column within the Track$points to plot on the Y-axis (options are x, y, t, frame, vx, vy and potentially vxs and vys smoothed velocities, X-axis is always 't' with this function)\n
							@param fun function A function to transform the slot data from each track before plotting (e.g., to scale between pixels per second velocity to microns per second velocity)\n
							@param rel boolean Whether to plot the slot relative to its mean"

							xRanges <- getProp(fun=function(track){track$range('t')})
							xRanges <- matrix(unlist(xRanges), ncol=2, byrow=TRUE)
							yRanges <- getProp(fun=function(track){track$range(slot)})
							yRanges <- matrix(unlist(yRanges), ncol=2, byrow=TRUE)

							Xmin <- min(xRanges[,1], na.rm=TRUE)
							Xmax <- max(xRanges[,2], na.rm=TRUE)
							Ymin <- min(yRanges[,1], na.rm=TRUE)
							Ymax <- max(yRanges[,2], na.rm=TRUE)

							args <- list(...)
							if(is.null(args$xlim))
							{
								xlim <- c(Xmin, Xmax)
							}
							else
							{
								xlim <- args$xlim
								args$xlim <- NULL
							}
							if(is.null(args$ylim))
							{
								ylim <- c(Ymin, Ymax)
							}
							else
							{
								ylim <- args$ylim
								args$ylim <- NULL
							}
							print(args)
							first <- TRUE
							for(track in tracks)
							{
								if(first)
								{
									do.call(track$plotTrack, c(list(slotY=slot, funY=fun, relY=rel, add=F, xlim=xlim, ylim=ylim, withTitle=FALSE, main='All Tracks'), args))
									first = FALSE
								}
								else
								{
									do.call(track$plotTrack, c(list(slotY=slot, funY=fun, relY=rel, add=T), args))
								}
							}
						},
						getTrack = function(id)
						{
							"Get the track associated with the provided ID\n
							@param id numeric the ID number of the track of interest"

							return(tracks[[as.character(id)]])
						},
						setTrack = function(newTrack)
						{
							"Set the given track within the TrackList. Any track with the same ID
							is replaced\n
							@param newTrack Track the track to put into the list/set"

							# newTrack$.parent <- .self // Only need to setMeta once after all tracks are added
							tracks[[as.character(newTrack$id)]] <<- newTrack
						},
						removeTrack = function(id)
						{
							"Remove the track associated with the provided ID\n
							@param id numeric the ID number of the track of interest"

							tracks[[as.character(id)]] <<- NULL
						},
						getProp = function(fun=function(x){return(x$length())}, ...)
						{
							"Apply the given function to all tracks, summarizing the results in a list\n
							@param id numeric the ID number of the track of interest"

							temp <- lapply(tracks, FUN=fun, ...)
							names(temp) <- names(tracks)
							return(temp)
							#    ret <- list()
							#    for(.track in tracks)
							#    {
							#         ret[[as.character(.track$id)]] <- fun(.track, ...)
							#    }
							#    return(ret)
						},
						addTrackPoint = function(id, x, y, frame)
						{
							"Add the given point (x, y, frame) to the track associated with the provided ID
							This is useful for building a TrackList one point at a time instead via initialization with a file\n
							@param id numeric The ID number of the track of interest\n
							@param x numeric The x position of the point\n
							@param y numeric The y position of the point\n
							@param frame numeric The frame number of the point"

							track <- getTrack(id)
							if(is.null(track))
							{
								track <- new('Track')
								track$id <- id
								# track$.parent <- .self # setTrack resets .parent
							}
							track$addPoint(x=x, y=y, t=(frame-meta$t0_Frame)*meta$timePerFrame, frame=frame)
							setTrack(track)
						},
						calculateAllFrames = function()
						{
							"Calculate the complete ordered list of frames that exist within the TrackList
						     and set that information in the metadata object as 'allFrames'"

							possibleFrames <- c()
							for(.track in tracks)
							{
								possibleFrames <- c(possibleFrames, .track$getSlot('frame'))
							}
							meta$allFrames <<- sort(unique(possibleFrames))
						},
						calculateTAll = function()
						{
							"Calculate the times associated with each frame that exists within the TrackList
						     and set that information in the metadata object as 'tAll'"

							print('Calculating times for all frames in TrackList')
							meta$tAll <<- (meta$allFrames - meta$t0_Frame) * meta$timePerFrame
						},
						calculateVelocities = function()
						{
							"Call the 'calulateVelocities' method on all Track objects in the TrackList
						     This results in calculating the point-to-point instantaneous velocities of each track (i.e., no smoothing)"

							callTrackFun('calculateVelocities')
						},
						smoothVelocities = function(fit, dist, maxWidth)
						{
							"Call the 'smoothVelocities' method on all the Track objects in the TrackList\n
							@param fit The fit results returned from 'getBulkPhaseShift' fitting function\n
							@param dist numeric The number of pixels of motion that should be observed before estimating the velocity"

							tempAllWidths <- getWindowWidths(fit=fit, dist=dist, maxWidth=maxWidth)
							callTrackFun('smoothVelocities', allWidths=tempAllWidths, allFrames=meta$allFrames)
						},
						calculateValidTimes = function()
						{
							"Calculate the times associated with the valid frames of this TrackList
						     and return that information directly"

							return(meta$tAll[meta$allFrames %in% meta$validFrames])
						},
						applyToTracks = function(fun, ...)
						{
							"Apply the provided function to alter each track, replacing the exisitng tracks with the altered tracks\n
							@param fun function The function to apply\n
							@param ... additional arguments to pass to fun"

							for(track in tracks)
							{
								track <- fun(track, ...)
							}
						},
						callTrackFun = function(funName, ...)
						{
							"Call a Track method/function on each track in the TrackList\n
							@param funName string The name of the Track function to call on each Track\n
							@param ... additional arguments to pass to the called function"

							tot <- length()
							myCount <- 0
							for(track in tracks)
							{
							     myCount <- myCount + 1
								cat("Calling", funName, "on track", myCount, "of", tot, "\n")
								# Have to do eval(parse()) because track[[funName]] is NULL while track$parsedFunName is not NULL, don't know why
								# Now that the function is loaded we can call it using the [[]] method
								theCall <- paste("track$'", funName, "'", sep="")
								theFunc <- eval(parse(text=theCall))
								if(is.null(theFunc))
								{
									stop(cat("Couldn't find function with name",funName))
								}
								do.call(theFunc, list(...))
							}
						},
						filterTracks = function(fun, ...)
						{
							"Filter tracks using the provided function that returns T or F when provided a track as the first argument\n
							@param fun function The filtering function that returns T or F for each Track individually\n
							@param ... additional arguments to fun"

							tracks <<- Filter(function(x){fun(x,...)}, tracks)
							if(base::length(tracks) == 0)
							{
								message("No tracks fit filter, resulting tracklist is of length 0!")
							}
						},
						sortTracks = function(fun=function(x){return(x$length())}, decreasing=TRUE, ...)
						{
							"Sort the track according to the function provided. The function should
							take a track as its first argument and return a sortable value for sorting
							The tracks are then sorted in the TrackList based on these values\n
							@param fun function The function to determine a value for sorting tracks - default=function(x){return(x$length())} (thus sorting by length)\n
							@param decreasing boolean whether to sort in increasing or decreasing order"

							ret <- getProp(fun=fun, ...)
							sorted <- sort(unlist(ret), index.return=T, decreasing=decreasing)
							tracks <<- tracks[sorted$ix]
						},
						getMatrix = function(slot='vx', validOnly=FALSE, rel=FALSE)
						{
							"Get a matrix of the tracklist data. Track id's are rows while
							frame numbers are columns. Use the matrix form of TrackList to do
							bulk operations such as determining the sse of all the data to a
							single curve.\n
							@param slot string The name of the column in Track$points for which to make a matrix of data for\n
							@param validOnly boolean whether to make the matrix for valid frames only - defualt=FALSE"

							if(!validOnly)
							{
								frames <- meta$allFrames
							}
							else
							{
								frames <- meta$validFrames
							}
							ids <- names(tracks)
							data <- matrix(NA, base::length(ids), base::length(frames), dimnames=list(id=names(tracks), frame=as.character(frames)))
							for(track in tracks)
							{
								data[as.character(track$id), as.character(track$getSlot(slot='frame', rel=rel, validOnly=validOnly))] <- track$getSlot(slot=slot, rel=FALSE, validOnly=validOnly)
							}
							return(data)
						},
						getWindowWidths = function(fit, dist, maxWidth)
						{
							"Calculate the appropriate window widths for calculating velocity\n
							@param fit The fit results returned by 'getBulkPhaseShift'\n
							@param dist The number of pixels that we would like to see (if possible) the particle/cell move before estimating velocity\n
							@param maxWidth numeric The maximum number of FRAMES that should be averaged together to estimate velocity (i.e., averaging is good up to a point where it just increases compuation time)\n\n
							Underlying Calculations: v represents the velocity calculated over a single period. Thus, the cell should travel
							4*A units of distance per period. This is used to estimate the velocity of the cell
							using the fit. Therefore, the expected distance the cell will travel between frames is
							v*dt where dt is the time between the frames of interest. Thus, if we would like to
							quantify velocity over a travel distance of 10 pixels (i.e., dist=10), then we can take
							dist/v*dt to get the number of frames (i.e., dt's) needed to cover 10 pixels. However, if maxWidth
							is exceeded, maxWidth is returned. dt could be variable over time depending on frame sampling rates etc.
							approximate a dt for each index based on diff's. Repeat last value at end to create a dt
							vector that is the same as the tracks$tAll vector"
							dt <- diff(meta$tAll)
							dt <- c(dt, last(dt))
							v <- 4*fit$par[['amplitude']]*(meta$fi*(meta$ff/meta$fi)^(meta$tAll/last(meta$tAll)))
							widths <- ceiling(dist/(v*dt))
							widths[widths > maxWidth] <- maxWidth
							return(widths)
						},
						calculateTrackAverage = function(slot = 'vx', rel=FALSE, validOnly=FALSE)
						{
						     "Calculate the average value of a track slot such as x, y, vx, vxs, vy, and vys
						     return it as a vector."

						     return(colMeans(getMatrix(slot=slot, rel=rel, validOnly=validOnly)))
						},
						calculateValidFrames = function(fit, validStart=0.01, validEnd=0.99)
						{
							"Takes the fit and determines where the cells switch directions\n
                                   'validStart' and 'validEnd' represent precentages. In other words, after the switch in direction,
							valid points start at 'validStart' % of the way to the next switch in direction while 'validEnd'
							occurs 'validEnd' % of the way to the that same next switch in direction.\n

							This filter is good for removing inaccurate values of the velocity when using a triangle waveform because
							the estimate of velocity will be artificially be lower due to sample aliasing of particle motion.
							Possible values are 0, 1, 2, and 3 for inflectionPtsToFilter. Multiple at once can be provided.\n

							0 Represents the upward zero-crossing of the position and max positive velocity\n
							1 Represents the max position and downward zero-crossing of the velocity\n
							2 Represents the downward zero-crossing of the position and min (i.e most negative) velocity\n
							4 Represents the min (i.e., most negative) position and upward zero-crossing of the velocity\n

                                   'inflectionPtsToFilter' refers to the first point of the numbered sections (i.e., 1 refers to the max point)
							Thus if 1 is included in 'inflectionPtsToFilter', points adjacent to this point will be considered for wobble and filtered.\n

                                   We use 'inflectionPtsToFilter' of c(1,3) to mark the points at which the flow switches direction"

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

							sweep <- getSweep(amplitude=fit$par['amplitude'], phaseShift=fit$par['phaseShift'], offset=0, sin=meta$sin, fi=meta$fi, ff=meta$ff, sweepDuration=meta$sweepDuration, t=meta$tAll, guess=NULL)
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
							validFrames0 <- meta$allFrames
							validFrames0 <- validFrames0[-indicesToRemove]

							validFrames1 <- numeric(0)
							for(tIndex in 1:base::length(meta$tAll))
							{
								# Get the nearest inflection at or beyond this time

								infIndex <- which.max((sweep$inflections >= meta$tAll[tIndex]) & inflectionsToAddress)
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
								if( (meta$tAll[tIndex] >= (infT1 + validStart*dInfT)) & (meta$tAll[tIndex] <= (infT1 + validEnd*dInfT)) )
								{
									# If it is within the startValid and endValid bounds, add it to the list of the valid frames
									validFrames1 <- c(validFrames1, meta$allFrames[tIndex]) # (tIndex-1) = frame because frames are indices that start at 0
								}
							}

							meta$validFrames <<- sort(intersect(validFrames1, validFrames0))

							for(track in tracks)
							{
								track$setValidFrames(meta$validFrames)
							}
						},
						getPercentAdhered = function(velocityThreshold=3)
						{
							"Calculate the percent adhered at each timepoint using the provided velocity threshold
                                   @param velocityThrehsold numeric The pixels per second below which (exclusive) a cell is considered adhered - default=3 [pixels/second]\n
							@return Return a dataframe with columns of 'time' and 'percentAdhered'"

							trackMatrix <- getMatrix(slot='vx', validOnly=TRUE, rel=FALSE)
							ret <- list()
							cellCount <- getCellCount()
							frames <- colnames(trackMatrix)
							for(frame in frames)
							{
								velocities <- trackMatrix[,frame]
								velocities <- abs(velocities[!is.na(velocities)])
								if(!isempty(velocities))
								{
									adhered <- sum(velocities < velocityThreshold)/cellCount
									ret[[frame]] <- adhered
								}
								else
								{
									ret[[frame]] <- 0
								}
							}
							percents <- 100*as.numeric(ret)
							times <- meta$tAll[meta$allFrames %in% frames]
							return(data.frame(time=times, percentAdhered=percents))
						},
						getCellCount = function()
						{
						     "Return the number cells/tracks in the TrackList object"

							trackMatrix <- getMatrix(slot='vx', validOnly=TRUE)
							if(base::length(trackMatrix) == 0 || nrow(trackMatrix) == 0)
							{
								return(0)
							}
							else
							{
								ret <- list()
								frames <- colnames(trackMatrix)
								lastFrame <- last(frames)
								return(sum(!is.na(trackMatrix[,lastFrame])))
							}
						},
						save = function(objectName, file) {
							"Save the current object on the file in R external object format."

							assign(objectName, .self)
							base::save(list=c(objectName), file = file)
						},
						refreshTracks = function()
						{
						     "Refresh the TrackList object. This causes a 'copy' of the object
						     which then creates a new object with the latest definition of the
						     TrackList class, thus updating the old version to a new version
						     if possible"

							for(.track in tracks)
							{
								tracks[.track$id] <<- .track$copy()
							}
						}
					)
)

#' Get the last element of a vector
#'
#' @param x vector
#'
#' @export
last <- function(x)
{
	return(x[numel(x)])
}
