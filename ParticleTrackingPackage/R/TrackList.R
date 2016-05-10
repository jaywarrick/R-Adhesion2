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
#' allTimes numeric vector The vector of times associate with 'allFrames'
#' allFrames numeric vector The vector containing a complete list of frame numbers existing in the data set
#' selectedFrames numeric vector The subset of 'allFrames' which are set as valid using 'setSelectedFrames'
TrackList <- setRefClass('TrackList',
					fields = list(tracks='list', meta='list'), #sin='logical', fi='numeric', ff='numeric', phaseShift='numeric', t0_Frame='numeric', timePerFrame='numeric', sweepDuration='numeric', allTimes='numeric', allFrames='numeric', selectedFrames='numeric'
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
							}
						},
						setSelectedFrames = function(selectedFrames)
						{
						     "Set which frames are 'selected'. Many functions have the ability to operate on only
						     the 'selected' frames of each track etc. This can be used to set which frames those are.\n
						     @param selectedFrames numeric vector of frame numbers that are 'valid'"

						     meta$selectedFrames <<- selectedFrames
						},
						getSelectedFrames = function()
						{
						     "Set which frames are 'selected'. Many functions have the ability to operate on only
						     the 'selected' frames of each track etc. This can be used to set which frames those are.
                                   If the selected frames have not been set, then the function returns the 'allFrames'
                                   value stored in the 'meta' list field.\n
						     @param selectedFrames numeric vector of frame numbers that are 'valid'"

						     if(is.null(meta$selectedFrames))
						     {
						          return(meta$allFrames)
						     }
						     else
						     {
						          return(meta$selectedFrames)
						     }
						},
						trackCount = function()
						{
							"Return the number of tracks in this TrackList"

							return(base::length(tracks))
						},
						plotTrackList = function(slotX='frame', slotY='x', fun=NULL, rel=FALSE, ...)
						{
							"Plot all the values of the provided 'slot' for all tracks on a single plot.\n
							@param slotX string Name of the column within the Track$pts to plot on the X-axis (e.g., x, y, t, frame, vx, vy and potentially vxs and vys smoothed velocities)\n
                                   @param slotY string Name of the column within the Track$pts to plot on the Y-axis (e.g., x, y, t, frame, vx, vy and potentially vxs and vys smoothed velocities)\n
							@param fun function A function to transform the slot data from each track before plotting (e.g., to scale between pixels per second velocity to microns per second velocity)\n
							@param rel boolean Whether to plot the slot relative to its mean"

							xRanges <- applyFun_Return(fun=function(track){track$range(slot=slotX)})
							xRanges <- matrix(unlist(xRanges), ncol=2, byrow=TRUE)
							yRanges <- applyFun_Return(fun=function(track){track$range(slot=slotY)})
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
							# print(args)
							first <- TRUE
							for(.track in tracks)
							{
								if(first)
								{
									do.call(.track$plotTrack, c(list(slotX=slotX, slotY=slotY, funY=fun, relY=rel, add=F, xlim=xlim, ylim=ylim, withTitle=FALSE, main='All Tracks'), args))
									first = FALSE
								}
								else
								{
									do.call(.track$plotTrack, c(list(slotX=slotX, slotY=slotY, funY=fun, relY=rel, add=T), args))
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
						addTrackPoint = function(id, frame, ...)
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
							track$addPoint(frame=frame, ...)
							setTrack(track)
						},
						applyFun_Return = function(fun=function(x){return(x$frameCount())}, ...)
						{
						     "Apply the given function to all tracks, summarizing the results in a list by Track id\n
							@param id numeric the ID number of the track of interest
						     @param ... additional args passed to fun"

						     temp <- lapply(tracks, FUN=fun, ...)
						     names(temp) <- names(tracks)
						     return(temp)
						},
						applyFun_Replace = function(fun, ...)
						{
							"Apply the provided function to alter each track, replacing the existing tracks with the altered tracks\n
							@param fun function The function to apply\n
							@param ... additional arguments to pass to fun"

							for(.track in tracks)
							{
								.track <- fun(.track, ...)
							}
						},
						applyFun_Void = function(funName, ...)
						{
						     "Apply the provided function to each track. Intended for 'void' functions that will produce
                                   and action such as adding the information for each track to a plot.\n
							@param fun function The function to apply\n
							@param ... additional arguments to pass to fun"

						     temp <- lapply(tracks, FUN=fun, ...)
						},
						callTrackFun = function(funName, ...)
						{
						     "Call a Track method/function on each track in the TrackList (e.g., plotting functions)\n
							@param funName string The name of the Track function to call on each Track\n
							@param ... additional arguments to pass to the called function"

						     tot <- trackCount()
						     myCount <- 0
						     for(.track in tracks)
						     {
						          myCount <- myCount + 1
						          cat("Calling", funName, "on track", myCount, "of", tot, "\n")
						          # Have to do eval(parse()) because track[[funName]] is NULL while track$parsedFunName is not NULL, don't know why
						          # Now that the function is loaded we can call it using the [[]] method
						          theCall <- paste(".track$'", funName, "'", sep="")
						          theFunc <- eval(parse(text=theCall))
						          if(is.null(theFunc))
						          {
						               stop(cat("Couldn't find function with name",funName))
						          }
						          do.call(theFunc, list(...))
						     }
						},
						updateFramesAndTimes = function(t0_Frame, timePerFrame)
						{
						     "Calculate the complete ordered list of frames that exist within the TrackList
						     and set that information in the metadata object as 'allFrames'.\n
						     Store the t0_Frame and timePerFrame information in the 'meta' list field of the TrackList.
						     The t0_Frame is the frame for which the time is equal to 0 and the amount
						     of time between each frame is the timePerFrame. From this information, the time
						     associated with each frame is determined as well as the time derivatives
						     of each measure.\n
						     Given the frames in the TrackList, t0_Frame, and timePerFrame, we can calculate
						     the times associated with each frame that exists within the TrackList
						     and set that information in the metadata object as 'allTimes'. This also
						     causes a call to 'updateTimes' on each track with the appropriate"

						     possibleFrames <- c()
						     for(.track in tracks)
						     {
						          possibleFrames <- c(possibleFrames, .track$getSlot('frame'))
						     }
						     meta$t0_Frame <<- t0_Frame
						     meta$timePerFrame <<- timePerFrame
						     meta$allFrames <<- sort(unique(possibleFrames))
							meta$allTimes <<- (meta$allFrames - meta$t0_Frame) * meta$timePerFrame
							callTrackFun('updateTimes', t0_Frame=t0_Frame, timePerFrame=timePerFrame)
						},
						calculateDerivatives = function(slots, withRespectTo='t', prefix='v')
						{
							"Call the 'calulateDerivatives' method on all Track objects in the TrackList"

							callTrackFun('calculateDerivatives', slots=slots, withRespectTo=withRespectTo, prefix=prefix)
						},
						calculateSmoothedData = function(windowWidths, slots, suffix='s')
						{
						     "Smooth the specified slots in each track using the specified windowWidths
                                   saving the calculations into each track using the slot name appended by the
                                   the supplied suffix.\n
                                   @param windowWidths data.frame(frame=numeric, width=numeric) containing an window width for each frame that exists in the list of Track objects or a single number that will be applied to all frames for all Track objects.\n
                                   @param slots character vector of all the columns in the 'pts' data.frame field for which smoothed versions of the data should be calculated."

						     for(.track in tracks)
						     {
						          .track$calculateSmoothedData(windowWidths=windowWidths, slots=slots, suffix=suffix)
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
						sortTracks = function(fun=function(x){return(x$frameCount())}, decreasing=TRUE, ...)
						{
							"Sort the track according to the function provided. The function should
							take a track as its first argument and return a sortable value for sorting
							The tracks are then sorted in the TrackList based on these values\n
                                   @param fun function The function to determine a value for sorting tracks
                                   [default=function(x){return(x$frameCount())} (thus sorting by length)]\n
							@param decreasing boolean whether to sort in increasing or decreasing order"

							ret <- applyFun_Return(fun=fun, ...)
							sorted <- sort(unlist(ret), index.return=T, decreasing=decreasing)
							tracks <<- tracks[sorted$ix]
						},
						getMatrix = function(slot='vx', selectedOnly=FALSE, rel=FALSE)
						{
							"Get a matrix of the tracklist data. Track id's are rows while
							frame numbers are columns. Use the matrix form of TrackList to do
							bulk operations such as determining the sse of all the data to a
							single curve.\n
							@param slot string The name of the column in Track$pts for which to make a matrix of data for\n
							@param selectedOnly boolean whether to make the matrix for valid frames only - defualt=FALSE"

							if(selectedOnly)
							{
							     frames <- meta$selectedFrames
							}
							else
							{
							     frames <- meta$allFrames
							}
							ids <- names(tracks)
							data <- matrix(NA, base::length(ids), base::length(frames), dimnames=list(id=names(tracks), frame=as.character(frames)))
							for(track in tracks)
							{
							     if(selectedOnly)
							     {
							          data[as.character(track$id), as.character(track$getSlot(slot='frame', rel=FALSE, selectedFrames=frames))] <- track$getSlot(slot=slot, rel=rel, selectedFrames=frames)
							     }
							     else
							     {
							          data[as.character(track$id), as.character(track$getSlot(slot='frame', rel=FALSE, selectedFrames=NULL))] <- track$getSlot(slot=slot, rel=rel, selectedFrames=NULL)
							     }
							}
							return(data)
						},
						calculateAvgPerFrameAcrossTracks = function(slot = 'vx', rel=FALSE, selectedOnly=FALSE)
						{
						     "Calculate the average value of a track slot such as x, y, vx, vxs, vy, and vys
						     return it as a vector."
                                   if(selectedOnly)
                                   {
                                        return(colMeans(getMatrix(slot=slot, rel=rel, selectedFrames=meta$selectedFrames)))
                                   }
						     else
						     {
						          return(colMeans(getMatrix(slot=slot, rel=rel, selectedFrames=meta$allFrames)))
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
								tracks[as.character(.track$id)] <<- .track$copy()
							}
						}
					)
)

#' Get the last element of a vector
#'
#' @param x vector
last <- function(x)
{
	return(x[numel(x)])
}

# calculateTrackCounts = function()
# {
#      "Return the number of tracks in the TrackList object"
#
#      trackMatrix <- getMatrix(slot='vx', selectedOnly=TRUE)
#      if(base::length(trackMatrix) == 0 || nrow(trackMatrix) == 0)
#      {
#           return(0)
#      }
#      else
#      {
#           ret <- list()
#           frames <- colnames(trackMatrix)
#           lastFrame <- last(frames)
#           return(sum(!is.na(trackMatrix[,lastFrame])))
#      }
# }
