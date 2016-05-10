#' @import methods
NULL

#' PointSetList
#'
#' A class for storing and working with list of PointSet objects. This class is the
#' class that actually does the point tracking over time.
#'
#' Note that the data is stored internally as a list called 'pointSetList'. This pointSetList
#' is a list object that stores each pointSet with the list element name equal to the frame
#' number. Index is used to refer to the index within that list, not the frame number / name.
#'
#' @field pointSetList list of PointSet objects
#' @field maxID numeric value indicating the maximum point ID encountered in all the PointSets in the PointSetList
#' @field trackingStartFrame numeric value indicating the frame at which cells were tracked backward from (i.e., usually a late frame)
#' @field trackingEndFrame numeric value indicating the frame at which cells were tracked back to (i.e., usually an early frame)
#' @field meta list for storing metadata about the PointSetList
PointSetList <- setRefClass('PointSetList',
                          fields = list(pointSetList='list', maxID='numeric', trackingStartFrame='numeric', trackingEndFrame='numeric', meta='list'),
                          methods = list(
                               initializeWithDataFrame = function(data, frameNumberCol, idCol=NULL)
                               {
                                    "Initialize with a data frame where the rows describe the points. The data.frame must have a column to indicate the frame number of the points. IdCol is optional, used if provided, otherwise set arbitrarily incrementally for each frame."

                                    for(frame in unique(data[,frameNumberCol]))
                                    {
                                         newPointSet <- new('PointSet')
                                         newPointSet$initializeWithROI(frame=pointSetListFile2[r,frameDimName], ptString=pointSetListFile2[r,'polygonPts'])
                                         temp[[as.character(pointSetListFile2[r,frameDimName])]] <- newPointSet
                                    }
                               },
                               initializeWithJEXROIFile = function(path, frameDimName='Time')
                               {
                                    "Initialize with a file (eg., .jxd ROI file from JEX) that stores all the point sets as semicolon separated comma separated point information (see PointSet documentation)."

                                    pointSetListFile <- read.arff(path)
                                    initializeWithROIDataFrame(roiTable=pointSetListFile, frameDimName=frameDimName)
                               },
                               initializeWithROIDataFrame = function(roiTable, frameDimName='Time')
                               {
                                    "Initialize with a data.frame of an ROI (eg., .jxd ROI file from JEX) that stores all the point sets as semicolon separated comma separated point information (see PointSet documentation)."

                                    meta$frameDimName <<- frameDimName

                                    # require(foreign)
                                    if(!is.null(roiTable))
                                    {
                                         pointSetListFile2 <- reorganize(roiTable, measurementCols='Metadata')
                                         pointSetListFile2[,frameDimName] <- as.numeric(as.character(pointSetListFile2[,frameDimName]))
                                         temp <- list()
                                         for(r in 1:nrow(pointSetListFile2))
                                         {
                                              newPointSet <- new('PointSet')
                                              newPointSet$initializeWithROI(frame=pointSetListFile2[r,frameDimName], ptString=pointSetListFile2[r,'polygonPts'])
                                              temp[[as.character(pointSetListFile2[r,frameDimName])]] <- newPointSet
                                         }
                                         temp <- temp[order(as.numeric(names(temp)), decreasing=FALSE)]
                                         pointSetList <<- temp
                                    }
                                    else
                                    {
                                         pointSetList <<- list()
                                    }

                                    maxID <<- -1 # Used for appending pointSetList
                               },
                               getIndexOfFrame = function(frame)
                               {
                                    "Get the index of a frame with a particular number"

                                    return(which(names(pointSetList)==as.character(frame)))
                               },
                               getFrameOfIndex = function(index)
                               {
                                    "Get the index of a frame with a particular number"

                                    return(getPointSet(index=index)$frame)
                               },
                               getPointSet = function(index=NULL, frame=NULL)
                               {
                                    "Get the PointSet at a particular frame or index. If both frame and index
                                    are specified, then NULL is returned. If no PointSet exists for the
                                    specified frame, NULL is returned. If index is out of bounds, NULL is returned."

                                    if(is.null(index) && is.null(frame))
                                    {
                                         return(NULL)
                                    }
                                    if(!is.null(frame))
                                    {
                                         return(pointSetList[[as.character(frame)]])
                                    }
                                    if(!is.null(index))
                                    {
                                         if(index < 1 | index > base::length(pointSetList))
                                         {
                                                 return(NULL)
                                         }
                                         return(pointSetList[[index]])
                                    }
                               },
                               getPreviousPointSet = function(frame)
                               {
                                    "Get the pointSetList that preceeds the given frame (good for back-tracking)"

                                    toGet <- getIndexOfFrame(frame)
                                    return(getPointSet(index=toGet-1))
                               },
                               getNextPointSet = function(frame)
                               {
                                    "Get the pointSetList that preceeds the given frame (good for forward-tracking)"

                                    toGet <- getIndexOfFrame(frame)
                                    return(getPointSet(index=toGet+1))
                               },
                               setPointSet = function(newPointSet)
                               {
                                    "Set this PointSet into the PointSetList using the frame number
                                    from the provided point set to put it into the list. This object
                                    is set by reference so changes to newPointSet are reflected in
                                    the object stored in this PointSetList."

                                    pointSetList[[as.character(newPointSet$frame)]] <<- newPointSet
                               },
                               pointSetCount = function()
                               {
                                    "Return the number of sets of pointSetList in this PointSetList object"

                                    return(base::length(pointSetList))
                               },
                               track = function(trackingStart='first', trackingEnd='last', trim=T, assignmentFunction, ...)
                               {
                                    "Perform point tracking going from start to end. When supplied as a numeric value, start and/or end
                                    are interpreted as indicies of the PointSetList, when supplied as a character value, they represent
                                    frame values. Thus, you can specify either as a frame number or index. One can be numeric and the
                                    character and they will be interppreted individually. One can also specify 'first' or 'last' to
                                    indicate the first or last index/frame in the PointSetList. Tracking can be done forwards or backwards.
                                    The trim parameter determines whether the PointSetList should be trimmed after tracking to include
                                    only those frames which were tracked. Trimming can help avoid confusion in point id's later on but
                                    is not technically necessary.

                                    Supply a tracking function that takes two PointSet objects and a variable length list of parameters
                                    for tracking. The two PointSet args should be ordered appropriately. The first arg ('before') will be linked to
                                    the second arg ('after'); however, that linking might be occuring forward or backward in time depending upon
                                    the start and end frame arguments. The function should return at least two objects (list(links, newTracks)).

                                    links = data.frame(id0=numeric, id1=numeric) where id0 column holds the ids from the 'before' PointSet
                                    object which are linked to the ids of the 'after' PointSet object, column id1.

                                    newTracks = numeric vector; ids in the 'after' PointSet object which are not linked to any ids in the
                                    'before' PointSet, thus representing the start of a new track.

                                    An example call to this function would be...

                                    psl$track(start='first', end='last', trim=T, directionalLinearAssignment, digits=1, maxDist=150, direction=c(1,0), perpendicularPenaltyFactor=1, uniformityDistThresh=-1)"

                                    # Check the arguments
                                    if(is.numeric(trackingStart))
                                    {
                                         trackingStart <- getFrameOfIndex(trackingStart)
                                    }
                                    if(is.numeric(trackingEnd))
                                    {
                                         trackingEnd <- getFrameOfIndex(trackingEnd)
                                    }
                                    if(trackingStart == 'first')
                                    {
                                         trackingStart <- getFrameOfIndex(1)
                                    }
                                    if(trackingStart == 'last')
                                    {
                                         trackingStart <- getFrameOfIndex(pointSetCount())
                                    }
                                    if(trackingEnd == 'first')
                                    {
                                         trackingEnd <- getFrameOfIndex(1)
                                    }
                                    if(trackingEnd == 'last')
                                    {
                                         trackingEnd <- getFrameOfIndex(pointSetCount())
                                    }
                                    startFrame <- trackingStart
                                    endFrame <- trackingEnd

                                    # trackingStart and trackingEnd have been interpreted/converted to startFrame and endFrame
                                    if(startFrame == endFrame)
                                    {
                                         # Then do nothing
                                         stop("Start and end frame are the same. Aborting tracking.")
                                    }
                                    if(is.null(getIndexOfFrame(startFrame)) || is.null(getIndexOfFrame(endFrame)))
                                    {
                                         stop("Start frame or end frame is not a valid frame in the pointSetList list.")
                                    }
                                    if(getIndexOfFrame(startFrame) > pointSetCount() | getIndexOfFrame(startFrame) < 0 | getIndexOfFrame(endFrame) > pointSetCount() | endFrame >= startFrame | getIndexOfFrame(endFrame) < 0)
                                    {
                                         stop("Start frame or end frame out of bounds. Start frame must be >=0 and > (yes >) end frame because we work backward. End frame must be < last frame of data set.")
                                    }

                                    frames <- startFrame:endFrame

                                    # Initialize the loop
                                    trackingStartFrame <<- startFrame
                                    trackingEndFrame <<- endFrame
                                    currentPointSet <- getPointSet(frame=startFrame)
                                    if(startFrame < endFrame)
                                    {
                                         nextPointSet <- getNextPointSet(startFrame)
                                    }
                                    else
                                    {
                                         nextPointSet <- getPreviousPointSet(startFrame)
                                    }
                                    maxID <<- max(currentPointSet$pts$id)

                                    # Loop over the specified frames and track
                                    while(!is.null(nextPointSet) && nextPointSet$frame %in% frames)
                                    {
                                         cat("Linking frame: ", currentPointSet$frame, " to ", nextPointSet$frame, "\n", sep="")

                                         assignments <- assignmentFunction(currentPointSet, nextPointSet, ...)
                                         # print(assignments)
                                         newPointSet <- new('PointSet')
                                         newPointSet$frame <- nextPointSet$frame
                                         if(nrow(assignments$links) > 0)
                                         {
                                              # Add linked points
                                              for(i in 1:nrow(assignments$links))
                                              {
                                                   p0 <- currentPointSet$getPoint(id=assignments$links$id0[i])
                                                   p1 <- nextPointSet$getPoint(id=assignments$links$id1[i])
                                                   newPointSet$addPoint(id=assignments$links$id0[i], p1[,!names(p0) %in% 'id'])
                                              }
                                         }
                                         if(base::length(assignments$id1_newTracks) > 0)
                                         {
                                              # Add points that represent the start of new tracks
                                              for(i in 1:base::length(assignments$id1_newTracks))
                                              {
                                                   maxID <<- maxID + 1
                                                   # print(paste0("maxID = ", maxID))
                                                   p1 <- nextPointSet$getPoint(id=assignments$id1_newTracks[i])
                                                   newPointSet$addPoint(id=maxID, p1[,!names(p0) %in% 'id'])
                                              }
                                         }

                                         # Overwrite the PointSet with the newPointSet that contains the updated
                                         # ids that form links between PointSets
                                         setPointSet(newPointSet)

                                         # Increment the loop
                                         currentPointSet <- getPointSet(frame=newPointSet$frame)
                                         if(startFrame < endFrame)
                                         {
                                              nextPointSet <- getNextPointSet(frame=newPointSet$frame)
                                         }
                                         else
                                         {
                                              nextPointSet <- getPreviousPointSet(frame=newPointSet$frame)
                                         }
                                    }
                                    cat("Done tracking.")

                                    if(trim)
                                    {
                                         cat("Trimming Data to start and end frame")
                                         allFrameNames <- as.numeric(names(pointSetList))
                                         indicesToRemove <- which(!(allFrameNames %in% endFrame:startFrame))
                                         pointSetList[indicesToRemove] <<- NULL
                                    }
                               },
                               getTrackList = function(exportStartFrame=0, exportEndFrame=-1)
                               {
                                    "Create a TrackList object from this PointSetList object (essentially a list indexed/grouped by cell id instead of time)\n
                                    @param exportStartFrame numeric value with two numbers lowest frame number that can be included in the exported TrackList. If the
                                    value is larger than the exportEndFrame, the smallest tracked frame, as recorded during the most recent call to 'track',
                                    will be used. Default is 0.\n
                                    @param exportEndFrame numeric value indicating the highest frame number that can be included in the exported TrackList. If the value is
                                    < exportStartFrame, the maximum tracked frame, as recorded during the most recent call to 'track', will be used instead. (default is -1)"

                                    if(exportEndFrame < exportStartFrame)
                                    {
                                         exportFrames <- names(pointSetList)[names(pointSetList) %in% as.character(trackingStartFrame:trackingEndFrame)]
                                    }
                                    else
                                    {
                                         exportFrames <- names(pointSetList)[names(pointSetList) %in% as.character(exportStartFrame:exportEndFrame)]
                                    }
                                    trackList <- new('TrackList')
                                    myCount = 0
                                    total = base::length(exportFrames)
                                    for(f in exportFrames)
                                    {
                                         .pointSet <- getPointSet(frame=f)
                                         if(.pointSet$pointCount() > 0)
                                         {
                                              for(i in 1:.pointSet$pointCount())
                                              {
                                                   data <- .pointSet$pts[i,]
                                                   trackList$addTrackPoint(id=data$id, frame=.pointSet$frame, x=data$x, y=data$y)
                                              }
                                              myCount <- myCount + 1
                                              cat("Generating TrackList: ", round(100*myCount/total, digits=2), "%\n", sep="")
                                         }
                                    }
                                    temp <- applyFun_Return(fun=function(.pointSet){return(.pointSet$pointCount())})
                                    temp <- data.frame(frame=as.numeric(as.character(names(temp))), val=as.vector(unlist(temp)))
                                    trackList$meta$pointCounts <- temp
                                    return(trackList)
                               },
                               applyFun_Return = function(fun=function(.pointSet){return(base::range(.pointSet$pts$x))}, ...)
                               {
                                    "Apply the given function to all PointSets, summarizing the results in a list by frame number\n
							@param id numeric the ID number of the track of interest
						     @param ... additional args passed to fun"

                                    temp <- lapply(pointSetList, FUN=fun, ...)
                                    names(temp) <- names(pointSetList)
                                    return(temp)
                               },
                               applyFun_Replace = function(fun, ...)
                               {
                                    "Apply the provided function to alter each track, replacing the existing tracks with the altered tracks\n
                                    @param fun function The function to apply\n
                                    @param ... additional arguments to pass to fun"

                                    for(.pointSet in pointSetList)
                                    {
                                         .pointSet <- fun(.pointSet, ...)
                                    }
                               },
                               applyFun_Void = function(funName, ...)
                               {
                                    "Apply the provided function to each track. Intended for 'void' functions that will produce
                                    and action such as adding the information for each track to a plot.\n
                                    @param fun function The function to apply\n
                                    @param ... additional arguments to pass to fun"

                                    temp <- lapply(pointSetList, FUN=fun, ...)
                               },
                               callPointSetFun = function(funName, ...)
                               {
                                    "Call a Track method/function on each track in the TrackList (e.g., plotting functions)\n
							@param funName string The name of the Track function to call on each Track\n
							@param ... additional arguments to pass to the called function"

                                    tot <- pointSetCount()
                                    myCount <- 0
                                    for(.pointSet in pointSetList)
                                    {
                                         myCount <- myCount + 1
                                         cat("Calling", funName, "on PointSet", myCount, "of", tot, "\n")
                                         # Have to do eval(parse()) because track[[funName]] is NULL while track$parsedFunName is not NULL, don't know why
                                         # Now that the function is loaded we can call it using the [[]] method
                                         theCall <- paste(".pointSet$'", funName, "'", sep="")
                                         theFunc <- eval(parse(text=theCall))
                                         if(is.null(theFunc))
                                         {
                                              stop(cat("Couldn't find function with name",funName))
                                         }
                                         do.call(theFunc, list(...))
                                    }
                               },
                               generatePointSetPlots = function(path=NULL, baseName='Frame_')
                               {
                                    "Produce plots of the pointSetList locations over time and output them to a folder"

                                    dir.create(path, recursive=TRUE)
                                    wd <- getwd()
                                    setwd(path)
                                    xlim <- base::range(unlist(applyFun_Return(fun=function(.pointSet){return(base::range(.pointSet$pts$x))})))
                                    ylim <- base::range(unlist(applyFun_Return(fun=function(.pointSet){return(base::range(.pointSet$pts$y))})))

                                    makeThePlots <- function()
                                    {
                                         for(.pointSet in pointSetList)
                                         {
                                              cat("Plotting Frame = ", .pointSet$frame, "\n", sep="")
                                              pdf(file=paste(path, '/', baseName, .pointSet$frame, '.pdf', sep=''), width=9, height=4)
                                              .pointSet$plotPointSet(IDs=TRUE, xlim=xlim, ylim=ylim, main=as.character(.pointSet$frame))
                                              dev.off()
                                         }
                                    }

                                    tryCatch(expr=makeThePlots(), finally=setwd(wd))
                               },
                               saveROI = function(file)
                               {
                                    "Export this PointSetList object as an ROI object for JEX (https://github.com/jaywarrick/JEX)"
                                    # print(as.factor(as.numeric(names(pointSetList))))
                                    temp1 <- data.frame(Time=as.numeric(names(pointSetList)), Metadata=rep(as.factor(c('Type','patternPts','polygonPts')), each=base::length(names(pointSetList))))
                                    temp1$Value[temp1$Metadata == 'Type'] <- '5'
                                    temp1$Value[temp1$Metadata == 'patternPts'] <- '0,0,0'
                                    for(.pointSet in pointSetList)
                                    {
                                         toSave <- ''
                                         jn <- nrow(.pointSet$pts)
                                         if(jn > 1)
                                         {
                                              for(j in 1:(jn-1))
                                              {
                                                   toSave <- paste0(toSave, paste0(.pointSet$pts$x[j], ",", .pointSet$pts$y[j], ",", as.character(.pointSet$pts$id[j]), ";"))
                                              }
                                         }
                                         if(jn > 0)
                                         {
                                              # use of "as.character" on id is important to get the right cell label because the factor levels are per pointSetList table instead of spanning pointSetList tables
                                              toSave <- paste0(toSave, paste0(.pointSet$pts$x[jn], ",", .pointSet$pts$y[jn], ",", as.character(pointSetList$pts$id[jn])))
                                         }
                                         temp1$Value[temp1$Metadata == 'polygonPts' & temp1$Time == .pointSet$frame] <- toSave
                                    }
                                    temp1$Time <- as.factor(temp1$Time);
                                    write.arff(temp1, file=file, relation='ROI-PointSet')
                               }
                          )
)
