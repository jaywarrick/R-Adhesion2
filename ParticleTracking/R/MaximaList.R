#' @import methods
NULL

#' A class for storing and working with list of Maxima objects. This class is the
#' class that actually does the point tracking over time.
#' @field maxima A list of Maxima objects
#' @field maxID The maximum ID encountered in all the Maxima in the MaximaList
#' @field trackBackStart A single numeric value indicating the frame at which cells were tracked backward from (i.e., usually a late frame)
#' @field trackBackEnd A single numeric value indicating the frame at which cells were tracked back to (i.e., usually an early frame)
#' @field frameDimName A string (character vector) value indicating the frame at which cells were tracked back to (i.e., usually an early frame)
MaximaList <- setRefClass('MaximaList',
                          fields = list(maxima='list', maxID='numeric', trackBackStart='numeric', trackBackEnd='numeric', frameDimName='character'),
                          methods = list(
                               initializeWithJEXROIFile = function(path, frameDimName='Time')
                               {
                                    "Initialize with a file (eg., .jxd ROI file from JEX) that stores all the point sets as semicolon separated comma separated point information (see Maxima documentation)."

                                    frameDimName <<- frameDimName
                                    maximaFile <- read.arff(path)
                                    initializeWithROIDataFrame(roiTable=maximaFile, frameDimName=frameDimName)
                               },
                               initializeWithROIDataFrame = function(roiTable, frameDimName='Time')
                               {
                                    "Initialize with a data.frame of an ROI (eg., .jxd ROI file from JEX) that stores all the point sets as semicolon separated comma separated point information (see Maxima documentation)."

                                    # require(foreign)
                                    if(!is.null(roiTable))
                                    {
                                         maximaFile2 <- reorganize(roiTable, measurementCols='Metadata')
                                         maximaFile2[,frameDimName] <- as.numeric(as.character(maximaFile2[,frameDimName]))
                                         temp <- list()
                                         for(r in 1:nrow(maximaFile2))
                                         {
                                              newMaxima <- new('Maxima')
                                              newMaxima$initializeWithROI(frame=maximaFile2[r,frameDimName], polygon=maximaFile2[r,'polygonPts'])
                                              temp[[as.character(maximaFile2[r,frameDimName])]] <- newMaxima
                                         }
                                         temp <- temp[order(as.numeric(names(temp)), decreasing=FALSE)]
                                         maxima <<- temp
                                    }
                                    else
                                    {
                                         maxima <<- list()
                                    }

                                    maxID <<- -1 # Used for appending maxima
                               },
                               setFrame = function(maximaObject, frameNumber)
                               {
                                    "Set the maxima for a given frame number (as opposed to index)"

                                    maxima[[as.character(frameNumber)]] <<- newMaxima
                               },
                               setFrameAtIndex = function(maximaObject, index)
                               {
                                    "Set the maxima for a given index (as opposed to frame number)"

                                    maxima[[index]] <<- newMaxima
                               },
                               getIndexOfFrame = function(frame)
                               {
                                    "Get the index of a frame with a particular number"

                                    return(which(names(maxima)==as.character(frame)))
                               },
                               getMaximaAtIndex = function(index)
                               {
                                    "Get the maxima at a particular index (vs frame number)"

                                    return(maxima[[index]])
                               },
                               getMaxima = function(frame)
                               {
                                    "Get the maxima at a particular frame number (vs index)"

                                    return(maxima[[as.character(frame)]])
                               },
                               getPreviousMaxima = function(frame)
                               {
                                    "Get the maxima that preceeds the given frame (good for backtracking)"

                                    toGet <- getIndexOfFrame(frame)
                                    if(is.null(toGet) || toGet <= 1)
                                    {
                                         return(NULL)
                                    }
                                    return(maxima[[toGet-1]])
                               },
                               setMaxima = function(newMaxima)
                               {
                                    "Put a maxima into the MaximaList object"

                                    maxima[[as.character(newMaxima$frame)]] <<- newMaxima
                               },
                               length = function()
                               {
                                    "Return the number of sets of maxima in this MaximaList object"

                                    return(base::length(maxima))
                               },
                               trackBack = function(startFrame, endFrame, maxDist=150, direction=c(1,0,0), directionality=10, uniformityDistThresh=3, digits=1)
                               {
                                    "Track backwards from startFrame (usually from a later point in time of the acquistion) to endFrame (earlier in the acquisition)"

                                    if(is.null(getIndexOfFrame(startFrame)) || is.null(getIndexOfFrame(endFrame)))
                                    {
                                         stop("Start frame or end frame is not a valid frame in the maxima list.")
                                    }
                                    if(getIndexOfFrame(startFrame) > length() | getIndexOfFrame(startFrame) < 0 | getIndexOfFrame(endFrame) > length() | endFrame >= startFrame | getIndexOfFrame(endFrame) < 0)
                                    {
                                         stop("Start frame or end frame out of bounds. Start frame must be >=0 and > (yes >) end frame because we work backward. End frame must be < last frame of data set.")
                                    }
                                    trackBackStart <<- startFrame
                                    trackBackEnd <<- endFrame
                                    currentMaxima <- getMaxima(startFrame)
                                    previousMaxima <- getPreviousMaxima(startFrame)
                                    maxID <<- max(currentMaxima$points$id)
                                    while(!is.null(previousMaxima) && previousMaxima$frame >= endFrame)
                                    {
                                         cat("Linking frame: ", currentMaxima$frame, "\n", sep="")
                                         t1 <- currentMaxima$getXYZ()
                                         t0 <- previousMaxima$getXYZ() # maxima to prepend
                                         results <- DirectionalLinearAssignment(points=list(t0=t0, t1=t1), maxDist=maxDist, direction=direction, directionality=directionality, uniformityDistThresh=uniformityDistThresh, digits=digits)
                                         px <- results$px
                                         newMaxima <- new('Maxima')
                                         newMaxima$frame <- previousMaxima$frame
                                         #### Can potentially vectorize this process to make it faster
                                         for(i in 1:base::length(px))
                                         {
                                              if(px[i] <= nrow(t0) & i <= nrow(t1))
                                              {
                                                   #print("A")
                                                   # Then we have a link between t0 and t1 and the points should share ids
                                                   p0 <- previousMaxima$getPoint(index=px[i])
                                                   p1 <- currentMaxima$getPoint(index=i)
                                                   newMaxima$addPoint(id=p1$id, x=p0$x, y=p0$y)
                                              }
                                              else if(px[i] <= nrow(t0) & i > nrow(t1))
                                              {
                                                   #print("B")
                                                   # This point in t0 represents the start of a new track
                                                   # increment the maxID and add the point to the newMaxima
                                                   maxID <<- maxID + 1
                                                   p0 <- previousMaxima$getPoint(index=px[i])
                                                   newMaxima$addPoint(id=maxID, x=p0$x, y=p0$y)
                                              }
                                              else
                                              {
                                                   #print("C")
                                                   # This point in t1 represents the end of a track or an auxillary match to enable other types of matches
                                                   # We don't need to do anything
                                              }
                                         }
                                         setMaxima(newMaxima)
                                         currentMaxima <- newMaxima
                                         previousMaxima <- getPreviousMaxima(newMaxima$frame)
                                    }
                                    cat("Done tracking. Trimming Data to start and end frame")
                                    allFrameNames <- as.numeric(names(maxima))
                                    indicesToRemove <- which(!(allFrameNames %in% endFrame:startFrame))
                                    maxima[indicesToRemove] <<- NULL
                               },
                               getTrackList = function(t0_Frame, timePerFrame, frameRange=c(0,-1))
                               {
                                    "Create a TrackList object from this MaximaList object (essentially a list indexed/grouped by cell id instead of time)"

                                    if(frameRange[1] < 0 | frameRange[2] > trackBackStart | frameRange[2] < frameRange[1])
                                    {
                                         exportFrames <- names(maxima)[names(maxima) %in% as.character(trackBackEnd:trackBackStart)]
                                    }
                                    else
                                    {
                                         exportFrames <- names(maxima)[names(maxima) %in% as.character(frameRange[1]:frameRange[2])]
                                    }
                                    trackList <- new('TrackList')
                                    trackList$setStandardMeta(t0_Frame=t0_Frame, timePerFrame=timePerFrame)
                                    myCount = 0
                                    total = base::length(exportFrames)
                                    for(f in exportFrames)
                                    {
                                         .maxima <- getMaxima(f)
                                         for(i in 1:.maxima$length())
                                         {
                                              data <- .maxima$points[i,]
                                              trackList$addTrackPoint(id=data$id, x=data$x, y=data$y, frame=.maxima$frame)
                                         }
                                         myCount <- myCount + 1
                                         cat("Generating TrackList: ", round(100*myCount/total, digits=2), "%\n", sep="")
                                    }
                                    trackList$setStandardMeta(t0_Frame=t0_Frame, timePerFrame=timePerFrame) # Set again to ensure information filters down to all tracks
                                    return(trackList)
                               },
                               getProp = function(fun=function(.maxima){return(base::range(.maxima$points$x))})
                               {
                                    "Get the property of a maxima by applying the user-specified function on each maxima object"

                                    ret <- list()
                                    for(.maxima in maxima)
                                    {
                                         ret[[as.character(.maxima$frame)]] <- fun(.maxima)
                                    }
                                    return(ret)
                               },
                               generateMaximaPlots = function(path=NULL, baseName='Frame_')
                               {
                                    "Produce plots of the maxima locations over time and output them to a folder"

                                    dir.create(path, recursive=TRUE)
                                    wd <- getwd()
                                    setwd(path)
                                    xlim <- base::range(unlist(getProp(fun=function(.maxima){return(base::range(.maxima$points$x))})))
                                    ylim <- base::range(unlist(getProp(fun=function(.maxima){return(base::range(.maxima$points$y))})))

                                    makeThePlots <- function()
                                    {
                                         for(.maxima in maxima)
                                         {
                                              cat("Plotting Frame = ", .maxima$frame, "\n", sep="")
                                              pdf(file=paste(path, '/', baseName, .maxima$frame, '.pdf', sep=''), width=9, height=4)
                                              .maxima$plotMaxima(IDs=TRUE, xlim=xlim, ylim=ylim, main=as.character(.maxima$frame))
                                              dev.off()
                                         }
                                    }

                                    tryCatch(expr=makeThePlots(), finally=setwd(wd))
                               },
                               offsetFrames = function(offset=0)
                               {
                                    "Add the offset to all the frame numbers, adjusting trackBaskStart and trackBackEnd accordingly"

                                    names(maxima) <<- as.character(as.numeric(names(maxima))+offset)
                                    for(.maxima in maxima)
                                    {
                                         .maxima$frame <-.maxima$frame+offset
                                    }
                                    trackBackStart <<- trackBackStart + offset
                                    trackBackEnd <<- trackBackEnd + offset
                               },
                               saveROI = function(file)
                               {
                                    "Export this MaximaList object as an ROI object for JEX (https://github.com/jaywarrick/JEX)"
                                    print(as.factor(as.numeric(names(maxima))))
                                    temp1 <- data.frame(Time=as.numeric(names(maxima)), Metadata=rep(as.factor(c('Type','patternPts','polygonPts')), each=base::length(names(maxima))))
                                    temp1$Value[temp1$Metadata == 'Type'] <- '5'
                                    temp1$Value[temp1$Metadata == 'patternPts'] <- '0,0,0'
                                    for(.maxima in maxima)
                                    {
                                         toSave <- ''
                                         jn <- nrow(.maxima$points)
                                         if(jn > 1)
                                         {
                                              for(j in 1:(jn-1))
                                              {
                                                   toSave <- paste0(toSave, paste0(.maxima$points$x[j], ",", .maxima$points$y[j], ",", as.character(.maxima$points$id[j]), ";"))
                                              }
                                         }
                                         if(jn > 0)
                                         {
                                              # use of "as.character" on id is important to get the right cell label because the factor levels are per maxima table instead of spanning maxima tables
                                              toSave <- paste0(toSave, paste0(.maxima$points$x[jn], ",", .maxima$points$y[jn], ",", as.character(maxima$points$id[jn])))
                                         }
                                         temp1$Value[temp1$Metadata == 'polygonPts' & temp1$Time == .maxima$frame] <- toSave
                                    }
                                    temp1$Time <- as.factor(temp1$Time);
                                    write.arff(temp1, file=file, relation='ROI-Maxima')
                               }
                          )
)
