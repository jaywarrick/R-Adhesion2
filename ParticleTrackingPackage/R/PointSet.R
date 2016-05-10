#' @import methods
NULL

#' PointSet
#'
#' Class for storing and working with a set of points. The points are stored as a data.frame field called 'pts'.
#' The first column of the pionts data.frame is the 'id' column representing the id for each point
#' within the PointSet. The remaining columns represent the 'position' of the point in the n-dimensional space.
#' Technically, at this point, each dimension does not have to necessarily be numeric in nature.
#'
#' @field frame A length one numeric vector
#' @field pts A data.frame storing the id, x, and y information of each point in the set of pts
PointSet <- setRefClass("PointSet",
                      fields = list(frame="numeric", pts="data.frame"),
                      methods = list(
                           initializeEmpty = function(frame)
                           {
                                "Create an empty PointSet with a NULL pts field. Every PointSet must be given a frame number (an assumption needed by the PointSetList class)."

                                frame <<- frame
                           },
                           initializeWithDataFrame = function(frame, data)
                           {
                                "Initialize this PointSet with the given frame number and a pts data.frame where each row represents an n-dimensional point. If the pts object is not a data.frame, an attempt will be made to convert to a data.frame using as.data.frame (data.table objects are also converted). If there is no column labeled 'id', one is prepended in the first column position and is equivalent to row number in the pts data.frame. All other columns are treated as positions in the arbitrary n-dimensional point-space."

                                frame <<- frame

                                temp <- data
                                if(!is.data.frame(pts))
                                {
                                     temp <- as.data.frame(pts)
                                }
                                if('id' %in% names(temp))
                                {
                                     # ensure the id column as the first column
                                     temp <- cbind(temp$id, temp[,!(names(temp) %in% 'id')])
                                }
                                else
                                {
                                     # make an id column
                                     temp <- cbind(data.frame(id=1:nrow(temp)), temp)
                                }

                                pts <<- temp
                           },
                           initializeWithROI = function(frame, ptString=NULL)
                           {
                                "initialize this with information from a semicolon separated list of comma separated point information (i.e., 'x,y,id;x,y,id;x,y,id;...')"

                                frame <<- as.numeric(as.character(frame))
                                if(!is.null(ptString))
                                {
                                     pairs <- strsplit(ptString,';')[[1]]
                                     x <- numeric(0)
                                     x0 <- numeric(0)
                                     y <- numeric(0)
                                     y0 <- numeric(0)
                                     index <- numeric(0)
                                     first <- TRUE
                                     for(pair in pairs)
                                     {
                                          nums <- strsplit(pair,',')[[1]]
                                          x <- append(x, as.numeric(nums[1]))
                                          y <- append(y, as.numeric(nums[2]))
                                          index <- append(index,as.numeric(nums[3]))
                                     }

                                     pts <<- data.frame(id=index, x=x, y=y)
                                }
                                else
                                {
                                     initialize(frame=frame)
                                }
                           },
                           pointCount = function()
                           {
                                "Return the number of pts in this PointSet object"

                                return(nrow(pts))
                           },
                           addPoint = function(id, ...)
                           {
                                "Add a point to this PointSet object. The list of supplied objects must represent a row of the pts data.frame with names that match the non-id columns of the pts data.frame exactly. If the pts, data.frame is initially empty, it will be created."

                                if(nrow(pts) > 0)
                                {
                                     if((ncol(pts) - 1) > base::length(...))
                                     {
                                          stop(paste0("The number of arguments does not match the number needed to properly define a point based upon the data already existing in the PoinSet.\n\nCurrent column names (including 'id'): ", paste(names(pts)[!names(pts) %in% 'id'])))
                                     }
                                }

                                if(id %in% pts$id)
                                {
                                     stop(paste0("Cannot add the point. The PointSet already contains a point with this id: ", paste(names(pts), collapse=','), ' -> ', paste(getPoint(id=id), collapse=',')))
                                }

                                # Note that rbind behaves well if pts is initially NULL
                                pts <<- rbind(pts, data.frame(id=id, list(...)))
                           },
                           plotPointSet = function(IDs=TRUE, pch=20, x='x', y='y', label='id', xlab='X [pixel]', ylab='Y [pixel]', cex=1, add=FALSE, ...)
                           {
                                "By default, it attempts to plot the x, y locations of the pts in this PointSet object with point ids as labels.
                                However, you can specify which two attributes of the pts to plot by specifying the attribute names for x and y.
                                Likewise, you can override the xlab, ylab, and cex args as well as a variabel arugment list which is passed to the
                                plot command. The label argument defines the column used for labeling each point. Thus, you could label each point
                                with the 'id' of the point (default behavior) or choose a different attribute by name."

                                if(add)
                                {
                                     points(pts[,x], pts[,y], pch=20, cex=cex, ...)
                                }
                                else
                                {
                                     plot(pts[,x], pts[,y], pch=20, xlab=xlab, ylab=ylab, cex=cex, ...)
                                }
                                if(IDs)
                                {
                                     text(x=pts[,x], y=pts[,y], labels=as.character(pts[,label]), cex=0.7, adj=c(0,0))
                                }
                           },
                           getPointsWithoutId = function()
                           {
                                "Convenience function to return the positions of the pts as a dataframe without the id column."

                                return(pts[,!(names(pts) %in% 'id')])
                           },
                           getPoint = function(index=NULL, id=NULL)
                           {
                                "Get a point at a specific index or with a specific id value as a data.frame.
                                Only one can be specified. Otherwise NULL is returned."

                                if(is.null(index) & is.null(id))
                                {
                                     stop("Can't return a point with parameters index and id as NULL.")
                                }
                                else if(!is.null(index))
                                {
                                     return(pts[index,])
                                }
                                else
                                {
                                     return(pts[pts$id==id,])
                                }
                           }
                      )
)
