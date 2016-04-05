#' @import methods
NULL

#' Class for storing and working with a set of points
#' @field frame A length one numeric vector
#' @field points A data.frame storing the id, x, and y information of each point in the set of points
Maxima <- setRefClass("Maxima",
                      fields = list(frame="numeric", points="data.frame"),
                      methods = list(
                           initializeWithROI = function(frame, polygon=NULL)
                           {
                                "initialize this with information from a semicolon separated list of comma separated point information (i.e., 'x,y,id;x,y,id;x,y,id;...')"

                                frame <<- as.numeric(as.character(frame))
                                if(!is.null(polygon))
                                {
                                     pairs <- strsplit(polygon,';')[[1]]
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

                                     points <<- data.frame(id=index, x=x, y=y)
                                }
                                else
                                {
                                     points <<- data.frame(id=numeric(0), x=numeric(0), y=numeric(0))
                                }
                           },
                           length = function()
                           {
                                "Return the number of points in this Maxima object"

                                return(nrow(points))
                           },
                           addPoint = function(id, x, y)
                           {
                                "Add a point to this Maxima object"

                                points <<- rbind(points, data.frame(id=id, x=x, y=y))
                           },
                           plotMaxima = function(IDs=TRUE, pch=20, xlab='X [pixel]', ylab='Y [pixel]', cex=0.5, ...)
                           {
                                "Plot the x, y locations of the points in this Maxima object with point ids as labels"

                                plot(points$x, points$y, pch=20, xlab=xlab, ylab=ylab, cex=cex, ...)
                                if(IDs)
                                {
                                     text(x=points$x, y=points$y, labels=as.character(points$id), cex=0.5, adj=c(0,0))
                                }
                           },
                           getXYZ = function()
                           {
                                "Return the X, Y, and Z locations for all the points in this Maxima object. If Z is not given, returns 0 for Z."

                                data <- cbind(points, data.frame(z=0))
                                return(data[,2:4])
                           },
                           getPoint = function(index=NULL, id=NULL)
                           {
                                "Get a point at a specific index or with a specific id value.
                                Only one can be specified."

                                if(is.null(index) & is.null(id))
                                {
                                     return(NULL)
                                }
                                else if(!is.null(index))
                                {
                                     return(points[index,])
                                }
                                else
                                {
                                     return(points[points$id==id,])
                                }
                           }
                      )
)
