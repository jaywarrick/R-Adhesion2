#' @title Track particles and calculate metrics of tracks.
#'
#' @details Information about particle locations and identities typically start out as 'points'.
#' These points are typically stored by frame (i.e., an image at a specific time). The points
#' in a specific frame are termed 'maxima'. Thus the base object of the package is called 'Maxima'.
#' A list of Maxima objects is called a 'MaximaList'. The MaximaList object holds the ability to
#' link points from one frame to another, also known as tracking. The package uses an algorithm
#' for the linear assigment problem (LAP) provided by the 'GraphAlignment' package in the
#' 'Bioconductor' package repository. See this paper (\url{http://www.nature.com/nmeth/journal/v5/n8/full/nmeth.1237.html}) for details.
#' Only point linking (not track linking) is implemented. However, this package provides additional
#' parameters for tracking particles amidst bulk motion of the group (see 'trackBack' method of
#' the 'MaximaList' class.
#'
#' Once tracked, the MaximaList object can be used to generate
#' a 'TrackList' object which reorganizes the data to be organized as a list of 'Track' objects.
#' Each Track object is just a table of the data for points belonging to the same track. The
#' TrackList object provides the ability to filter and calculate metrics of 'Track' objects easily.
#'
#' Additional functions are included in this package that are utilized for fitting tracks to
#' specific driving functions. These will likely be less useful to most people but hopefully
#' the base classes are useful to others doing particle tracking.
#'
#' @seealso \code{\link{PointSet}}
#' @seealso \code{\link{PointSetList}}
#' @seealso \code{\link{Track}}
#' @seealso \code{\link{TrackList}}
"_PACKAGE"
#> [1] "_PACKAGE"
