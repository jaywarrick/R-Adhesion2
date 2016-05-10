#' Filter for testing tracks based on a min and max length
#'
#' @param track Track The track to filter
#' @param min numeric The minimum length of a track (inclusive, >=)
#' @param max numeric The maximum length of a track (inclusive, <=)
#'
#' @return boolean whether the track meets the specified conditions
#'
#' @export
trackLengthFilter <- function(track, min=0, max=1000000)
{
     l <- track$frameCount()
     if(l >= min & l <=max)
     {
          return(TRUE)
     }
     else
     {
          return(FALSE)
     }
}

#' trackMinRangeFilter
#'
#' Filter for testing tracks based on a min displacement [pixels]
#'
#' @param track Track The track to filter
#' @param slot character value indicating the slot/column in the pts data.frame to analyze
#' @param min numeric value indicating the minimum range (max-min) of the values in the specified slot allowed (inclusive, <=)
#'
#' @return boolean whether the track meets the specified conditions
#'
#' @export
trackMinRangeFilter <- function(track, slot, min=4)
{
     l <- track$range(slot=slot, rel=FALSE)
     l <- l[2]-l[1]
     if(l >= min)
     {
          return(TRUE)
     }
     else
     {
          return(FALSE)
     }
}

#' trackMaxRangeFilter
#'
#' Filter for testing tracks based on a min displacement [pixels]
#'
#' @param track Track The track to filter
#' @param slot character value indicating the slot/column in the pts data.frame to analyze
#' @param max numeric value indicating the maximum range (max-min) of the values in the specified slot allowed (inclusive, <=)
#'
#' @return boolean whether the track meets the specified conditions
#'
#' @export
trackMaxRangeFilter <- function(track, slot, max=4)
{
     l <- range(track$pts[,slot])
     l <- l[2]-l[1]
     if(l <= max)
     {
          return(TRUE)
     }
     else
     {
          return(FALSE)
     }
}

#' trackFrameFilter
#'
#' Filter for testing if tracks exist within indicated frame bounds.
#'
#' @param track Track The track to be tested
#' @param startMin numeric The min (inclusive) frame for the starting frame of the track - default=0
#' @param startMax numeric The max (inclusive) frame for the starting frame of the track - default=1000000
#' @param endMin numeric The min (inclusive) frame for the ending frame of the track - default=0
#' @param endMax numeric The max (inclusive) frame for the ending frame of the track - default=1000000
#' @return boolean whether the track meets the specified conditions
#'
#' @export
trackFrameFilter <- function(track, startMin=0, startMax=1000000, endMin=0, endMax=1000000)
{
     startFrame <- track$pts$frame[1]
     endFrame <- last(track$pts$frame)

     if(startFrame >= startMin & startFrame <= startMax & endFrame >= endMin & endFrame <= endMax)
     {
          return(TRUE)
     }
     else
     {
          return(FALSE)
     }
}
