#' Filter for testing tracks based on a min and max length
#' @param track Track The track to filter
#' @param min numeric The minimum length of a track (inclusive, >=)
#' @param max numeric The maximum length of a track (inclusive, <=)
#' @return boolean whether the track meets the specified conditions
trackLengthFilter <- function(track, min=0, max=1000000)
{
     l <- track$length()
     if(l >= min & l <=max)
     {
          return(TRUE)
     }
     else
     {
          return(FALSE)
     }
}

#' Filter for testing if tracks exist within indicated frame bounds.
#' @param track Track The track to be tested
#' @param startMin numeric The min (inclusive) frame for the starting frame of the track - default=0
#' @param startMax numeric The max (inclusive) frame for the starting frame of the track - default=1000000
#' @param endMin numeric The min (inclusive) frame for the ending frame of the track - default=0
#' @param endMax numeric The max (inclusive) frame for the ending frame of the track - default=1000000
#' @return boolean whether the track meets the specified conditions
trackFrameFilter <- function(track, startMin=0, startMax=1000000, endMin=0, endMax=1000000)
{
     startFrame <- track$points$frame[1]
     endFrame <- last(track$points$frame)

     if(startFrame >= startMin & startFrame <= startMax & endFrame >= endMin & endFrame <= endMax)
     {
          return(TRUE)
     }
     else
     {
          return(FALSE)
     }
}
