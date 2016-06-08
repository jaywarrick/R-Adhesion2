# Need to make a pointlist where each track is described by a "point"
#
# The point has coordinates of "x, y, t/frame".
#
# Again, movement in the x direction should be allowed while y is penalized. Also, the tl$
#
test1 <- function(tl, slots)
{
     before <- new('PointSet')
     before$initializeEmpty(frame=1)
     after <- new('PointSet')
     after$initializeEmpty(frame=2)
     slots <- slots[slots %in% names(tl$tracks[[1]]$pts)]
     for(track in tl$tracks)
     {
          do.call(before$addPoint, c(list(id=track$id), as.list(track$pts[1,slots]))) # first point
          do.call(after$addPoint, c(list(id=track$id), as.list(track$pts[track$frameCount(),slots]))) # last point
     }
     ret <- new('PointSetList')
     ret$setPointSet(before)
     ret$setPointSet(after)
     return(ret)
}

test2 <- function(track, fun, par0)
{
     getSweep(amplitude = 1, phaseShift = 0, offset = 0, sin = FALSE, ti = 0, fi = 2, ff = 0.1, sweepDuration = 300, t = seq(0, 300, 0.05), guess = NULL, calcVelocity = TRUE, flipped = TRUE)
}



