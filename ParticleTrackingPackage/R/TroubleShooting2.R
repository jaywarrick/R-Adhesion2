bf <- butter(3, 0.1)                        # 10 Hz low-pass filter
t <- trackList$meta$tAll
x <- colMeans(trackList$getMatrix(), na.rm=T)
y <- filtfilt(bf, x)
lines(t, y, col="blue", lwd=1)
#points(t, z, col="blue")
legend("bottomleft", legend = c("data", "filtfilt", "filter"),
       pch = 1, col = c("black", "red", "blue"), bty = "n")

t(y)
