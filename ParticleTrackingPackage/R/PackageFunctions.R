#' @import foreign
#' @import plyr
NULL

#' @title Take an arff file and reorganize it into a more standard 'table' format.
#' @description Specifically this is used to import an arff file from JEX as JEX
#' uses a column called 'Measurement' to define the type of measurment
#' or property being stored and 'Value', the value of that property.
#'
#' @param data An object that is the result of using foreign::read.arff(file) on an arff file
#' @param baseName An optional basename to add to whatever label is in the \code{nameCol} portion of each row entry
#' @param convertToNumeric An option to convert the columns of information within \code{data} leading up to
#' \code{nameCol} and \code{valueCol} to numeric or to leave as text. Default is to convert to numeric (i.e., TRUE)
#' @param nameCol The name of the column that describes the nature of the value in the \code{valueCol}
#' @param valueCol The name of the column with the values of the properties listed in the \code{nameCol}
reorganizeTable <- function(data, baseName=NA, convertToNumeric=TRUE, nameCol='Measurement', valueCol='Value')
{
     #require(plyr)
     idCols <- names(data)
     idCols <- idCols[-which(idCols %in% c(nameCol,valueCol))]
     newData <- data.frame(stringsAsFactors=FALSE)
     measurements <- unique(data[,nameCol])

     for(m in measurements)
     {
          if(is.na(baseName))
          {
               newColName <- m
               newColName <- gsub(' ','.',newColName, fixed=TRUE) # Get rid of extraneous spaces
          }else
          {
               newColName <- paste(baseName,'.',m, sep='')
               newColName <- gsub(' ','.',newColName, fixed=TRUE) # Get rid of extraneous spaces
          }

          temp <- data[data[,nameCol]==m,]
          temp2 <- temp[,c(idCols,valueCol)]
          if (length(idCols) == 0) {
               temp2 <- data.frame(ReallyRandomNameYo = temp2)
               names(temp2) <- newColName
          }
          else {
               names(temp2)[names(temp2) == valueCol] <- newColName
          }
          if(nrow(newData) == 0)
          {
               newData <- temp2
          }else
          {
               newData <- merge(newData, temp2, by=idCols)
          }
     }

     if(convertToNumeric)
     {
          for(n in idCols)
          {
               newData[,n] <- as.numeric(as.character(newData[,n]))
          }
     }

     return(newData)
}

#' Return information and data for a log frequency sweep with the given parameters
#'
#' @param amplitude numeric
#' @param phaseShift numeric - default=0
#' @param offset numeric - defulat=0
#' @param sin boolean T=sinsoid and F=triangular - default=FALSE
#' @param fi numeric initial frequency of the sweep
#' @param ti numeric initial time at which the sweep starts
#' @param sweepDuration numeric duration of the sweep
#' @param t numeric time at which to calculate the function value
#' @param guess list of parameter values to override the inidividual parameter values
#' @param calcVelocity logical value whether to calculate the velocity or not
#' @param ff numeric final frequency of the sweep
getSweep <- function(amplitude=1, phaseShift=0, offset=0, sin=FALSE, ti=0, fi=2, ff=0.1, sweepDuration=300, t=seq(0,300,0.05), guess=NULL, calcVelocity=TRUE)
{
     tOriginal <- t
     t <- t-ti
     tf <- sweepDuration
     N <- log(ff/fi)/log(2)
     R <- N / (sweepDuration) # Usually tf-ti but time adjusting t by ti

     if(is.null(guess))
     {
          A <- amplitude
          phi <- phaseShift
          b <- offset
     }
     else
     {
          A <- guess['amplitude']
          phi <- guess['phaseShift']
          b <- guess['offset']
     }

     offsetInflections <- (  (-1*(phi/(pi/2))) %/% 1  )
     nf <- -((4*fi-4*ff)*pi*tf+2*log(ff/fi)*phi)/(log(ff/fi)*pi)
     nf <- nf + 1 # for good measure.
     ni <- seq(0,nf,1) + offsetInflections
     suppressWarnings(inflections <- (log((log(ff/fi)*phi)/(2*fi*pi*tf)+(log(ff/fi)*ni)/(4*fi*tf)+1)*tf)/log(ff/fi))
     inflections <- inflections[!is.nan(inflections)]
     inflectionNums <- (  seq(0,length(inflections)-1,1) + (  offsetInflections %% 4  )  ) %% 4
     startTimeI <- which(inflections >= t[1])[1]
     endTimeI <- which(inflections <= t[length(t)])
     endTimeI <- endTimeI[length(endTimeI)]
     if(is.na(endTimeI))
     {
          stop(paste(t[length(t)], "is greater than last time of last inflection point. Therefore, there must be an error in fi/ff or tAll as there are inflections guaranteed in during times ti to tf given fi and ff."))
     }
     inflections <- inflections[startTimeI:endTimeI]
     inflectionNums <- inflectionNums[startTimeI:endTimeI]

     v <- NULL
     if(sin)
     {
          predicted <- A*sin(2*pi*((fi*(-1+2^(R*t)))/(R*log(2))) - phi) + b
          if(calcVelocity)
          {
               v=getDerivative(x=predicted, t=t)
          }
          return(list(A=A, t=tOriginal, x=predicted, v=v, inflections=inflections+ti, inflectionNums=inflectionNums))
     }
     else
     {
          predicted <- (2*A/pi)*asin(sin(2*pi*((fi*(-1+2^(R*t)))/(R*log(2))) - phi)) + b
          if(calcVelocity)
          {
               v=getDerivative(x=predicted, t=t)
          }
          return(list(A=A, t=tOriginal, x=predicted, v=v, inflections=inflections+ti, inflectionNums=inflectionNums))
     }
}


#' @title Get a table of time vs frequency for a log descending frequency sweep
#'
#' @param t A numeric vector of times [s]
#' @param fi A numeric value of the initial frequency of the sweep [Hz]
#' @param ff A numeric value of the final frequency of the sweep [Hz]
#' @param duration A numeric value of the duration of the sweep [s]
#'
#' @export
getFrequencies <- function(t=seq(0,300,0.035), fi=1, ff=0.01, duration=300)
{
     return(data.frame(t=t, f=fi*(ff/fi)^(t/duration)))
}

#' @title Calculate shear stress for given parameters
#'
#' @param f A numeric value or vector of frequencies for which to calculate shear stress values
#' @param pixelAmplitude A numeric value of the maximum amplitude observed in an image [pixels] for a particl in the center streamline of a wide flat microchannel
#' @param h A numeric value of the microchannel height [meters]
#' @param mu A numeric value of the dynamic viscosity [kg / (m * s)]
#' @param cameraPixelSize A numeric value of the size of the pixels on the camera itself in [meters] (assumed square, default = 6.45e-6)
#' @param mag A numeric value of the maginification used on the microscope
#'
#' @export
getShearStress <- function(f, pixelAmplitude, h=200e-6, mu=0.00078, cameraPixelSize = 6.45e-6, mag=4)
{
     micronAmplitude <- (6.45e-6/mag)*pixelAmplitude
     uMax <- 4*micronAmplitude*f
     return((2/3)*uMax*(6*mu)/h)
}
