#' @import foreign
#' @import data.table
#' @importFrom pracma numel
NULL

#' Reorganize a table from long form to wide form.
#'
#' This function can be applied to data.frame objects or data.table objects. It
#' returns the same type as given.
#'
#' @param data A data.table or data.frame
#' @param idCols A character vector of id column names (e.g., 'Id')
#' @param measurementCols A character vector of column names that describe measurements (e.g., a column
#' called 'Measurement' with values such as 'Min', 'Max', 'Mean', etc.
#' @param valueCols A character vector of column names (typically one) that contains the numeric data
#' to reorganize into wide format (e.g., 'Value')
#'
#' @export
reorganize <- function(data, idCols=NULL, measurementCols='Measurement', valueCols='Value')
{
     isDataTable <- FALSE
     if(is.data.table(data))
     {
          isDataTable <- TRUE
     }
     else
     {
          data <- data.table(data)
     }

     # If idCols = NULL, then use all remaining cols except measurementCols and valueCols
     if(is.null(idCols))
     {
          idCols <- names(data)[!(names(data) %in% c(measurementCols, valueCols))]
     }

     formula <- as.formula(paste(paste(idCols, collapse='+'), " ~ ", paste(measurementCols, collapse='+')))
     print(formula)
     data <- dcast(data, as.formula(paste(paste(idCols, collapse='+'), " ~ ", paste(measurementCols, collapse='+'))), value.var = valueCols)
     if(isDataTable)
     {
          return(data)
     }
     else
     {
          return(data.frame(data))
     }
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
#' @param flipped logical indicating whether to "flip" the curve using a 180 degree phase shift
#'
#' @export
getSweep <- function(amplitude=1, phaseShift=0, offset=0, sin=FALSE, ti=0, fi=2, ff=0.1, sweepDuration=300, t=seq(0,300,0.05), guess=NULL, calcVelocity=TRUE, flipped=TRUE)
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
          if(flipped)
          {
               phi <- phi + pi # Flipped to match arduino actuation
          }
     }
     else
     {
          A <- guess['amplitude']
          phi <- guess['phaseShift']
          b <- guess['offset']
          if(flipped)
          {
               phi <- phi + pi # Flipped to match arduino actuation
          }
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
          stop(paste(t[length(t)], "is greater than last time of last inflection point. Therefore, there must be an error in fi/ff or allTimes as there are inflections guaranteed in during times ti to tf given fi and ff."))
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


#' Get a table of time vs frequency
#'
#' Return the t vs f table for a descending log-frequency sweep
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

#' Calculate shear stress for given parameters
#'
#' Based on laminar flow in a microchannel where the height is << width.
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

#' lseq
#'
#' Generate a sequence on a log scale
#'
#' @param from numeric starting value
#' @param to numeric ending value
#' @param length.out numeric number of numbers in the sequence between from and to
#'
#' @export
lseq <- function (from, to, length.out)
{
     exp(seq(log(from), log(to), length.out = length.out))
}
