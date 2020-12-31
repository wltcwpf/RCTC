#' A Sinc interpolation Function
#'
#' This function applies Sinc interpolation to reduce the time step and
#' returns a time series with more data points than in \code{data}
#' @param data An acceleration time series array
#' @param interpolation_factor An interpolated factor, specifies the degree of interpolation
#' and should be a power of two.
#' @return The interpolated time series of \code{data} by factor \code{interpolation_factor}
#' @keywords internal
#' @export
#' @examples
#' Interpft(c(0.01,0.03,...0.5), 64)

########################################
######## inputs: time series, interpolation multiplier
######## outputs: Sinc interpolated time series
Interpft <- function(data, interpolation_factor){
  ## check if required packages are installed
  list.of.packages <- c("stats", "pracma")
  new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  if(length(new.packages)) install.packages(new.packages)
  ## check if libraries are loaded
  if(!('package:stats' %in% search())){library(stats)}
  if(!('package:pracma' %in% search())){library(pracma)}
  npts <- length(data)
  npw2_exp = floor(log(npts)/log(2))
  npw2_temp = 2^npw2_exp
  if(npw2_temp == npts){
    npw2 = npts
  } else{
    npw2 <- 2^(npw2_exp+1)
    data <- c(data, zeros((npw2-npts),1))
  }
  ## ensure length(data) is a power of 2
  nyqst = npw2/2 + 1
  m = length(data)
  a = fft(data)
  b = c(a[1:nyqst], zeros((m*interpolation_factor-m),1), a[(nyqst+1):m])
  b[nyqst] = b[nyqst]/2
  b[nyqst+m*interpolation_factor-m] = b[nyqst]
  y = ifft(b)
  y = y*interpolation_factor
  y = Re(y)
  y <- y[1:(npts*interpolation_factor)]
  return(y)
}




