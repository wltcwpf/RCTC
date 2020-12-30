#' A Function to Read PEER (corrected) files
#'
#' This function reads PEER format corrected data and output a list of date, station code, channel code,
#' time step, and acceleration of the record.
#' @param file_path A directory specifies the location of the PEER format files,
#' @return A list is returned which includes the date, station code, channel code,
#' time step, and acceleration time series of the record
#' @keywords read_peer
#' @export

read_peer <- function( file_path ){

  head_meta <- scan( file_path, sep = ',', skip = 1, nlines = 1, what = character() )

  Date <- as.Date( head_meta[ 2 ], format = '%m/%d/%Y' )

  StationCode <- gsub( ' ', '', head_meta[ 3 ] )

  ChannelCode <- gsub( ' ', '', head_meta[ 4 ] )

  head_dt <- scan( file_path, sep = ",", skip = 3, nlines = 1, what = character() )

  dt <- as.numeric( gsub( "[^0-9\\.]", "", head_dt[ 2 ] ) )

  data <- scan( file_path, sep = "", skip = 4, what = numeric() )

  res <- list()

  res$Date <- Date

  res$StationCode <- StationCode

  res$ChannelCode <- ChannelCode

  res$dt <- dt

  res$data <- data

  return( res )

}
