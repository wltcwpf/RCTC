#' A Function to Read COSMOS V2 (corrected) files
#'
#' This function reads SMC format corrected data and output a list of date time
#' network, station, channel, time step, and acceleration (in the unit of cm/sec2) of the record.
#' @param file_path A directory specifies the location of the SMC format files,
#' @param n_comp A number, can be 3 or 1. 3 indicates the SMC file contains all three components together,
#' while 1 indicates the SMC file only contains one component
#' @return A list is returned which includes the date time, network, station, channel, time step, and
#' acceleration time series (in the unit of cm/sec2) of the record
#' @keywords read_cosmosV2
#' @export

read_cosmosV2 <- function( file_path, n_comp = 3 ){

  if( n_comp == 3 ){

    head_time <- scan( file_path, sep ="", skip = 4, nlines = 1, what = character() )

    ymd <- gsub( ',', '', head_time[ 4 ] )

    ymd <- gsub( '-', '', as.Date( ymd, format = '%m/%d/%y' ) )

    hms <- strsplit( head_time[ 5 ], split = ':' )[[ 1 ]]

    hms <- paste0( c( hms[ 1:2 ], round( as.numeric( hms[ 3 ] ) ) ), collapse = '' )

    head_sensor <- scan( file_path, sep = "", nlines = 1, what = character() )

    ifelse( grepl( 'accel', head_sensor[ 2 ] ), sensor <- 'HN', sensor <- 'HH' )

    head_cha <- scan( file_path, sep = "", skip = 4, nlines = 1, what = character() )

    head_cha <- head_cha[ 1 ]

    head_sta <- scan( file_path, sep = "", skip = 5, nlines = 1, what = character() )

    head_sta <- head_sta[ 3 ]

    if( is.na( as.numeric( head_sta ) ) ){

      Idx <- gregexpr( head_sta, head_cha )

      netsta <- paste0( c( substr( head_cha, Idx[[ 1 ]][ 1 ] - 3, Idx[[ 1 ]][ 1 ] - 2 ),
                           substr( head_cha, Idx[[ 1 ]][ 1 ], Idx[[ 1 ]][ 1 ] + nchar( head_sta ) - 1 ) ),
                        collapse = '_' )
    }else{

      netsta <- paste0( 'CE', '_', head_sta )

    }

    nsch <- c( paste0( c( netsta, paste0( sensor, 'E', sep = '' ) ), collapse = '_' ),
               paste0( c( netsta, paste0( sensor, 'N', sep = '' ) ), collapse = '_' ),
               paste0( c( netsta, paste0( sensor, 'Z', sep = '' ) ), collapse = '_' ) )

    head_dt1 <- scan( file_path, sep = "", skip = 45, nlines = 1, what = character() )

    npts1 <- as.numeric( head_dt1[ 1 ] )

    num_ind <- which( !is.na( as.numeric( head_dt1 ) ) )

    dt1 <- as.numeric( head_dt1[ num_ind[ 2 ] ] )

    if( dt1 > 1 ) dt1 <- 1 / dt1

    len_h1 <- 45 + ( ceiling( npts1 / 8 ) + 1 ) * 3 + 1

    head_dt2 <- scan( file_path, sep = "", skip = len_h1 + 45, nlines = 1, what = character() )

    npts2 <- as.numeric( head_dt2[ 1 ] )

    num_ind <- which( !is.na( as.numeric( head_dt2 ) ) )

    dt2 <- as.numeric( head_dt2[ num_ind[ 2 ] ] )

    if( dt2 > 1 ) dt2 <- 1 / dt2

    len_h2 <- 45 + ( ceiling( npts2 / 8 ) + 1 ) * 3 + 1

    head_dt3 <- scan( file_path, sep = "", skip = len_h1 + len_h2 + 45, nlines = 1, what = character() )

    npts3 <- as.numeric( head_dt3[ 1 ] )

    num_ind <- which( !is.na( as.numeric( head_dt3 ) ) )

    dt3 <- as.numeric( head_dt3[ num_ind[ 2 ] ] )

    if( dt3 > 1 ) dt3 <- 1 / dt3

    if( ( dt1 != dt2 ) | ( dt1 != dt3 ) | ( dt2 != dt3 ) ){

      print( 'Time steps between three components are not the same!!' )

      stop()
    }

    line1 <- readLines( file_path, n = ceiling( npts1 / 8 ) + 46 )

    data1 <- c()

    for( k in seq( 47, ceiling( npts1 / 8 ) + 46 ) ){

      if( length( data1 ) == 0 ){

        data1 <- as.numeric( substring( line1[ k ], c( 1, 11, 21, 31, 41, 51, 61, 71 ),
                                        c( 10, 20, 30, 40, 50, 60, 70, 80 ) ) )

      } else{

        data1 <- c( data1, as.numeric( substring( line1[ k ], c( 1, 11, 21, 31, 41, 51, 61, 71 ),
                                                  c( 10, 20, 30, 40, 50, 60, 70, 80 ) ) ) )
      }
    }

    data1 <- data1[ 1:npts1 ]

    line2 <- readLines( file_path, n = len_h1 + ceiling( npts2 / 8 ) + 46 )

    data2 <- c()

    for( k in seq( len_h1 + 47, len_h1 + ceiling( npts2 / 8 ) + 46 ) ){

      if( length( data2 ) == 0 ){

        data2 <- as.numeric( substring( line2[ k ],
                                        c( 1, 11, 21, 31, 41, 51, 61, 71 ),
                                        c( 10, 20, 30, 40, 50, 60, 70, 80 ) ) )
      } else{

        data2 <- c( data2, as.numeric( substring( line2[ k ],
                                                  c( 1, 11, 21, 31, 41, 51, 61, 71 ),
                                                  c( 10, 20, 30, 40, 50, 60, 70, 80 ) ) ) )
      }
    }

    data2 <- data2[ 1:npts2 ]

    line3 <- readLines( file_path, n = len_h1 + len_h2 + ceiling( npts3 / 8 ) + 46 )

    data3 <- c()

    for( k in seq( len_h1 + len_h2 + 47, len_h1 + len_h2 + ceiling( npts3 / 8 ) + 46 ) ){

      if( length( data3 ) == 0 ){

        data3 <- as.numeric( substring( line3[ k ],
                                        c( 1, 11, 21, 31, 41, 51, 61, 71 ),
                                        c( 10, 20, 30, 40, 50, 60, 70, 80 ) ) )
      } else{

        data3 <- c( data3, as.numeric( substring( line3[ k ],
                                                  c( 1, 11, 21, 31, 41, 51, 61, 71 ),
                                                  c( 10, 20, 30, 40, 50, 60, 70, 80 ) ) ) )
      }
    }

    data3 <- data3[ 1:npts3 ]

    res <- list()

    res$datetime <- paste0( ymd, hms, collapse = '' )

    res$netstacha <- nsch

    res$dt <- dt1

    # NOTE: unit is in cm/sec2
    res$data1 <- data1

    res$data2 <- data2

    res$data3 <- data3

    return( res )

  }else if( n_comp == 1 ){

    n_comment1 <- scan( file_path, sep = "", skip = 13, nlines = 1, what = character() )

    n_comment1 <- as.numeric( n_comment1[ 6 ] )

    n_comment2 <- scan( file_path, sep = "", skip = 14 + n_comment1, nlines = 1, what = character() )

    n_comment2 <- as.numeric( n_comment2[ 6 ] )

    n_comment3 <- scan( file_path, sep = "", skip = 14 + n_comment1 + 1 + n_comment2, nlines = 1, what = character() )

    n_comment3 <- as.numeric( n_comment2[ 1 ] )

    head_time <- scan( file_path, sep = "", skip = 7, nlines = 1, what = character() )

    ymd <- gsub( '/', '', head_time[ 4 ] )

    hms <- strsplit( head_time[ 5 ], split = ':' )[[ 1 ]]

    hms <- paste0( c( hms[ 1:2 ], round( as.numeric( hms[ 3 ] ) ) ), collapse = '' )

    head_netstacha <- scan( file_path, sep = "", skip = 14 + n_comment1 + 1 + n_comment2 + 3, nlines = 1, what = character() )

    nsch <- gsub( '[|]<SCNL>', '', head_netstacha[ 1 ] )

    nsch <- paste0( strsplit( nsch, split = '[.]' )[[ 1 ]][ c( 3, 1, 2 ) ], collapse = '_' )

    head_pt <- scan( file_path, sep = "", skip = 14 + n_comment1 + 1 + n_comment2 + 1 + n_comment3, nlines = 1, what = character() )

    npts <- as.numeric( head_pt[ 1 ] )

    head_dt <- scan( file_path, sep = "", skip = 14 + n_comment1 + 1 + n_comment2 + 1 + 4, nlines = 1, what = character() )

    num_ind <- which( !is.na( as.numeric( head_dt ) ) )

    dt <- as.numeric( head_dt[ num_ind[ 1 ] ] )

    if( dt > 1 ) dt <- 1 / dt

    data <- scan( file_path, sep = "", skip = 14 + n_comment1 + 1 + n_comment2 + 1 + n_comment3 + 1, nlines = npts, what = numeric() )

    res <- list()

    res$datetime <- paste0( ymd, hms, collapse = '' )

    res$netstacha <- nsch

    res$dt <- dt

    # NOTE: unit is in cm/sec2
    res$data <- data

    return( res )
  }
}
