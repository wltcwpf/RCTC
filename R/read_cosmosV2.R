#' A Function to Read COSMOS V2 (corrected) files
#'
#' This function reads SMC format corrected data and output a list with time step and acceleration.
#' @param file_path A directory specifies the location of the SMC format files,
#' @param n_comp A number, can be 3 or 1. 3 indicates the SMC file contains all three components together,
#' while 1 indicates the SMC file only contains one component

read_cosmosV2 <- function( file_path, n_comp = 3 ){

  if( n_comp == 3 ){

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

    res$dt <- dt1

    # NOTE: unit is in cm/sec2
    res$data1 <- data1

    res$data2 <- data2

    res$data3 <- data3

    return( res )

  }else if( n_comp == 1 ){

    head_pt <- scan( file_path, sep = "", skip = 54, nlines = 1, what = character() )

    npts <- as.numeric( head_pt[ 1 ] )

    head_dt <- scan( file_path, sep = "", skip = 50, nlines = 1, what = character() )

    num_ind <- which( !is.na( as.numeric( head_dt ) ) )

    dt <- as.numeric( head_dt[ num_ind[ 1 ] ] )

    if( dt > 1 ) dt <- 1 / dt

    data <- scan( file_path, sep = "", skip = 55, nlines = npts, what = numeric() )

    res <- list()

    res$dt <- dt

    # NOTE: unit is in cm/sec2
    res$data <- data

    return( res )
  }
}
