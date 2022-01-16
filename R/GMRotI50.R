#' A Function for GMRotI50
#'
#' This function computes Orientation-Independent Measures of Ground Motion -- GMRotI50
#' @param GMRotD Geometric Mean of rotated gorund motions. It is from the output list of \code{GM_RotD_cal}.
#' You can get it by: \code{result} = \code{GM_RotD_cal(...)}, \code{GMRotD} = \code{result[[2]]}.
#' It is a matrix with length of number of oscillator periods,
#' width of number of rotation angles (from 1 to 90 with 1 degree increment)
#'
#' @param period_t A array of oscillator periods, it is associated with periods of \code{GMRotD}
#' @param tmax The maximum usable period for penalty calculation, its default value is 10s
#' @param tmin The minimum usable period for penalty calculation, its defualt value is 0s
#' @return A list is returned with GMRotI50, GMRotI50_ang, penalty, the index of tmax in \code{period_t}
#' and the index of tmin in \code{period_t}
#' @keywords GMRotI50_cal
#' @importFrom stats median
#' @export

GMRotI50_cal <- function(GMRotD, period_t, tmax=10, tmin=0){
  # find the index for tmax and tmin
  if(tmax > max(period_t)){
    tmax_index <- length(period_t)
  }else{
    for(i in seq(1, length(period_t))){
      if(period_t[i] == tmax){
        tmax_index <- i
        break
      }
      if(period_t[i] > tmax){
        tmax_index <- i -1
        break
      }
    }
  }
  for(i in seq(length(period_t), 1)){
    if(period_t[i] == tmin){
      tmin_index <- i
      break
    }
    if(period_t[i] < tmin){
      tmin_index <- i + 1
      break
    }
    if(tmin == 0){
      tmin_index <- 1
      break
    }
  }

  # a vector that GM_PSA VS period
  GMRotI <- seq(1,length(period_t))
  penalty <- seq(1,90) * 0
  for (theta in seq(1,90)){
    GM_ratio <- 0
    for (per_index in seq(tmin_index, tmax_index)){
      GM_ratio <- GM_ratio+(GMRotD[theta, per_index]/stats::median(GMRotD[,per_index])-1)^2
    }
    penalty[theta] <- GM_ratio/(tmax_index-tmin_index+1)
  }
  GMRotI50_ang <- which.min(penalty)
  GMRotI <- GMRotD[GMRotI50_ang, ]

  res <- list()
  res$GMRotI <- GMRotI
  res$GMRotI50_ang <- GMRotI50_ang
  res$penalty <- penalty
  res$tmax_index <- tmax_index
  res$tmin_index <- tmin_index
  return(res)
}
