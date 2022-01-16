#' A function to calculate standardized Cumulative Absolute Velocity
#'
#' This function computes standardized CAV for the given acceleration ground motion. Note that if the number of
#' data points in every second time interval is not integer (i.e., 1/dt is not an integer),
#' then interpolation is applied at the integer second. Then the returned values
#' are the results after interpolation (so the number of data points could be greater than original).
#' @param acc The acceleration time series, in g
#' @param dt Time step, in second.
#' @param thres The threshold for the standardized CAV, in cm/sec^2. Default is 5 cm/sec^2, denoted as CAV5 for the most
#' commonly used metric
#' @return The calculated standardized CAV and the corresponding time are returned.
#' @examples CAVsd_cal(acc = h1$data, dt = h1$dt)
#' @references Campbell, K.W. and Bozorgnia, Y. (2011). Prediction equations for the standardized version of cumulative
#' absolute velocity as adapted for use in the shutdown of U.S. nuclear power plants. Nuclear Engineering and Design.
#' 241. 2558-2569.
#' @export
#' @importFrom pracma cumtrapz
CAVsd_cal <- function(acc, dt, thres = 5) {

  if (((1/dt) %% 1) == 0) {
    # there are integer number of data points within 1 sec interval, so no need to interpolate
    acc_use <- acc
    time_use <- seq(0, length(acc) - 1) * dt
  } else {
    # the number of data points within 1 sec is not exactly integer, so linear interpolation is applied
    # at each integer sec
    time <- seq(0, length(acc) - 1) * dt
    time_use <- sort(unique(c(time, seq(0, max(time), by = 1))))
    acc_use <- rep(0, length(time_use))
    idx_org <- which(time_use %in% time)
    idx_intp <- which(!(time_use %in% time))
    acc_use[idx_org] <- acc
    acc_use[idx_intp] <- approx(x = time, y = acc, xout = time_use[idx_intp])$y
  }

  acc_cav <- rep(0, length(acc_use))
  for (i_0 in seq(0, max(time_use))) {
    idx <- which(time_use >= i_0 & time_use <= i_0 + 1)
    if (length(idx) > 1) {
      # if there are more than one data points within [i_0, i_0 + 1] time interval
      pga <- max(abs(acc_use[idx]))
      if (pga*9.81*100 < thres) {
        # exclude the data points less than thres
        if (i_0 == 0) {
          acc_cav[idx] <- 0
        } else {
          acc_cav[idx] <- acc_cav[min(idx)]
        }
      } else {
        # include the data points
        if (i_0 == 0) {
          acc_cav[idx] <- cumtrapz(x = time_use[idx], y = abs(acc_use[idx]))
        } else {
          acc_cav[idx] <- acc_cav[min(idx)] + cumtrapz(x = time_use[idx], y = abs(acc_use[idx]))
        }
      }
    } else {
      # if there is only one data point within [i_0, i_0 + 1] time interval
      acc_cav[idx] <- acc_cav[idx - 1]
    }
  }

  return(data.frame(time = time_use, CAVsd = acc_cav))
}


#' A function to calculate Cumulative Absolute Velocity (CAV)
#'
#' This function computes CAV for the given acceleration ground motion
#' @param acc The acceleration time series, in g
#' @param dt Time step, in second.
#' @return The calculated CAV and the corresponding time are returned.
#' @examples CAV_cal(acc = h1$data, dt = h1$dt)
#' @references Campbell, K. W. and Bozorgnia, Y., (2012). Comparison of Ground Motion Prediction
#' Equations for Arias Intensity and Cumulative Absolute Velocity Developed Using a
#' Consistent Database and Functional Form. Earthquake Spectra, 28(3), 931-941
#' @export
#' @importFrom pracma cumtrapz
CAV_cal <- function(acc, dt) {
  return(data.frame(time = seq(0, length(acc) - 1) * dt,
                    CAV = cumtrapz(x = seq(0, length(acc) - 1) * dt, y = abs(acc))))
}
