#' A function to calculate significant duration of rotated combination of two horizontals
#'
#' This function computes rotated combination of two horizontals for significant duration based on AI
#' @param data1 The first horizontal component of acceleration ground motion, in g
#' @param data2 The second horizontal component of acceleration ground motion, in g
#' @param dt Time step of two time series, in second. We assume that the two time series have the same time step.
#' @param combine_index It specifies the value of xx for RotDxx. xx is in percent, e.g., 0 is used for RotD00,
#' 50 is for RotD50, and 100 is for RotD100. The default is 50.
#' @param start An array of the percentile of start time for the significant duration calculation
#' @param end An array of the percentile of end time for the significant duration calculation. The length should be
#' the same as \code{start}
#' @return A data frame is returned which contains three elements, the start time for each percentiles defined by \code{start}
#' the end time for each percentiles defined by \code{end}, and the corresponding significant duration
#' @examples RotD_cal_SigDur(data1 = h1$data, data2 = h2$data, dt = h1$dt, combine_index = 50)
#' @references Bozorgnia, Y., and Campbell, K. W., (2004). Engineering characterization of ground motion,
#' in Earthquake Engineering, From Engineering Seismology to Performance-Based Engineering, Y. Bozorgnia and
#' V. V. Bertero (eds.), CRC Press, Boca Rotan, FL, pp. 5-1-5-74.
#' @export
#' @importFrom stats quantile

RotD_cal_SigDur <- function(data1, data2, dt, combine_index, start = c(0.05, 0.05, 0.1, 0.1, 0.15),
                            end = c(0.95, 0.9, 0.9, 0.8, 0.75)) {
  start_rot <- matrix(0, nrow = 180, ncol = length(start))
  end_rot <- matrix(0, nrow = 180, ncol = length(start))

  for(theta in seq(1,90)){
    # acceleration
    Rot_a1 <- data1*cos(theta/180*pi) + data2*sin(theta/180*pi)
    Rot_a2 <- -data1*sin(theta/180*pi) + data2*cos(theta/180*pi)

    start_rot[theta, ] <- AI_times(Rot_a1, dt, percentiles = start)
    start_rot[theta + 90, ] <- AI_times(Rot_a1, dt, percentiles = start)

    end_rot[theta, ] <- AI_times(Rot_a1, dt, percentiles = end)
    end_rot[theta + 90, ] <- AI_times(Rot_a1, dt, percentiles = end)
  }

  sig_dur_start <- apply(start_rot, 2, function(x) {as.numeric(quantile(x, probs = combine_index/100))})
  sig_dur_end <- apply(end_rot, 2, function(x) {as.numeric(quantile(x, probs = combine_index/100))})
  sig_dur <- ifelse(sig_dur_end >= sig_dur_start, sig_dur_end - sig_dur_start, NA)

  return(data.frame(start = sig_dur_start, end = sig_dur_end, duration = sig_dur))
}

#' A function to calculate Arias Intensity
#'
#' This function computes Arias Intensity array
#' @param acc The acceleration time series, in g
#' @param dt The time step, in second
#' @return An array of calculated Arias Intensity
#' @export
#' @references Bozorgnia, Y., and Campbell, K. W., (2004). Engineering characterization of ground motion,
#' in Earthquake Engineering, From Engineering Seismology to Performance-Based Engineering, Y. Bozorgnia and
#' V. V. Bertero (eds.), CRC Press, Boca Rotan, FL, pp. 5-1-5-74.
#' @importFrom pracma cumtrapz

AI_cal <- function(acc, dt) {
  return(pi / (2 * 9.81) * cumtrapz(x = seq(0, length(acc) - 1) * dt, y = acc^2))
}

#' A function to calculate significant duration
#'
#' This function computes significant duration based on Arias Intensity (AI)
#' @param acc The acceleration time series, in g
#' @param dt The time step, in second
#' @param percentiles The percentiles of normalized AI at which the time is calculated. Default is NA, the funcation will
#' calculate the time for the default percentiles of normalized AI for seq(0.05, 0.95, by = 0.05)
#' @return An array of the time corresponds to each of the percentiles
#' @export
#' @references Bozorgnia, Y., and Campbell, K. W., (2004). Engineering characterization of ground motion,
#' in Earthquake Engineering, From Engineering Seismology to Performance-Based Engineering, Y. Bozorgnia and
#' V. V. Bertero (eds.), CRC Press, Boca Rotan, FL, pp. 5-1-5-74.
#' @importFrom stats approx

AI_times <- function(acc, dt, percentiles = NA) {
  AI_val <- AI_cal(acc, dt)
  AI_norm <- AI_val / max(AI_val)
  time_t <- seq(0, length(acc) - 1) * dt
  times <- approx(x = AI_norm, y = time_t, xout = percentiles, rule = 2, ties = min)$y
  return(times)
}



