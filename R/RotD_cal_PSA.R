#' A function to calculate rotated combination of two horizontals for pseudo spectra acceleration
#'
#' This function computes rotated combination of two horizontals for pseudo spectra acceleration
#' @param data1 The first horizontal component of acceleration ground motion, in g
#' @param data2 The second horizontal component of acceleration ground motion, in g
#' @param period A array of oscillator periods. -1 is for PGV, 0 is for PGA.
#' The default is NA, which would use NGA-West2 111 + 2 (PGV and PGA) periods.
#' @param damping Damping ratio.
#' @param combine_index It specifies the value of xx for RotDxx. xx is in percent, e.g., 0 is used for RotD00,
#' 50 is for RotD50, and 100 is for RotD100. The default is 50.
#' @param dt Time step of two time series. We assume that the two time series have the same time step.
#' @param Interpolation_factor The interpolation factor.
#' The specific value depends on sampling rate, the detailed explanation is
#' described in the PEER report. Users can give their desired value, or can use defult "auto", which
#' will compute \code{Interpolation_factor} according to time step automatically.
#' @param fraction A fraction for subset selection to speed up calculation. 0 is using all data.
#' 0.7 is recommended for PSA RotD50 if the input data are Sinc interpolated (default);
#' 0.5 is recommended if PSA GMRotI50 is also interested; 0.0 should be used regarding for PSA RotD00
#' @return A list is returned with RotDxx for the specified \code{combine_index} percentile.
#' The units for returned PGV and PGA/PSA are cm/s and g, respectively.
#' @keywords RotD_cal_PSA
#' @importFrom stats quantile
#' @export

RotD_cal_PSA <- function(data1, data2, period = NA, damping = 0.05, combine_index = 50, dt, Interpolation_factor = 'auto',
                         fraction = 0.7) {

  if (is.na(period[1])) {
    period <- c(-1, 0, 0.010, 0.020, 0.022, 0.025, 0.029, 0.030, 0.032, 0.035, 0.036, 0.040, 0.042,
                0.044, 0.045, 0.046, 0.048, 0.050, 0.055, 0.060, 0.065, 0.067, 0.070, 0.075,
                0.080, 0.085, 0.090, 0.095, 0.100, 0.110, 0.120, 0.130, 0.133, 0.140, 0.150,
                0.160, 0.170, 0.180, 0.190, 0.200, 0.220, 0.240, 0.250, 0.260, 0.280, 0.290,
                0.300, 0.320, 0.340, 0.350, 0.360, 0.380, 0.400, 0.420, 0.440, 0.450, 0.460,
                0.480, 0.500, 0.550, 0.600, 0.650, 0.667, 0.700, 0.750, 0.800, 0.850, 0.900,
                0.950, 1.000, 1.100, 1.200, 1.300, 1.400, 1.500, 1.600, 1.700, 1.800, 1.900,
                2.000, 2.200, 2.400, 2.500, 2.600, 2.800, 3.000, 3.200, 3.400, 3.500, 3.600,
                3.800, 4.000, 4.200, 4.400, 4.600, 4.800, 5.000, 5.500, 6.000, 6.500, 7.000,
                7.500, 8.000, 8.500, 9.000, 9.500, 10.000, 11.000, 12.000, 13.000, 14.000,
                15.000, 20.000)
  }

  res <- GM_RotD_cal(data1 = data1, data2 = data2, period_t = period, damping = damping, time_dt = dt, fraction = fraction,
                     Interpolation_factor = Interpolation_factor)
  PSA_h1 <- c()
  PSA_h2 <- c()
  PSA_RotDxx <- c()

  for (i in 1:length(period)) {

    if (period[i] == -1) {
      # for PGV
      PSA_h1 <- c(PSA_h1, res$pgv_rot[180])
      PSA_h2 <- c(PSA_h2, res$pgv_rot[90])
      PSA_RotDxx <- c(PSA_RotDxx, as.numeric(quantile(res$pgv_rot, probs = combine_index/100)))
    } else if (period[i] == 0) {
      # for PGA
      PSA_h1 <- c(PSA_h1, res$pga_rot[180])
      PSA_h2 <- c(PSA_h2, res$pga_rot[90])
      PSA_RotDxx <- c(PSA_RotDxx, as.numeric(quantile(res$pga_rot, probs = combine_index/100)))
    } else {
      # for PSA
      PSA_h1 <- c(PSA_h1, res$RotD180[180, i])
      PSA_h2 <- c(PSA_h2, res$RotD180[90, i])
      PSA_RotDxx <- c(PSA_RotDxx, as.numeric(quantile(res$RotD180[, i], probs = combine_index/100)))
    }
  }

  output <- data.frame(period = period, PSA_h1 = PSA_h1, PSA_h2 = PSA_h2, PSA_RotDxx = PSA_RotDxx)
  return(output)
}
