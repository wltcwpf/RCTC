#' A Function for Combination of Two Rotated Ground Motions
#'
#' This function computes RotD with 180 rotated angles and Geometric Mean of RotD with 90 angles
#' @param data1 The first horizontal component of acceleration ground motion, in g
#' @param data2 The second horizontal component of acceleration ground motion, in g
#' @param period_t A array of oscillator periods
#' @param damping Damping ratio
#' @param time_dt Time step of two time series. These two datasets have the same time step
#' @param fraction A fraction for subset selection to speed up calculation. 0 is using all dataset.
#' 0.7 is recommended for RotD50 if the input data are Sinc interpolated (default);
#' 0.5 is recommended if GMRotI50 is also interested; 0.0 should be used regarding for RotD00
#' @param Interpolation_factor The specific value depends on sampling rate, the detailed explanation is
#' described in the PEER report. Users can give their desired value, or can use defult "auto", which
#' will compute \code{Interpolation_factor} according to time step automatically.
#' @return A list is returned with RotD, GMRotD, number of selected subsets points, threshold values, PGARotD (in g), PGVRotD (in cm/s), and PGDRotD (in cm)
#' @keywords GM_RotD_cal
#' @importFrom pracma cumtrapz
#' @export
GM_RotD_cal <- function(data1, data2, period_t, damping, time_dt, fraction = 0.7, Interpolation_factor = "auto"){

  ### Compute for PGA, PGV, PGD
  vel_1 <- pracma::cumtrapz(time_dt*seq(1,length(data1)),data1) * 981 # unit, cm/s
  vel_2 <- pracma::cumtrapz(time_dt*seq(1,length(data2)),data2) * 981 # unit, cm/s
  disp_1 <- pracma::cumtrapz(time_dt*seq(1,length(vel_1)),vel_1) # unit, cm
  disp_2 <- pracma::cumtrapz(time_dt*seq(1,length(vel_2)),vel_2) # unit, cm
  pga_rot <- rep(0, 180)
  pga_gmrot <- rep(0, 90)
  pgv_rot <- rep(0, 180)
  pgv_gmrot <- rep(0, 90)
  pgd_rot <- rep(0, 180)
  pgd_gmrot <- rep(0, 90)
  for(theta in seq(1,90)){
    # acceleration
    Rot_a1 <- data1*cos(theta/180*pi) + data2*sin(theta/180*pi)
    Rot_a2 <- -data1*sin(theta/180*pi) + data2*cos(theta/180*pi)
    pga_rot[theta] <- max(abs(Rot_a1))
    pga_rot[theta+90] <- max(abs(Rot_a2))
    pga_gmrot[theta] <- sqrt(max(abs(Rot_a1))*max(abs(Rot_a2)))
    # velocity
    Rot_v1 <- vel_1*cos(theta/180*pi) + vel_2*sin(theta/180*pi)
    Rot_v2 <- -vel_1*sin(theta/180*pi) + vel_2*cos(theta/180*pi)
    pgv_rot[theta] <- max(abs(Rot_v1))
    pgv_rot[theta+90] <- max(abs(Rot_v2))
    pgv_gmrot[theta] <- sqrt(max(abs(Rot_v1))*max(abs(Rot_v2)))
    # displacement
    Rot_d1 <- disp_1*cos(theta/180*pi) + disp_2*sin(theta/180*pi)
    Rot_d2 <- -disp_1*sin(theta/180*pi) + disp_2*cos(theta/180*pi)
    pgd_rot[theta] <- max(abs(Rot_d1))
    pgd_rot[theta+90] <- max(abs(Rot_d2))
    pgd_gmrot[theta] <- sqrt(max(abs(Rot_d1))*max(abs(Rot_d2)))
  }

  ### Compute for PSA
  ## Sinc interpolation
  if(Interpolation_factor == 'auto'){
    # interp_factor <- ceiling(log(time_dt/0.001)/log(2))

    # 1) the interpolation factor is based on f_nyq first -- we only care T_nyq,
    # because T less than T_nyq is not affected by interpolation,
    # therefore, dt/dt' = dt/(T_nyq/10) = dt/(dt/5) = 5, -> IF = 8
    # 2) if time_dt is less than the above calcualtion, then we use smaller interpolation
    interp_factor <- min(ceiling(log(time_dt/0.001)/log(2)), 3)
    # if time_dt even smaller than the diresed time step (not common), we do not interpolate
    interp_factor <- max(interp_factor, 0)
    interp_factor <- 2^interp_factor
  }else if(log(Interpolation_factor)%%log(2) != 0){
    print("The Interpolation factor is not a number of power of 2, please try correct it!")
    stop()
  }else{
    interp_factor <- Interpolation_factor
  }
  time_dt <- time_dt/interp_factor # compute new time step
  data1 <- Interpft(data1, interp_factor) # new data after Sinc interpolation
  data2 <- Interpft(data2, interp_factor)
  number <- min(length(data1), length(data2)) # if two data sizes diff, we set both size as the smaller
  data1 <- data1[1:number]
  data2 <- data2[1:number]
  vel_1 <- pracma::cumtrapz(time_dt*seq(1,length(data1)),data1) * 981 # unit, cm/s
  vel_2 <- pracma::cumtrapz(time_dt*seq(1,length(data2)),data2) * 981 # unit, cm/s
  disp_1 <- pracma::cumtrapz(time_dt*seq(1,length(vel_1)),vel_1) # unit, cm
  disp_2 <- pracma::cumtrapz(time_dt*seq(1,length(vel_2)),vel_2) # unit, cm

  # row is period sequence for a given rotated angle, colomn is angle sequence for a given period
  GMRotD <- matrix(nrow = length(seq(1,90)), ncol = length(period_t))
  RotD180 <- matrix(nrow = length(seq(1,180)), ncol = length(period_t))
  # record the number of selected points
  len_period <- length(period_t)
  Num_points <- seq(1,len_period)
  rd_alevel <- seq(1,len_period)
  RD1 <- PS_cal_cpp(data1, period_t, damping, time_dt, 2) # matrix [period by data points], in cm
  RD2 <- PS_cal_cpp(data2, period_t, damping, time_dt, 2) # matrix [period by data points], in cm
  omega = 2*pi/period_t
  for(per_index in seq(1,len_period)){
    length_min <- min(length(RD1[per_index,]), length(RD2[per_index,]))
    rd_subset <- subset_select(RD1[per_index,], RD2[per_index,], fraction, length_min, time_dt, 1)
    rd_alevel[per_index] <- subset_select(RD1[per_index,], RD2[per_index,], fraction, length_min, time_dt, 2)
    Num_points[per_index] <- length(rd_subset)/3
    Rot_rd1 <- rd_subset[1,]
    Rot_rd2 <- rd_subset[2,]
    rot_angle4rot <- seq(1,180)
    for(theta in seq(1,90)){
      RS1 = Rot_rd1*cos(theta/180*pi) + Rot_rd2*sin(theta/180*pi)  # in cm
      RS2 = -Rot_rd1*sin(theta/180*pi) + Rot_rd2*cos(theta/180*pi)  # in cm
      RotD180[theta, per_index] = max(abs(RS1))*omega[per_index]^2 / 981 # convert to g
      RotD180[theta+90, per_index] = max(abs(RS2))*omega[per_index]^2 / 981  # convert to g
      GMRotD[theta,per_index] = sqrt(max(abs(RS1))*max(abs(RS2)))*omega[per_index]^2 / 981 # convert to g
    }
  }

  gm <- list()
  gm$RotD180 <- RotD180
  gm$GMRotD <- GMRotD
  gm$Num_points <- Num_points
  gm$rd_level <- rd_alevel
  gm$pga_rot <- pga_rot
  gm$pgv_rot <- pgv_rot
  gm$pgd_rot <- pgd_rot
  gm$pga_gmrot <- pga_gmrot
  gm$pgv_gmrot <- pgv_gmrot
  gm$pgd_gmrot <- pgd_gmrot
  return(gm)
}
