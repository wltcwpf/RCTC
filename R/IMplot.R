#' A Function for Intensity Measures Calculation and Plots
#'
#' This is the main function to compute response spectra and generate spreadsheets and plots.
#' The function computes RotD00, RotD50, RotD100, GMRotI50, and other related period independent
#' variables. The outputs will be two cvs files and plots for each station. A summary csv file is also generated.
#' @param inputpath A directory specifies the location of the Inputdata folder where all data are saved.
#' If the data are not downloaded from NGA-West2 database,
#' the file names should be standardized by function \code{nametransfer}
#' @param datatype A string, containing options "ngaw2", "nga", "peer_format", "cosmos", "smc",
#' or "timeseries".
#' "ngaw2" represents the downloaded data from NGA-West2 database,
#'  while "nga" stands for the PEER formatted NGA-Sub and NGA-East data.
#' The reason we distinguish them is that flatfiles for NGA-Sub and NGA-East data
#' are not presently included in the RCTC package and their azimuths and usable frequencies
#' cannot be recognized by RCTC automatically.
#' "peer_format" is for the data follows PEER format but not from NGA database.
#' "cosmos" and "smc" are data formats from
#' Strong-Motion Virtual Data Center and USGS, respectively.
#' "timeseries" consists of a one column array.
#' The array has the time step in the first row of the column and
#' acceleration values below. The unit of "cosmos" and "smc" is cm/sec/sec, while the
#' units are defaulted as g for other data fromats. If the units are not consistent with here,
#' please convert them before using RCTC.
#' Data that are not downloaded from PEER should have files names adjusted as
#' described in Section 3.1, which can be automated by using the nametransfer function.
#' The operation of this function is described further below. If an error occurs that
#' is related to data format, the easiest remedy is to format the data as a column vector
#' ("timeseries" option).
#' @param tmax4penalty_in The maximum usable period for computation of penalty, default value is 10.
#' The specification of the value is only needed if input data are non-NGA-West2.
#' @param tmin4penalty_in The minimum usable period for computation of penalty, default value is 0
#' The specification of the value is only needed if input data are non-NGA-West2.
#' @param combine_index It specifies the value of xx for RotDxx. 0 is used for RotD00,
#' 50 is for RotD50, and 100 is for RotD100. The default is 50;
#' @param ang1 The as-recorded azimuth angle of the first horizontal component.
#' It can be obtained automatically from NGA-West2 flatfile if the input data are downloaded
#' from NGA-West2 database. Otherwise, users need to provide this angle.
#' The second component is assumed to be 90 degrees clockwise from this value.
#' If this value is not specified, a default value of zero is used
#' @param damping It is the fraction of critical damping for which the oscillator response is computed
#' (expressed as decimal). The default is 0.05.
#' @param fraction 0.7 is default for RotD50; 0.5 is recommended if GMRotI50 is also interested;
#' 0.0 should be used regarding for RotD00
#' @param Interpolation_factor The specific value depends on sampling rate, the detailed explanation is
#' described in the PEER report. Users can give their desired value, or can use defult "auto", which
#' will compute \code{Interpolation_factor} according to time step automatically.
#' @useDynLib RCTC
#' @importFrom Rcpp sourceCpp
#' @return It generates two csv files and plots for each input station recordings. Two folders are generated in the
#' \code{inputpath} directory, one is for "Outputdata" and another is for "Outputplot"
#' @keywords IMplot
#' @export

IMplot <- function(inputpath='/Users/PFW/Desktop/test/Inputdata', datatype = "ngaw2", tmax4penalty_in=10, tmin4penalty_in=0, combine_index=50,
                   ang1 = 0, damping=0.05, fraction=0.7, Interpolation_factor='auto'){
  ## check if required packages are installed
  list.of.packages <- c("Rcpp")
  new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  if(length(new.packages)) install.packages(new.packages)
  ## check if libraries are loaded
  if(!('package:Rcpp' %in% search())){library(Rcpp)}

  setwd(inputpath)
  setwd("..")
  maindir <- getwd()
  dir.create(file.path(maindir, 'Outputdata')) # create a folder for Outputdata
  outputdatadir <- file.path(maindir, 'Outputdata')
  dir.create(file.path(maindir, 'Outputplot')) # creata a folder for Outputplot
  outputplotdir <- file.path(maindir, 'Outputplot')


  ############### periods ###################
  period_t = c(0.010, 0.020, 0.022, 0.025, 0.029, 0.030, 0.032, 0.035, 0.036, 0.040, 0.042,
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

  ########### penalty calculation function ############
  # input parameters: 1, tmax4penalty -> default maximum penalty period, 10s, only applied if lowest usable frequency from Flatfile is not available
  #                   2, tmin4penalty -> default minimum penalty period, 0s, only applied if highest usable frequency from Flatfile is not available
  #                   3, f_lowest -> the lowest usable frequency, obtaining from NGA Flatfile
  #                   4, f_highest -> the highest usable frequency, obtaining from NGA Flatfile
  # output: 1, tmax4penalty -> the final maximum period for penalty;
  #         2, tmin4penalty -> the final mimimum period for penalty;
  #         3, tmax_index -> the index of tmax4penalty in period_t vector
  #         4, tmin_index -> the index of tmin4penalty in period_t vector
  penalty_fun <- function(tmax4penalty, tmin4penalty, f_lowest, f_highest){
    if(f_lowest != -999){
      tmax4penalty <- 1/f_lowest
    }
    if(f_highest != -999){
      tmin4penalty <- 1/f_highest
    }
    if(max(period_t) < tmax4penalty){
      tmax_index <- length(period_t)
    }else{
      for(i in seq(1, length(period_t))){
        if(period_t[i] == tmax4penalty){
          tmax_index <- i
          break
        }
        if(period_t[i] > tmax4penalty){
          tmax_index <- i -1
          tmax4penalty <- period_t[tmax_index]
          break
        }
      }
    }
    for(i in seq(length(period_t), 1)){
      if(period_t[i] == tmin4penalty){
        tmin_index <- i
        break
      }
      if(period_t[i] < tmin4penalty){
        tmin_index <- i + 1
        tmin4penalty <- period_t[tmin_index]
        break
      }
      if(tmin4penalty == 0){
        tmin_index <- 1
        break
      }
    }
    return(c(tmax4penalty, tmin4penalty, tmax_index, tmin_index))
  }
  ###############################

  ############### process downloaded NGA-West2 data #################################
  if(datatype == 'ngaw2'){
    setwd(inputpath)
    temp = list.files(pattern = "*.AT2")
    if(!length(temp)%%2){            # judge if it is even number of groups
      ## define combined spreadsheet variables ##
      Record_Sequence_Number <- seq(1, length(temp)/2)
      EQID <- seq(1, length(temp)/2)
      pga_combine <- seq(1, length(temp)/2)*0
      pgv_combine <- seq(1, length(temp)/2)*0
      pgd_combine <- seq(1, length(temp)/2)*0
      rotd_combine <- matrix(seq(1, (length(temp)/2)*length(period_t)) ,nrow = length(temp)/2, ncol = length(period_t))*0
      highest_freq_combine <- seq(1, length(temp)/2)
      lowest_freq_combine <- seq(1, length(temp)/2)

      for(i in 1:(length(temp)/2)){
        setwd(inputpath)
        name1 <- basename(temp[2*i-1])
        temp_name1 <- strsplit(name1, split = "[_]")
        name2 <- basename(temp[2*i])
        temp_name2 <- strsplit(name2, split = "[_]")
        if(temp_name1[[1]][1] != temp_name2[[1]][1]){      # judge if there is a coupled group for name1
          print(c("You should put the coupled-group in the Inputdata folder!", name1, "is not coupled. Quit!"))
          stop()
        }
        ## get the azimuth for data1
        temp_rsn <- as.numeric(gsub("[^0-9\\.]","", temp_name1[[1]][1]))
        temp_ind <- which(NGAW2$Record.Sequence.Number == temp_rsn)
        if(grepl(tail(temp_name1[[1]],1), as.character(NGAW2[temp_ind, 9]))){
          ang1 <- as.numeric(NGAW2[temp_ind, 11])
        } else if(grepl(tail(temp_name1[[1]],1), as.character(NGAW2[temp_ind, 10]))){
          ang1 <- as.numeric(NGAW2[temp_ind, 12])
        } else{
          ang1 <- 0
        }

        assign(temp[2*i-1], head_dt1 <- scan(temp[2*i-1], sep = ",", skip = 3, nlines = 1, what = character()))
        dt1 <- as.numeric(gsub("[^0-9\\.]","", head_dt1[2]))
        assign(temp[2*i], head_dt2 <- scan(temp[2*i], sep = ",", skip = 3, nlines = 1, what = character()))
        dt2 <- as.numeric(gsub("[^0-9\\.]","", head_dt2[2]))
        if(dt1 != dt2){
          print("Samples per second are different! Quit!")
          stop()
        }
        assign(temp[2*i-1], data1 <- scan(temp[2*i-1], sep = "", skip = 4))
        assign(temp[2*i], data2 <- scan(temp[2*i], sep = "", skip = 4))

        m <- gregexpr('[0-9]+',temp_name1[[1]][1])
        number <- regmatches(temp_name1[[1]][1],m)
        Record_Sequence_Number[i] <- as.numeric(number[[1]][1])
        EQID[i] <- NGAW2[match(Record_Sequence_Number[i], NGAW2[,1]), 2]

        #################### usable frequency ##############################
        ### combined lowest usable frequenct is the higher LUFreq of two components
        lowest_usable_freq <- pmax(NGAW2[match(Record_Sequence_Number[i], NGAW2[,1]), 6],
                                   NGAW2[match(Record_Sequence_Number[i], NGAW2[,1]), 7])
        lowest_freq_combine[i] <- lowest_usable_freq
        ### combined highest usable frequenct is the lower HUFreq of two components
        higher <- pmax(NGAW2[match(Record_Sequence_Number[i], NGAW2[,1]), 3],
                       NGAW2[match(Record_Sequence_Number[i], NGAW2[,1]), 4])
        lower <- pmin(NGAW2[match(Record_Sequence_Number[i], NGAW2[,1]), 3],
                      NGAW2[match(Record_Sequence_Number[i], NGAW2[,1]), 4])
        if(higher == -999 || higher == 0){
          highest_usable_freq <- -999
        } else if(lower == -999 || lower == 0){
          highest_usable_freq <- higher/NGAW2[match(Record_Sequence_Number[i], NGAW2[,1]), 5]
        } else{
          highest_usable_freq <- lower/NGAW2[match(Record_Sequence_Number[i], NGAW2[,1]), 5]
        }
        highest_freq_combine[i] <- highest_usable_freq

        penalty_res <- penalty_fun(tmax4penalty_in, tmin4penalty_in, lowest_usable_freq, highest_usable_freq)
        tmax4penalty <- penalty_res[1]
        tmin4penalty <- penalty_res[2]
        tmax_index <- penalty_res[3]
        tmin_index <- penalty_res[4]

        ############################# Process data Begin ###################################
        ######### Sinc interpolate data ###############
        if(Interpolation_factor == 'auto'){
          # interp_factor <- ceiling(log(dt1/0.001)/log(2))

          # 1) the interpolation factor is based on f_nyq first -- we only care T_nyq,
          # because T less than T_nyq is not affected by interpolation,
          # therefore, dt/dt' = dt/(T_nyq/10) = dt/(dt/5) = 5, -> IF = 8
          # 2) if dt1 is less than the above calcualtion, then we use smaller interpolation
          interp_factor <- min(ceiling(log(dt1/0.001)/log(2)), 3)
          # if dt1 even smaller than the diresed time step (not common), we do not interpolate
          interp_factor <- max(interp_factor, 0)
          interp_factor <- 2^interp_factor
        }else if(log(Interpolation_factor)%%log(2) != 0){
          print("The Interpolation factor is not a number of power of 2, please try correct it!")
          stop()
        }else{
          interp_factor <- Interpolation_factor
        }
        # dt <- dt1  # original time step
        # data1_1 <- data1
        # data2_1 <- data2
        dt1 <- dt1/interp_factor # compute new time step
        dt2 <- dt2/interp_factor
        data1 <- Interpft(data1, interp_factor) # new data after Sinc interpolation
        data2 <- Interpft(data2, interp_factor)
        number <- min(length(data1), length(data2)) # if two data sizes diff, we set both size as the smaller
        data1 <- data1[1:number]
        data2 <- data2[1:number]

        ############ Output PGA, PGV, PGD ###############
        ######### unit of original accel is g
        pga_1_ar <- max(abs(data1)) # unit, g
        pga_2_ar <- max(abs(data2)) # unit, g
        vel_1 <- cumtrapz(dt1*seq(1,length(data1)),data1) * 981 # unit, cm/s
        vel_2 <- cumtrapz(dt2*seq(1,length(data2)),data2) * 981 # unit, cm/s
        pgv_1_ar <- max(abs(vel_1))
        pgv_2_ar <- max(abs(vel_2))
        disp_1 <- cumtrapz(dt1*seq(1,length(vel_1)),vel_1) # unit, cm
        disp_2 <- cumtrapz(dt2*seq(1,length(vel_2)),vel_2) # unit, cm
        pgd_1_ar <- max(abs(disp_1))
        pgd_2_ar <- max(abs(disp_2))

        ############## Rotation values calculation ##############
        # subset slection
        length_min <- min(length(data1), length(data2))
        a_subset <- subset_select(data1, data2, fraction, length_min, dt1, 1)
        v_subset <- subset_select(vel_1, vel_2, fraction, length_min, dt1, 1)
        d_subset <- subset_select(disp_1, disp_2, fraction, length_min, dt1, 1)
        # define vectors
        pga_rot <- seq(1, 180)
        pga_gmrot <- seq(1,90)
        pgv_rot <- seq(1, 180)
        pgv_gmrot <- seq(1,90)
        pgd_rot <- seq(1, 180)
        pgd_gmrot <- seq(1,90)
        # calculation
        for(theta in seq(1,90)){
          # acceleration
          Rot_a1 <- a_subset[1,]*cos(theta/180*pi) + a_subset[2,]*sin(theta/180*pi)
          Rot_a2 <- -a_subset[1,]*sin(theta/180*pi) + a_subset[2,]*cos(theta/180*pi)
          pga_rot[theta] <- max(abs(Rot_a1))
          pga_rot[theta+90] <- max(abs(Rot_a2))
          pga_gmrot[theta] <- sqrt(max(abs(Rot_a1))*max(abs(Rot_a2)))
          # velocity
          Rot_v1 <- v_subset[1,]*cos(theta/180*pi) + v_subset[2,]*sin(theta/180*pi)
          Rot_v2 <- -v_subset[1,]*sin(theta/180*pi) + v_subset[2,]*cos(theta/180*pi)
          pgv_rot[theta] <- max(abs(Rot_v1))
          pgv_rot[theta+90] <- max(abs(Rot_v2))
          pgv_gmrot[theta] <- sqrt(max(abs(Rot_v1))*max(abs(Rot_v2)))
          # displacement
          Rot_d1 <- d_subset[1,]*cos(theta/180*pi) + d_subset[2,]*sin(theta/180*pi)
          Rot_d2 <- -d_subset[1,]*sin(theta/180*pi) + d_subset[2,]*cos(theta/180*pi)
          pgd_rot[theta] <- max(abs(Rot_d1))
          pgd_rot[theta+90] <- max(abs(Rot_d2))
          pgd_gmrot[theta] <- sqrt(max(abs(Rot_d1))*max(abs(Rot_d2)))
        }
        ## acceleration
        pga_rot00 <- min(pga_rot)
        pga_rot50 <- median(pga_rot)
        pga_rot100 <- max(pga_rot)
        pga_rot00ang <- which.min(pga_rot)
        pga_rot100ang <- which.max(pga_rot)
        pga_gmrot50 <- median(pga_gmrot)
        pga_gmrot100 <- max(pga_gmrot)
        pga_gmrot100ang <- which.max(pga_gmrot)
        ## velocity
        pgv_rot00 <- min(pgv_rot)
        pgv_rot50 <- median(pgv_rot)
        pgv_rot100 <- max(pgv_rot)
        pgv_rot00ang <- which.min(pgv_rot)
        pgv_rot100ang <- which.max(pgv_rot)
        pgv_gmrot50 <- median(pgv_gmrot)
        pgv_gmrot100 <- max(pgv_gmrot)
        pgv_gmrot100ang <- which.max(pgv_gmrot)
        ## displacement
        pgd_rot00 <- min(pgd_rot)
        pgd_rot50 <- median(pgd_rot)
        pgd_rot100 <- max(pgd_rot)
        pgd_rot00ang <- which.min(pgd_rot)
        pgd_rot100ang <- which.max(pgd_rot)
        pgd_gmrot50 <- median(pgd_gmrot)
        pgd_gmrot100 <- max(pgd_gmrot)
        pgd_gmrot100ang <- which.max(pgd_gmrot)

        ################# Applying PS_cal function to get PSA,PSV #############
        ########################## plug in the parameters ########################
        # plug in the parameters for two components
        # PSA calculation is from textbook of "Dynamics of Structure" 3rd edition by Chopra on page 168-169.
        # PS_cal_c.cpp function, the last inpur parameter decides what type of results return: 1 for PSA; 2 for relative displacements
        ##### compute pseudo spectral acceleration
        PSA1 <- PS_cal_cpp(data1, period_t, damping, dt1, 1) # first row is absolute acceleration of oscillator, second row is psa
        PSA2 <- PS_cal_cpp(data2, period_t, damping, dt1, 1) # first row is absolute acceleration of oscillator, second row is psa
        PSA_gm_ar <- sqrt(PSA1[2,]*PSA2[2,])
        PSA_larger_ar <- seq(1, length(period_t))
        for(j in seq(1, length(PSA_larger_ar))){
          PSA_larger_ar[j] <- max(PSA1[2,j], PSA2[2,j])
        }

        ############## Applying GMRotD_cal function to get PSA for different periods for each rotated angle ####
        ########################## plug in the parameters ########################
        ##########################################################################
        results <- GM_RotD_cal(data1, data2, period_t, damping, dt1, fraction)  # the result is list, 1st is RotD180, 2nd is GMRotD
        RotD180 <- results[[1]]  # it is a list
        RotD180 <- matrix(unlist(RotD180), nrow = 180) # transform it into matrix
        GMRotD <- results[[2]]  # it is a list
        GMRotD <- matrix(unlist(GMRotD), nrow = 90) # transform it into matrix
        Num_points <- results[[3]]  # it is a list
        Num_points <- c(unlist(Num_points)) # transform it into matrix
        rd_alevel <- results[[4]]

        ##################### Applying GMRotI_cal function to get PSA for different periods for only one desired rotated angle
        ########################## plug in the parameters ########################
        ##########################################################################
        Result <- GMRotI50_cal(GMRotD, period_t, tmax4penalty, tmin4penalty)
        GMRotI50 <- Result[[1]]
        GMRotI50_ang <- Result[[2]]
        penalty <- Result[[3]]
        # acceleration
        pga_gmrotI50 <- pga_gmrot[GMRotI50_ang]
        # velocity
        pgv_gmrotI50 <- pgv_gmrot[GMRotI50_ang]
        # displacement
        pgd_gmrotI50 <- pgd_gmrot[GMRotI50_ang]

        ##################### Applying RotD_cal function to get PSA for different periods for each rotated angle
        ########################## plug in the parameters ########################
        ##########################################################################
        RotD00 <- apply(RotD180, 2, min)
        RotD100 <- apply(RotD180, 2, max)
        RotD50 <- apply(RotD180, 2, median)
        RotD_t <- t(RotD180) # transpose of RotD matrix, it is used to get the max index
        RotD100_ang <- max.col(RotD_t)
        RotD00_ang <- max.col(-RotD_t)

        ############################# Process data End ###################################
        ###############################################################################
        if(combine_index == 0){
          pga_combine[i] <- pga_rot00
          pgv_combine[i] <- pgv_rot00
          pgd_combine[i] <- pgd_rot00
          rotd_combine[i,] <- RotD00
        } else if(combine_index == 50){
          pga_combine[i] <- pga_rot50
          pgv_combine[i] <- pgv_rot50
          pgd_combine[i] <- pgd_rot50
          rotd_combine[i,] <- RotD50
        } else{
          pga_combine[i] <- pga_rot100
          pgv_combine[i] <- pgv_rot100
          pgd_combine[i] <- pgd_rot100
          rotd_combine[i,] <- RotD100
        }

        #################### Output calculation results ####################
        # ## check if required packages are installed
        # list.of.packages <- c("xlsx", "rJava", "xlsxjars")
        # new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
        # if(length(new.packages)) install.packages(new.packages)
        # ## check if libraries are loaded
        # if(!('package:xlsx' %in% search())){library(xlsx)}
        # if(!('package:rJava' %in% search())){library(rJava)}
        # if(!('package:xlsxjars' %in% search())){library(xlsxjars)}
        PSA_1 <- PSA1[2,]
        Absolute_acc_1 <- PSA1[1,] # absolute acceleration of oscillator
        PSA_2 <- PSA2[2,]
        Absolute_acc_2 <- PSA2[1,]
        npts1 <- length(data1)
        npts2 <- length(data2)
        n_a_subset <- length(a_subset[1,])
        n_v_subset <- length(v_subset[1,])
        n_d_subset <- length(d_subset[1,])
        Period <- c("T0.010S","T0.020S","T0.022S","T0.025S","T0.029S","T0.030S",
                    "T0.032S","T0.035S","T0.036S","T0.040S","T0.042S","T0.044S",
                    "T0.045S","T0.046S","T0.048S","T0.050S","T0.055S","T0.060S",
                    "T0.065S","T0.067S","T0.070S","T0.075S","T0.080S","T0.085S",
                    "T0.090S","T0.095S","T0.100S","T0.110S","T0.120S","T0.130S",
                    "T0.133S","T0.140S","T0.150S","T0.160S","T0.170S","T0.180S",
                    "T0.190S","T0.200S","T0.220S","T0.240S","T0.250S","T0.260S",
                    "T0.280S","T0.290S","T0.300S","T0.320S","T0.340S","T0.350S",
                    "T0.360S","T0.380S","T0.400S","T0.420S","T0.440S","T0.450S",
                    "T0.460S","T0.480S","T0.500S","T0.550S","T0.600S","T0.650S",
                    "T0.667S","T0.700S","T0.750S","T0.800S","T0.850S","T0.900S",
                    "T0.950S","T1.000S","T1.100S","T1.200S","T1.300S","T1.400S",
                    "T1.500S","T1.600S","T1.700S","T1.800S","T1.900S","T2.000S",
                    "T2.200S","T2.400S","T2.500S","T2.600S","T2.800S","T3.000S",
                    "T3.200S","T3.400S","T3.500S","T3.600S","T3.800S","T4.000S",
                    "T4.200S","T4.400S","T4.600S","T4.800S","T5.000S","T5.500S",
                    "T6.000S","T6.500S","T7.000S","T7.500S","T8.000S","T8.500S",
                    "T9.000S","T9.500S","T10.000S","T11.000S","T12.000S","T13.000S",
                    "T14.000S","T15.000S","T20.000S")

        # compute real azimuth and set angle in [0,360]
        RotD00_ang <- RotD00_ang + ang1
        RotD100_ang <- RotD100_ang + ang1
        RotD00_ang[which(RotD00_ang > 360)] <- RotD00_ang[which(RotD00_ang > 360)] - 360
        RotD00_ang[which(RotD00_ang < 0)] <- RotD00_ang[which(RotD00_ang < 0)] + 360
        RotD100_ang[which(RotD100_ang > 360)] <- RotD100_ang[which(RotD100_ang > 360)] - 360
        RotD100_ang[which(RotD100_ang < 0)] <- RotD100_ang[which(RotD100_ang < 0)] + 360

        frame.data1 <- data.frame(Period, PSA_1, Absolute_acc_1, PSA_2, Absolute_acc_2,
                                  PSA_gm_ar, PSA_larger_ar, GMRotI50, RotD00, RotD50,
                                  RotD100, RotD00_ang, RotD100_ang, Num_points)
        sub_frame.data1 <- t(frame.data1[,-1])
        frame.data1 <- rbind(t(frame.data1[1]), sub_frame.data1)
        Name <- c("npts1","npts2","tmax4penalty","tmin4penalty","damping","Lowest usable freq",
                  "Highest usable freq", "GMRotI50angle", "PGA_GMRotI50", "PGV_GMRotI50",
                  "PGD_GMRotI50", "PGA_GMRot50", "PGV_GMRot50", "PGA_GMRot100angle",
                  "PGA_GMRot100", "PGV_GMRot100angle", "PGV_GMRot100", "PGA_1", "PGA_2",
                  "n_a_subset","PGV_1", "PGV_2", "n_v_subset","PGD_1", "PGD_2", "n_d_subset",
                  "PGA_Rot00", "PGA_Rot50", "PGA_Rot100", "PGA_Rot00angle", "PGA_Rot100angle",
                  "PGV_Rot00", "PGV_Rot50", "PGV_Rot100", "PGV_Rot00angle", "PGV_Rot100angle",
                  "PGD_Rot00", "PGD_Rot50", "PGD_Rot100", "PGD_Rot00angle", "PGD_Rot100angle")

        # compute real azimuth and set angel in [0,360]
        GMRotI50_ang <- GMRotI50_ang + ang1
        if(GMRotI50_ang > 360) {GMRotI50_ang <- GMRotI50_ang - 360}
        if(GMRotI50_ang < 0) {GMRotI50_ang <- GMRotI50_ang + 360}
        pga_gmrot100ang <- pga_gmrot100ang + ang1
        if(pga_gmrot100ang > 360) {pga_gmrot100ang <- pga_gmrot100ang - 360}
        if(pga_gmrot100ang < 0) {pga_gmrot100ang <- pga_gmrot100ang + 360}
        pgv_gmrot100ang <- pgv_gmrot100ang + ang1
        if(pgv_gmrot100ang > 360) {pgv_gmrot100ang <- pgv_gmrot100ang - 360}
        if(pgv_gmrot100ang < 0) {pgv_gmrot100ang <- pgv_gmrot100ang + 360}
        pga_rot00ang <- pga_rot00ang + ang1
        if(pga_rot00ang > 360) {pga_rot00ang <- pga_rot00ang - 360}
        if(pga_rot00ang < 0) {pga_rot00ang <- pga_rot00ang + 360}
        pga_rot100ang <- pga_rot100ang + ang1
        if(pga_rot100ang > 360) {pga_rot100ang <- pga_rot100ang - 360}
        if(pga_rot100ang < 0) {pga_rot100ang <- pga_rot100ang + 360}
        pgv_rot00ang <- pgv_rot00ang + ang1
        if(pgv_rot00ang > 360) {pgv_rot00ang <- pgv_rot00ang - 360}
        if(pgv_rot00ang < 0) {pgv_rot00ang <- pgv_rot00ang + 360}
        pgv_rot100ang <- pgv_rot100ang + ang1
        if(pgv_rot100ang > 360) {pgv_rot100ang <- pgv_rot100ang - 360}
        if(pgv_rot100ang < 0) {pgv_rot100ang <- pgv_rot100ang + 360}
        pgd_rot00ang <- pgd_rot00ang + ang1
        if(pgd_rot00ang > 360) {pgd_rot00ang <- pgd_rot00ang - 360}
        if(pgd_rot00ang < 0) {pgd_rot00ang <- pgd_rot00ang + 360}
        pgd_rot100ang <- pgd_rot100ang + ang1
        if(pgd_rot100ang > 360) {pgd_rot100ang <- pgd_rot100ang - 360}
        if(pgd_rot100ang < 0) {pgd_rot100ang <- pgd_rot100ang + 360}

        Value <- c(npts1, npts2, tmax4penalty, tmin4penalty,damping, lowest_usable_freq, highest_usable_freq,
                   GMRotI50_ang, pga_gmrotI50, pgv_gmrotI50, pgd_gmrotI50, pga_gmrot50,
                   pgv_gmrot50, pga_gmrot100ang, pga_gmrot100, pgv_gmrot100ang, pgv_gmrot100,
                   pga_1_ar, pga_2_ar,n_a_subset, pgv_1_ar, pgv_2_ar, n_v_subset, pgd_1_ar,
                   pgd_2_ar, n_d_subset, pga_rot00, pga_rot50, pga_rot100, pga_rot00ang,
                   pga_rot100ang, pgv_rot00, pgv_rot50, pgv_rot100, pgv_rot00ang, pgv_rot100ang,
                   pgd_rot00, pgd_rot50, pgd_rot100, pgd_rot00ang, pgd_rot100ang)
        frame.data2 <- data.frame(Name, Value)
        filename1 <- paste(temp_name1[[1]][1], "_dep.csv", sep = "")
        filename2 <- paste(temp_name1[[1]][1], "_indep.csv", sep = "")
        setwd(outputdatadir)
        write.table(frame.data1, file = filename1, col.names = FALSE, sep = ',')
        write.csv(frame.data2, file = filename2, row.names = FALSE)
        # write.xlsx(frame.data1, file = filename, sheetName = "Calc Period-dependent", col.names = FALSE)
        # write.xlsx(frame.data2, file = filename, sheetName = "Calc Period-independent", row.names = FALSE, append = TRUE)


        #################### Plot curves  ########################################
        setwd(outputplotdir)
        png(filename = paste(temp_name1[[1]][1],"_PSA.png", sep = ""), width = 1280, height = 768)
        par(oma=c(2,2,2,2))
        par(mar=c(5,5,4,2) + 0.1)
        plot(c(0.01,20), c(0,max(max(PSA1,PSA2))), log = "x", type = "n", xlab = "Period(s)", ylab = "PSA(g)", main = "PSA versus Period", cex.main=2, cex.lab=2, cex.axis=2)
        lines(period_t, PSA_1, lty = 2, lwd = 2, col = "red")
        lines(period_t, PSA_2, lty = 3, lwd = 2, col = "blue")
        legend(5,max(max(PSA1,PSA2)), c("PSA_H1","PSA_H2"), lty=c(2,3), lwd=c(2,2),col=c("red","blue"), cex = 2)
        dev.off()
        png(filename = paste(temp_name1[[1]][1],"_GMRotI50.png", sep = ""), width = 1280, height = 768)
        par(oma=c(2,2,2,2))
        par(mar=c(5,5,4,2) + 0.1)
        plot(period_t, GMRotI50, log = "x", type = "l", lwd = 2, xlab = "Period(s)", ylab = "GMRotI50(g)", main = "GMRotI50 versus Period", cex.lab=2, cex.axis=2, cex.main=2)
        dev.off()
        png(filename = paste(temp_name1[[1]][1],"_RotD.png", sep = ""), width = 1280, height = 768)
        par(oma=c(2,2,2,2))
        par(mar=c(5,5,4,2) + 0.1)
        plot(c(0.01,20), c(0,max(RotD100)), log = "x", type = "n", xlab = "Period(s)", ylab = "RotDxx(g)", main = "RotDxx versus Period", cex.lab=2, cex.axis=2, cex.main=2)
        lines(period_t, RotD00,  lty = 2, lwd = 2, col = 2)
        lines(period_t, RotD50,  lty = 3, lwd = 2, col = 3)
        lines(period_t, RotD100, lty = 4, lwd = 2, col = 4)
        legend(5, max(RotD100), c("RotD00","RotD50","RotD100"), lty = c(2,3,4),  lwd = c(2,2,2), col = c(2,3,4), cex = 2)
        dev.off()

      }

    } else{
      print("You should put the couple groups in the Inputdata folder, we detect you have odd number of time series")
      stop()
    }

    ###################### Create a combinated spreadsheet ########################
    ###############################################################################
    name_com <- c('Record_Sequence_Number', 'EQID', 'Lowest_Usable_Freq', 'Highest_Usable_Freq')
    frame.data_com1 <- data.frame(Record_Sequence_Number, EQID, lowest_freq_combine, highest_freq_combine)
    frame.data_com1 <- rbind(name_com, frame.data_com1)
    name_com <- c(paste('PGA_Rot', combine_index, sep = ''), paste('PGV_Rot', combine_index, sep = ''), paste('PGD_Rot', combine_index, sep = ''), Period)
    value_com <- cbind(pga_combine, pgv_combine, pgd_combine, rotd_combine)
    frame.data_com2 <- t(data.frame(name_com, t(value_com)))
    frame.data_com <- cbind(frame.data_com1, frame.data_com2)
    filename <- paste('summary', ".csv", sep = "")
    setwd(outputdatadir)
    write.table(frame.data_com, file = filename, col.names = FALSE, row.names = FALSE, sep = ',')


  }


  ############### process PEER formated data #################################
  else if(datatype == 'peer_format'){
    setwd(inputpath)
    temp = list.files()
    if(!length(temp)%%2){            # judge if it is even number of groups
      ## define combined spreadsheet variables ##
      Station <- seq(1, length(temp)/2)
      pga_combine <- seq(1, length(temp)/2)*0
      pgv_combine <- seq(1, length(temp)/2)*0
      pgd_combine <- seq(1, length(temp)/2)*0
      rotd_combine <- matrix(seq(1, (length(temp)/2)*length(period_t)) ,nrow = length(temp)/2, ncol = length(period_t))*0

      for(i in 1:(length(temp)/2)){
        setwd(inputpath)
        name1 <- basename(temp[2*i-1])
        temp_name1 <- strsplit(name1, split = "_")
        name2 <- basename(temp[2*i])
        temp_name2 <- strsplit(name2, split = "_")
        if(temp_name1[[1]][1] != temp_name2[[1]][1]){      # judge if there is a coupled group for name1
          print(c("You should put the coupled-group in the Inputdata folder!", name1, "is not coupled. Quit!"))
          stop()
        }
        station_name <- paste(temp_name1[[1]][1], temp_name1[[1]][2], sep = "_")
        Station[i] <- station_name
        assign(temp[2*i-1], head_dt1 <- scan(temp[2*i-1], sep = ",", skip = 3, nlines = 1, what = character()))
        dt1 <- as.numeric(gsub("[^0-9\\.]","", head_dt1[2]))
        assign(temp[2*i], head_dt2 <- scan(temp[2*i], sep = ",", skip = 3, nlines = 1, what = character()))
        dt2 <- as.numeric(gsub("[^0-9\\.]","", head_dt2[2]))
        if(dt1 != dt2){
          print("Samples per second are different! Quit!")
          stop()
        }
        assign(temp[2*i-1], data1 <- scan(temp[2*i-1], sep = "", skip = 4))
        assign(temp[2*i], data2 <- scan(temp[2*i], sep = "", skip = 4))

        ############################# Process data Begin ###################################
        ######### Sinc interpolate data ###############
        if(Interpolation_factor == 'auto'){
          # interp_factor <- ceiling(log(dt1/0.001)/log(2))

          # 1) the interpolation factor is based on f_nyq first -- we only care T_nyq,
          # because T less than T_nyq is not affected by interpolation,
          # therefore, dt/dt' = dt/(T_nyq/10) = dt/(dt/5) = 5, -> IF = 8
          # 2) if dt1 is less than the above calcualtion, then we use smaller interpolation
          interp_factor <- min(ceiling(log(dt1/0.001)/log(2)), 3)
          # if dt1 even smaller than the diresed time step (not common), we do not interpolate
          interp_factor <- max(interp_factor, 0)
          interp_factor <- 2^interp_factor
        }else if(log(Interpolation_factor)%%log(2) != 0){
          print("The Interpolation factor is not a number of power of 2, please try correct it!")
          stop()
        }else{
          interp_factor <- Interpolation_factor
        }

        dt1 <- dt1/interp_factor # compute new time step
        dt2 <- dt2/interp_factor
        data1 <- Interpft(data1, interp_factor) # new data after Sinc interpolation
        data2 <- Interpft(data2, interp_factor)
        number <- min(length(data1), length(data2)) # if two data sizes diff, we set both size as the smaller
        data1 <- data1[1:number]
        data2 <- data2[1:number]

        ############ Output PGA, PGV, PGD ###############
        ########## unit of original accel is g
        pga_1_ar <- max(abs(data1)) # unit, g
        pga_2_ar <- max(abs(data2)) # unit, g
        vel_1 <- cumtrapz(dt1*seq(1,length(data1)),data1) * 981 # unit, cm/s
        vel_2 <- cumtrapz(dt2*seq(1,length(data2)),data2) * 981 # unit, cm/s
        pgv_1_ar <- max(abs(vel_1))
        pgv_2_ar <- max(abs(vel_2))
        disp_1 <- cumtrapz(dt1*seq(1,length(vel_1)),vel_1) # unit, cm
        disp_2 <- cumtrapz(dt2*seq(1,length(vel_2)),vel_2) # unit, cm
        pgd_1_ar <- max(abs(disp_1))
        pgd_2_ar <- max(abs(disp_2))

        ############## Rotation values calculation ##############
        ######## subset slection #######
        length_min <- min(length(data1), length(data2))
        a_subset <- subset_select(data1, data2, fraction, length_min, dt1, 1)
        v_subset <- subset_select(vel_1, vel_2, fraction, length_min, dt1, 1)
        d_subset <- subset_select(disp_1, disp_2, fraction, length_min, dt1, 1)
        # define vectors
        pga_rot <- seq(1, 180)
        pga_gmrot <- seq(1,90)
        pgv_rot <- seq(1, 180)
        pgv_gmrot <- seq(1,90)
        pgd_rot <- seq(1, 180)
        pgd_gmrot <- seq(1,90)
        # calculation
        for(theta in seq(1,90)){
          # acceleration
          Rot_a1 <- a_subset[1,]*cos(theta/180*pi) + a_subset[2,]*sin(theta/180*pi)
          Rot_a2 <- -a_subset[1,]*sin(theta/180*pi) + a_subset[2,]*cos(theta/180*pi)
          pga_rot[theta] <- max(abs(Rot_a1))
          pga_rot[theta+90] <- max(abs(Rot_a2))
          pga_gmrot[theta] <- sqrt(max(abs(Rot_a1))*max(abs(Rot_a2)))
          # velocity
          Rot_v1 <- v_subset[1,]*cos(theta/180*pi) + v_subset[2,]*sin(theta/180*pi)
          Rot_v2 <- -v_subset[1,]*sin(theta/180*pi) + v_subset[2,]*cos(theta/180*pi)
          pgv_rot[theta] <- max(abs(Rot_v1))
          pgv_rot[theta+90] <- max(abs(Rot_v2))
          pgv_gmrot[theta] <- sqrt(max(abs(Rot_v1))*max(abs(Rot_v2)))
          # displacement
          Rot_d1 <- d_subset[1,]*cos(theta/180*pi) + d_subset[2,]*sin(theta/180*pi)
          Rot_d2 <- -d_subset[1,]*sin(theta/180*pi) + d_subset[2,]*cos(theta/180*pi)
          pgd_rot[theta] <- max(abs(Rot_d1))
          pgd_rot[theta+90] <- max(abs(Rot_d2))
          pgd_gmrot[theta] <- sqrt(max(abs(Rot_d1))*max(abs(Rot_d2)))
        }
        ## acceleration
        pga_rot00 <- min(pga_rot)
        pga_rot50 <- median(pga_rot)
        pga_rot100 <- max(pga_rot)
        pga_rot00ang <- which.min(pga_rot)
        pga_rot100ang <- which.max(pga_rot)
        pga_gmrot50 <- median(pga_gmrot)
        pga_gmrot100 <- max(pga_gmrot)
        pga_gmrot100ang <- which.max(pga_gmrot)
        ## velocity
        pgv_rot00 <- min(pgv_rot)
        pgv_rot50 <- median(pgv_rot)
        pgv_rot100 <- max(pgv_rot)
        pgv_rot00ang <- which.min(pgv_rot)
        pgv_rot100ang <- which.max(pgv_rot)
        pgv_gmrot50 <- median(pgv_gmrot)
        pgv_gmrot100 <- max(pgv_gmrot)
        pgv_gmrot100ang <- which.max(pgv_gmrot)
        ## displacement
        pgd_rot00 <- min(pgd_rot)
        pgd_rot50 <- median(pgd_rot)
        pgd_rot100 <- max(pgd_rot)
        pgd_rot00ang <- which.min(pgd_rot)
        pgd_rot100ang <- which.max(pgd_rot)
        pgd_gmrot50 <- median(pgd_gmrot)
        pgd_gmrot100 <- max(pgd_gmrot)
        pgd_gmrot100ang <- which.max(pgd_gmrot)

        ################# Applying PS_cal function to get PSA,PSV #############
        ########################## plug in the parameters ########################
        ##########################################
        # plug in the parameters for two components
        # using C++ to compute to speed up calculation
        # PSA calculation is from textbook of "Dynamics of Structure" 3rd edition by Chopra on page 168-169.
        ##### compute pseudo spectral acceleration
        PSA1 <- PS_cal_cpp(data1, period_t, damping, dt1, 1) # first row is absolute acceleration of oscillator, second row is psa
        PSA2 <- PS_cal_cpp(data2, period_t, damping, dt1, 1) # first row is absolute acceleration of oscillator, second row is psa
        PSA_gm_ar <- sqrt(PSA1[2,]*PSA2[2,])
        PSA_larger_ar <- seq(1, length(period_t))
        for(j in seq(1, length(PSA_larger_ar))){
          PSA_larger_ar[j] <- max(PSA1[2,j], PSA2[2,j])
        }

        ############## Applying GMRotD_cal function to get PSA for different periods for each rotated angle ####
        ########################## plug in the parameters ########################
        ##########################################################################
        results <- GM_RotD_cal(data1, data2, period_t, damping, dt1, fraction)  # the result is list, 1st is RotD180, 2nd is GMRotD
        RotD180 <- results[[1]]  # it is a list
        RotD180 <- matrix(unlist(RotD180), nrow = 180) # transform it into matrix
        GMRotD <- results[[2]]  # it is a list
        GMRotD <- matrix(unlist(GMRotD), nrow = 90) # transform it into matrix
        Num_points <- results[[3]]  # it is a list
        Num_points <- c(unlist(Num_points)) # transform it into matrix
        rd_alevel <- results[[4]]

        ##################### Applying GMRotI_cal function to get PSA for different periods for only one desired rotated angle
        ########################## plug in the parameters ########################
        ##########################################################################
        Result <- GMRotI50_cal(GMRotD, period_t, tmax = tmax4penalty_in, tmin = tmin4penalty_in)
        GMRotI50 <- Result[[1]]
        GMRotI50_ang <- Result[[2]]
        penalty <- Result[[3]]
        # acceleration
        pga_gmrotI50 <- pga_gmrot[GMRotI50_ang]
        # velocity
        pgv_gmrotI50 <- pgv_gmrot[GMRotI50_ang]
        # displacement
        pgd_gmrotI50 <- pgd_gmrot[GMRotI50_ang]

        ##################### Applying RotD_cal function to get PSA for different periods for each rotated angle
        ########################## plug in the parameters ########################
        ##########################################################################
        RotD00 <- apply(RotD180, 2, min)
        RotD100 <- apply(RotD180, 2, max)
        RotD50 <- apply(RotD180, 2, median)
        RotD_t <- t(RotD180) # transpose of RotD matrix, it is used to get the max index
        RotD100_ang <- max.col(RotD_t)
        RotD00_ang <- max.col(-RotD_t)

        ############################# Process data End ###################################
        ###############################################################################
        if(combine_index == 0){
          pga_combine[i] <- pga_rot00
          pgv_combine[i] <- pgv_rot00
          pgd_combine[i] <- pgd_rot00
          rotd_combine[i,] <- RotD00
        } else if(combine_index == 50){
          pga_combine[i] <- pga_rot50
          pgv_combine[i] <- pgv_rot50
          pgd_combine[i] <- pgd_rot50
          rotd_combine[i,] <- RotD50
        } else{
          pga_combine[i] <- pga_rot100
          pgv_combine[i] <- pgv_rot100
          pgd_combine[i] <- pgd_rot100
          rotd_combine[i,] <- RotD100
        }

        #################### Output calculation results ####################
        ## check if required packages are installed
        PSA_1 <- PSA1[2,]
        Absolute_acc_1 <- PSA1[1,] # absolute acceleration of oscillator
        PSA_2 <- PSA2[2,]
        Absolute_acc_2 <- PSA2[1,]
        npts1 <- length(data1)
        npts2 <- length(data2)
        n_a_subset <- length(a_subset[1,])
        n_v_subset <- length(v_subset[1,])
        n_d_subset <- length(d_subset[1,])
        Period <- c("T0.010S","T0.020S","T0.022S","T0.025S","T0.029S","T0.030S",
                    "T0.032S","T0.035S","T0.036S","T0.040S","T0.042S","T0.044S",
                    "T0.045S","T0.046S","T0.048S","T0.050S","T0.055S","T0.060S",
                    "T0.065S","T0.067S","T0.070S","T0.075S","T0.080S","T0.085S",
                    "T0.090S","T0.095S","T0.100S","T0.110S","T0.120S","T0.130S",
                    "T0.133S","T0.140S","T0.150S","T0.160S","T0.170S","T0.180S",
                    "T0.190S","T0.200S","T0.220S","T0.240S","T0.250S","T0.260S",
                    "T0.280S","T0.290S","T0.300S","T0.320S","T0.340S","T0.350S",
                    "T0.360S","T0.380S","T0.400S","T0.420S","T0.440S","T0.450S",
                    "T0.460S","T0.480S","T0.500S","T0.550S","T0.600S","T0.650S",
                    "T0.667S","T0.700S","T0.750S","T0.800S","T0.850S","T0.900S",
                    "T0.950S","T1.000S","T1.100S","T1.200S","T1.300S","T1.400S",
                    "T1.500S","T1.600S","T1.700S","T1.800S","T1.900S","T2.000S",
                    "T2.200S","T2.400S","T2.500S","T2.600S","T2.800S","T3.000S",
                    "T3.200S","T3.400S","T3.500S","T3.600S","T3.800S","T4.000S",
                    "T4.200S","T4.400S","T4.600S","T4.800S","T5.000S","T5.500S",
                    "T6.000S","T6.500S","T7.000S","T7.500S","T8.000S","T8.500S",
                    "T9.000S","T9.500S","T10.000S","T11.000S","T12.000S","T13.000S",
                    "T14.000S","T15.000S","T20.000S")

        # compute real azimuth and set angle in [0,360]
        RotD00_ang <- RotD00_ang + ang1
        RotD100_ang <- RotD100_ang + ang1
        RotD00_ang[which(RotD00_ang > 360)] <- RotD00_ang[which(RotD00_ang > 360)] - 360
        RotD00_ang[which(RotD00_ang < 0)] <- RotD00_ang[which(RotD00_ang < 0)] + 360
        RotD100_ang[which(RotD100_ang > 360)] <- RotD100_ang[which(RotD100_ang > 360)] - 360
        RotD100_ang[which(RotD100_ang < 0)] <- RotD100_ang[which(RotD100_ang < 0)] + 360

        frame.data1 <- data.frame(Period, PSA_1, Absolute_acc_1, PSA_2, Absolute_acc_2,
                                  PSA_gm_ar, PSA_larger_ar, GMRotI50, RotD00, RotD50,
                                  RotD100, RotD00_ang, RotD100_ang, Num_points)
        sub_frame.data1 <- t(frame.data1[,-1])
        frame.data1 <- rbind(t(frame.data1[1]), sub_frame.data1)
        Name <- c("npts1","npts2","tmax4penalty","tmin4penalty","damping", "GMRotI50angle",
                  "PGA_GMRotI50", "PGV_GMRotI50",
                  "PGD_GMRotI50", "PGA_GMRot50", "PGV_GMRot50", "PGA_GMRot100angle",
                  "PGA_GMRot100", "PGV_GMRot100angle", "PGV_GMRot100", "PGA_1", "PGA_2",
                  "n_a_subset","PGV_1", "PGV_2", "n_v_subset","PGD_1", "PGD_2", "n_d_subset",
                  "PGA_Rot00", "PGA_Rot50", "PGA_Rot100", "PGA_Rot00angle", "PGA_Rot100angle",
                  "PGV_Rot00", "PGV_Rot50", "PGV_Rot100", "PGV_Rot00angle", "PGV_Rot100angle",
                  "PGD_Rot00", "PGD_Rot50", "PGD_Rot100", "PGD_Rot00angle", "PGD_Rot100angle")

        # compute real azimuth and set angel in [0,360]
        GMRotI50_ang <- GMRotI50_ang + ang1
        if(GMRotI50_ang > 360) {GMRotI50_ang <- GMRotI50_ang - 360}
        if(GMRotI50_ang < 0) {GMRotI50_ang <- GMRotI50_ang + 360}
        pga_gmrot100ang <- pga_gmrot100ang + ang1
        if(pga_gmrot100ang > 360) {pga_gmrot100ang <- pga_gmrot100ang - 360}
        if(pga_gmrot100ang < 0) {pga_gmrot100ang <- pga_gmrot100ang + 360}
        pgv_gmrot100ang <- pgv_gmrot100ang + ang1
        if(pgv_gmrot100ang > 360) {pgv_gmrot100ang <- pgv_gmrot100ang - 360}
        if(pgv_gmrot100ang < 0) {pgv_gmrot100ang <- pgv_gmrot100ang + 360}
        pga_rot00ang <- pga_rot00ang + ang1
        if(pga_rot00ang > 360) {pga_rot00ang <- pga_rot00ang - 360}
        if(pga_rot00ang < 0) {pga_rot00ang <- pga_rot00ang + 360}
        pga_rot100ang <- pga_rot100ang + ang1
        if(pga_rot100ang > 360) {pga_rot100ang <- pga_rot100ang - 360}
        if(pga_rot100ang < 0) {pga_rot100ang <- pga_rot100ang + 360}
        pgv_rot00ang <- pgv_rot00ang + ang1
        if(pgv_rot00ang > 360) {pgv_rot00ang <- pgv_rot00ang - 360}
        if(pgv_rot00ang < 0) {pgv_rot00ang <- pgv_rot00ang + 360}
        pgv_rot100ang <- pgv_rot100ang + ang1
        if(pgv_rot100ang > 360) {pgv_rot100ang <- pgv_rot100ang - 360}
        if(pgv_rot100ang < 0) {pgv_rot100ang <- pgv_rot100ang + 360}
        pgd_rot00ang <- pgd_rot00ang + ang1
        if(pgd_rot00ang > 360) {pgd_rot00ang <- pgd_rot00ang - 360}
        if(pgd_rot00ang < 0) {pgd_rot00ang <- pgd_rot00ang + 360}
        pgd_rot100ang <- pgd_rot100ang + ang1
        if(pgd_rot100ang > 360) {pgd_rot100ang <- pgd_rot100ang - 360}
        if(pgd_rot100ang < 0) {pgd_rot100ang <- pgd_rot100ang + 360}

        Value <- c(npts1, npts2, tmax4penalty_in, tmin4penalty_in, damping,
                   GMRotI50_ang, pga_gmrotI50, pgv_gmrotI50, pgd_gmrotI50, pga_gmrot50,
                   pgv_gmrot50, pga_gmrot100ang, pga_gmrot100, pgv_gmrot100ang, pgv_gmrot100,
                   pga_1_ar, pga_2_ar,n_a_subset, pgv_1_ar, pgv_2_ar, n_v_subset, pgd_1_ar,
                   pgd_2_ar, n_d_subset, pga_rot00, pga_rot50, pga_rot100, pga_rot00ang,
                   pga_rot100ang, pgv_rot00, pgv_rot50, pgv_rot100, pgv_rot00ang, pgv_rot100ang,
                   pgd_rot00, pgd_rot50, pgd_rot100, pgd_rot00ang, pgd_rot100ang)
        frame.data2 <- data.frame(Name, Value)
        filename1 <- paste(station_name, "_dep.csv", sep = "")
        filename2 <- paste(station_name, "_indep.csv", sep = "")
        setwd(outputdatadir)
        write.table(frame.data1, file = filename1, col.names = FALSE, sep = ',')
        write.csv(frame.data2, file = filename2, row.names = FALSE)

        #################### Plot curves  ########################################
        setwd(outputplotdir)
        png(filename = paste(station_name, "_PSA.png", sep = ""), width = 1280, height = 768)
        par(oma=c(2,2,2,2))
        par(mar=c(5,5,4,2) + 0.1)
        plot(c(0.01,20), c(0,max(max(PSA1,PSA2))), log = "x", type = "n", xlab = "Period(s)", ylab = "PSA(g)", main = "PSA versus Period", cex.main=2, cex.lab=2, cex.axis=2)
        lines(period_t, PSA_1, lty = 2, lwd = 2, col = "red")
        lines(period_t, PSA_2, lty = 3, lwd = 2, col = "blue")
        legend(5,max(max(PSA1,PSA2)), c("PSA_H1","PSA_H2"), lty=c(2,3), lwd=c(2,2),col=c("red","blue"), cex = 2)
        dev.off()
        png(filename = paste(station_name, "_GMRotI50.png", sep = ""), width = 1280, height = 768)
        par(oma=c(2,2,2,2))
        par(mar=c(5,5,4,2) + 0.1)
        plot(period_t, GMRotI50, log = "x", type = "l", lwd = 2, xlab = "Period(s)", ylab = "GMRotI50(g)", main = "GMRotI50 versus Period", cex.lab=2, cex.axis=2, cex.main=2)
        dev.off()
        png(filename = paste(station_name, "_RotD.png", sep = ""), width = 1280, height = 768)
        par(oma=c(2,2,2,2))
        par(mar=c(5,5,4,2) + 0.1)
        plot(c(0.01,20), c(0,max(RotD100)), log = "x", type = "n", xlab = "Period(s)", ylab = "RotDxx(g)", main = "RotDxx versus Period", cex.lab=2, cex.axis=2, cex.main=2)
        lines(period_t, RotD00,  lty = 2, lwd = 2, col = 2)
        lines(period_t, RotD50,  lty = 3, lwd = 2, col = 3)
        lines(period_t, RotD100, lty = 4, lwd = 2, col = 4)
        legend(5, max(RotD100), c("RotD00","RotD50","RotD100"), lty = c(2,3,4),  lwd = c(2,2,2), col = c(2,3,4), cex = 2)
        dev.off()

      }

    } else{
      print("You should put the couple groups in the Inputdata folder, we detect you have odd number of time series")
      stop()
    }

    ###################### Create a combinated spreadsheet ########################
    ###############################################################################
    name_com <- c('Station_Name')
    frame.data_com1 <- data.frame(Station)
    frame.data_com1 <- rbind(name_com, frame.data_com1)
    name_com <- c(paste('PGA_Rot', combine_index, sep = ''), paste('PGV_Rot', combine_index, sep = ''), paste('PGD_Rot', combine_index, sep = ''), Period)
    value_com <- cbind(pga_combine, pgv_combine, pgd_combine, rotd_combine)
    frame.data_com2 <- t(data.frame(name_com, t(value_com)))
    frame.data_com <- cbind(frame.data_com1, frame.data_com2)
    filename <- paste('summary', ".csv", sep = "")
    setwd(outputdatadir)
    write.table(frame.data_com, file = filename, col.names = FALSE, row.names = FALSE, sep = ',')


  }


  ############### process NGA-East/Sub formated data #################################
  else if(datatype == 'nga'){
    setwd(inputpath)
    temp = list.files(pattern = "*.AT2")
    if(!length(temp)%%2){            # judge if it is even number of groups
      ## define combined spreadsheet variables ##
      Record_Sequence_Number <- seq(1, length(temp)/2)
      pga_combine <- seq(1, length(temp)/2)*0
      pgv_combine <- seq(1, length(temp)/2)*0
      pgd_combine <- seq(1, length(temp)/2)*0
      rotd_combine <- matrix(seq(1, (length(temp)/2)*length(period_t)) ,nrow = length(temp)/2, ncol = length(period_t))*0

      for(i in 1:(length(temp)/2)){
        setwd(inputpath)
        name1 <- basename(temp[2*i-1])
        temp_name1 <- strsplit(name1, split = "[_]")
        name2 <- basename(temp[2*i])
        temp_name2 <- strsplit(name2, split = "[_]")
        if(temp_name1[[1]][1] != temp_name2[[1]][1]){      # judge if there is a coupled group for name1
          print(c("You should put the coupled-group in the Inputdata folder!", name1, "is not coupled. Quit!"))
          stop()
        }

        assign(temp[2*i-1], head_dt1 <- scan(temp[2*i-1], sep = ",", skip = 3, nlines = 1, what = character()))
        dt1 <- as.numeric(gsub("[^0-9\\.]","", head_dt1[2]))
        assign(temp[2*i], head_dt2 <- scan(temp[2*i], sep = ",", skip = 3, nlines = 1, what = character()))
        dt2 <- as.numeric(gsub("[^0-9\\.]","", head_dt2[2]))
        if(dt1 != dt2){
          print("Samples per second are different! Quit!")
          stop()
        }

        m <- gregexpr('[0-9]+',temp_name1[[1]][1])
        number <- regmatches(temp_name1[[1]][1],m)
        Record_Sequence_Number[i] <- as.numeric(number[[1]][1])
        assign(temp[2*i-1], data1 <- scan(temp[2*i-1], sep = "", skip = 4))
        assign(temp[2*i], data2 <- scan(temp[2*i], sep = "", skip = 4))

        ############################# Process data Begin ###################################
        ######### Sinc interpolate data ###############
        if(Interpolation_factor == 'auto'){
          # interp_factor <- ceiling(log(dt1/0.001)/log(2))
          # 1) the interpolation factor is based on f_nyq first -- we only care T_nyq,
          # because T less than T_nyq is not affected by interpolation,
          # therefore, dt/dt' = dt/(T_nyq/10) = dt/(dt/5) = 5, -> IF = 8
          # 2) if dt1 is less than the above calcualtion, then we use smaller interpolation
          interp_factor <- min(ceiling(log(dt1/0.001)/log(2)), 3)
          # if dt1 even smaller than the diresed time step (not common), we do not interpolate
          interp_factor <- max(interp_factor, 0)
          interp_factor <- 2^interp_factor
        }else if(log(Interpolation_factor)%%log(2) != 0){
          print("The Interpolation factor is not a number of power of 2, please try correct it!")
          stop()
        }else{
          interp_factor <- Interpolation_factor
        }

        dt1 <- dt1/interp_factor # compute new time step
        dt2 <- dt2/interp_factor
        data1 <- Interpft(data1, interp_factor) # new data after Sinc interpolation
        data2 <- Interpft(data2, interp_factor)
        number <- min(length(data1), length(data2)) # if two data sizes diff, we set both size as the smaller
        data1 <- data1[1:number]
        data2 <- data2[1:number]

        ############ Output PGA, PGV, PGD ###############
        ########## unit of original accel is g
        pga_1_ar <- max(abs(data1)) # unit, g
        pga_2_ar <- max(abs(data2)) # unit, g
        vel_1 <- cumtrapz(dt1*seq(1,length(data1)),data1) * 981 # unit, cm/s
        vel_2 <- cumtrapz(dt2*seq(1,length(data2)),data2) * 981 # unit, cm/s
        pgv_1_ar <- max(abs(vel_1))
        pgv_2_ar <- max(abs(vel_2))
        disp_1 <- cumtrapz(dt1*seq(1,length(vel_1)),vel_1) # unit, cm
        disp_2 <- cumtrapz(dt2*seq(1,length(vel_2)),vel_2) # unit, cm
        pgd_1_ar <- max(abs(disp_1))
        pgd_2_ar <- max(abs(disp_2))

        ############## Rotation values calculation ##############
        ######## subset slection #######
        length_min <- min(length(data1), length(data2))
        a_subset <- subset_select(data1, data2, fraction, length_min, dt1, 1)
        v_subset <- subset_select(vel_1, vel_2, fraction, length_min, dt1, 1)
        d_subset <- subset_select(disp_1, disp_2, fraction, length_min, dt1, 1)
        # define vectors
        pga_rot <- seq(1, 180)
        pga_gmrot <- seq(1,90)
        pgv_rot <- seq(1, 180)
        pgv_gmrot <- seq(1,90)
        pgd_rot <- seq(1, 180)
        pgd_gmrot <- seq(1,90)
        # calculation
        for(theta in seq(1,90)){
          # acceleration
          Rot_a1 <- a_subset[1,]*cos(theta/180*pi) + a_subset[2,]*sin(theta/180*pi)
          Rot_a2 <- -a_subset[1,]*sin(theta/180*pi) + a_subset[2,]*cos(theta/180*pi)
          pga_rot[theta] <- max(abs(Rot_a1))
          pga_rot[theta+90] <- max(abs(Rot_a2))
          pga_gmrot[theta] <- sqrt(max(abs(Rot_a1))*max(abs(Rot_a2)))
          # velocity
          Rot_v1 <- v_subset[1,]*cos(theta/180*pi) + v_subset[2,]*sin(theta/180*pi)
          Rot_v2 <- -v_subset[1,]*sin(theta/180*pi) + v_subset[2,]*cos(theta/180*pi)
          pgv_rot[theta] <- max(abs(Rot_v1))
          pgv_rot[theta+90] <- max(abs(Rot_v2))
          pgv_gmrot[theta] <- sqrt(max(abs(Rot_v1))*max(abs(Rot_v2)))
          # displacement
          Rot_d1 <- d_subset[1,]*cos(theta/180*pi) + d_subset[2,]*sin(theta/180*pi)
          Rot_d2 <- -d_subset[1,]*sin(theta/180*pi) + d_subset[2,]*cos(theta/180*pi)
          pgd_rot[theta] <- max(abs(Rot_d1))
          pgd_rot[theta+90] <- max(abs(Rot_d2))
          pgd_gmrot[theta] <- sqrt(max(abs(Rot_d1))*max(abs(Rot_d2)))
        }
        ## acceleration
        pga_rot00 <- min(pga_rot)
        pga_rot50 <- median(pga_rot)
        pga_rot100 <- max(pga_rot)
        pga_rot00ang <- which.min(pga_rot)
        pga_rot100ang <- which.max(pga_rot)
        pga_gmrot50 <- median(pga_gmrot)
        pga_gmrot100 <- max(pga_gmrot)
        pga_gmrot100ang <- which.max(pga_gmrot)
        ## velocity
        pgv_rot00 <- min(pgv_rot)
        pgv_rot50 <- median(pgv_rot)
        pgv_rot100 <- max(pgv_rot)
        pgv_rot00ang <- which.min(pgv_rot)
        pgv_rot100ang <- which.max(pgv_rot)
        pgv_gmrot50 <- median(pgv_gmrot)
        pgv_gmrot100 <- max(pgv_gmrot)
        pgv_gmrot100ang <- which.max(pgv_gmrot)
        ## displacement
        pgd_rot00 <- min(pgd_rot)
        pgd_rot50 <- median(pgd_rot)
        pgd_rot100 <- max(pgd_rot)
        pgd_rot00ang <- which.min(pgd_rot)
        pgd_rot100ang <- which.max(pgd_rot)
        pgd_gmrot50 <- median(pgd_gmrot)
        pgd_gmrot100 <- max(pgd_gmrot)
        pgd_gmrot100ang <- which.max(pgd_gmrot)

        ################# Applying PS_cal function to get PSA,PSV #############
        ########################## plug in the parameters ########################
        ##########################################
        # plug in the parameters for two components
        # using C++ to compute to speed up calculation
        # PSA calculation is from textbook of "Dynamics of Structure" 3rd edition by Chopra on page 168-169.
        ##### compute pseudo spectral acceleration
        PSA1 <- PS_cal_cpp(data1, period_t, damping, dt1, 1) # first row is absolute acceleration of oscillator, second row is psa
        PSA2 <- PS_cal_cpp(data2, period_t, damping, dt1, 1) # first row is absolute acceleration of oscillator, second row is psa
        PSA_gm_ar <- sqrt(PSA1[2,]*PSA2[2,])
        PSA_larger_ar <- seq(1, length(period_t))
        for(j in seq(1, length(PSA_larger_ar))){
          PSA_larger_ar[j] <- max(PSA1[2,j], PSA2[2,j])
        }

        ############## Applying GMRotD_cal function to get PSA for different periods for each rotated angle ####
        ########################## plug in the parameters ########################
        ##########################################################################
        results <- GM_RotD_cal(data1, data2, period_t, damping, dt1, fraction)  # the result is list, 1st is RotD180, 2nd is GMRotD
        RotD180 <- results[[1]]  # it is a list
        RotD180 <- matrix(unlist(RotD180), nrow = 180) # transform it into matrix
        GMRotD <- results[[2]]  # it is a list
        GMRotD <- matrix(unlist(GMRotD), nrow = 90) # transform it into matrix
        Num_points <- results[[3]]  # it is a list
        Num_points <- c(unlist(Num_points)) # transform it into matrix
        rd_alevel <- results[[4]]

        ##################### Applying GMRotI_cal function to get PSA for different periods for only one desired rotated angle
        ########################## plug in the parameters ########################
        ##########################################################################
        Result <- GMRotI50_cal(GMRotD, period_t, tmax = tmax4penalty_in, tmin = tmin4penalty_in)
        GMRotI50 <- Result[[1]]
        GMRotI50_ang <- Result[[2]]
        penalty <- Result[[3]]
        # acceleration
        pga_gmrotI50 <- pga_gmrot[GMRotI50_ang]
        # velocity
        pgv_gmrotI50 <- pgv_gmrot[GMRotI50_ang]
        # displacement
        pgd_gmrotI50 <- pgd_gmrot[GMRotI50_ang]

        ##################### Applying RotD_cal function to get PSA for different periods for each rotated angle
        ########################## plug in the parameters ########################
        ##########################################################################
        RotD00 <- apply(RotD180, 2, min)
        RotD100 <- apply(RotD180, 2, max)
        RotD50 <- apply(RotD180, 2, median)
        RotD_t <- t(RotD180) # transpose of RotD matrix, it is used to get the max index
        RotD100_ang <- max.col(RotD_t)
        RotD00_ang <- max.col(-RotD_t)

        ############################# Process data End ###################################
        ###############################################################################
        if(combine_index == 0){
          pga_combine[i] <- pga_rot00
          pgv_combine[i] <- pgv_rot00
          pgd_combine[i] <- pgd_rot00
          rotd_combine[i,] <- RotD00
        } else if(combine_index == 50){
          pga_combine[i] <- pga_rot50
          pgv_combine[i] <- pgv_rot50
          pgd_combine[i] <- pgd_rot50
          rotd_combine[i,] <- RotD50
        } else{
          pga_combine[i] <- pga_rot100
          pgv_combine[i] <- pgv_rot100
          pgd_combine[i] <- pgd_rot100
          rotd_combine[i,] <- RotD100
        }

        #################### Output calculation results ####################
        ## check if required packages are installed
        PSA_1 <- PSA1[2,]
        Absolute_acc_1 <- PSA1[1,] # absolute acceleration of oscillator
        PSA_2 <- PSA2[2,]
        Absolute_acc_2 <- PSA2[1,]
        npts1 <- length(data1)
        npts2 <- length(data2)
        n_a_subset <- length(a_subset[1,])
        n_v_subset <- length(v_subset[1,])
        n_d_subset <- length(d_subset[1,])
        Period <- c("T0.010S","T0.020S","T0.022S","T0.025S","T0.029S","T0.030S",
                    "T0.032S","T0.035S","T0.036S","T0.040S","T0.042S","T0.044S",
                    "T0.045S","T0.046S","T0.048S","T0.050S","T0.055S","T0.060S",
                    "T0.065S","T0.067S","T0.070S","T0.075S","T0.080S","T0.085S",
                    "T0.090S","T0.095S","T0.100S","T0.110S","T0.120S","T0.130S",
                    "T0.133S","T0.140S","T0.150S","T0.160S","T0.170S","T0.180S",
                    "T0.190S","T0.200S","T0.220S","T0.240S","T0.250S","T0.260S",
                    "T0.280S","T0.290S","T0.300S","T0.320S","T0.340S","T0.350S",
                    "T0.360S","T0.380S","T0.400S","T0.420S","T0.440S","T0.450S",
                    "T0.460S","T0.480S","T0.500S","T0.550S","T0.600S","T0.650S",
                    "T0.667S","T0.700S","T0.750S","T0.800S","T0.850S","T0.900S",
                    "T0.950S","T1.000S","T1.100S","T1.200S","T1.300S","T1.400S",
                    "T1.500S","T1.600S","T1.700S","T1.800S","T1.900S","T2.000S",
                    "T2.200S","T2.400S","T2.500S","T2.600S","T2.800S","T3.000S",
                    "T3.200S","T3.400S","T3.500S","T3.600S","T3.800S","T4.000S",
                    "T4.200S","T4.400S","T4.600S","T4.800S","T5.000S","T5.500S",
                    "T6.000S","T6.500S","T7.000S","T7.500S","T8.000S","T8.500S",
                    "T9.000S","T9.500S","T10.000S","T11.000S","T12.000S","T13.000S",
                    "T14.000S","T15.000S","T20.000S")

        # compute real azimuth and set angle in [0,360]
        RotD00_ang <- RotD00_ang + ang1
        RotD100_ang <- RotD100_ang + ang1
        RotD00_ang[which(RotD00_ang > 360)] <- RotD00_ang[which(RotD00_ang > 360)] - 360
        RotD00_ang[which(RotD00_ang < 0)] <- RotD00_ang[which(RotD00_ang < 0)] + 360
        RotD100_ang[which(RotD100_ang > 360)] <- RotD100_ang[which(RotD100_ang > 360)] - 360
        RotD100_ang[which(RotD100_ang < 0)] <- RotD100_ang[which(RotD100_ang < 0)] + 360

        frame.data1 <- data.frame(Period, PSA_1, Absolute_acc_1, PSA_2, Absolute_acc_2,
                                  PSA_gm_ar, PSA_larger_ar, GMRotI50, RotD00, RotD50,
                                  RotD100, RotD00_ang, RotD100_ang, Num_points)
        sub_frame.data1 <- t(frame.data1[,-1])
        frame.data1 <- rbind(t(frame.data1[1]), sub_frame.data1)
        Name <- c("npts1","npts2","tmax4penalty","tmin4penalty","damping", "GMRotI50angle",
                  "PGA_GMRotI50", "PGV_GMRotI50",
                  "PGD_GMRotI50", "PGA_GMRot50", "PGV_GMRot50", "PGA_GMRot100angle",
                  "PGA_GMRot100", "PGV_GMRot100angle", "PGV_GMRot100", "PGA_1", "PGA_2",
                  "n_a_subset","PGV_1", "PGV_2", "n_v_subset","PGD_1", "PGD_2", "n_d_subset",
                  "PGA_Rot00", "PGA_Rot50", "PGA_Rot100", "PGA_Rot00angle", "PGA_Rot100angle",
                  "PGV_Rot00", "PGV_Rot50", "PGV_Rot100", "PGV_Rot00angle", "PGV_Rot100angle",
                  "PGD_Rot00", "PGD_Rot50", "PGD_Rot100", "PGD_Rot00angle", "PGD_Rot100angle")

        # compute real azimuth and set angel in [0,360]
        GMRotI50_ang <- GMRotI50_ang + ang1
        if(GMRotI50_ang > 360) {GMRotI50_ang <- GMRotI50_ang - 360}
        if(GMRotI50_ang < 0) {GMRotI50_ang <- GMRotI50_ang + 360}
        pga_gmrot100ang <- pga_gmrot100ang + ang1
        if(pga_gmrot100ang > 360) {pga_gmrot100ang <- pga_gmrot100ang - 360}
        if(pga_gmrot100ang < 0) {pga_gmrot100ang <- pga_gmrot100ang + 360}
        pgv_gmrot100ang <- pgv_gmrot100ang + ang1
        if(pgv_gmrot100ang > 360) {pgv_gmrot100ang <- pgv_gmrot100ang - 360}
        if(pgv_gmrot100ang < 0) {pgv_gmrot100ang <- pgv_gmrot100ang + 360}
        pga_rot00ang <- pga_rot00ang + ang1
        if(pga_rot00ang > 360) {pga_rot00ang <- pga_rot00ang - 360}
        if(pga_rot00ang < 0) {pga_rot00ang <- pga_rot00ang + 360}
        pga_rot100ang <- pga_rot100ang + ang1
        if(pga_rot100ang > 360) {pga_rot100ang <- pga_rot100ang - 360}
        if(pga_rot100ang < 0) {pga_rot100ang <- pga_rot100ang + 360}
        pgv_rot00ang <- pgv_rot00ang + ang1
        if(pgv_rot00ang > 360) {pgv_rot00ang <- pgv_rot00ang - 360}
        if(pgv_rot00ang < 0) {pgv_rot00ang <- pgv_rot00ang + 360}
        pgv_rot100ang <- pgv_rot100ang + ang1
        if(pgv_rot100ang > 360) {pgv_rot100ang <- pgv_rot100ang - 360}
        if(pgv_rot100ang < 0) {pgv_rot100ang <- pgv_rot100ang + 360}
        pgd_rot00ang <- pgd_rot00ang + ang1
        if(pgd_rot00ang > 360) {pgd_rot00ang <- pgd_rot00ang - 360}
        if(pgd_rot00ang < 0) {pgd_rot00ang <- pgd_rot00ang + 360}
        pgd_rot100ang <- pgd_rot100ang + ang1
        if(pgd_rot100ang > 360) {pgd_rot100ang <- pgd_rot100ang - 360}
        if(pgd_rot100ang < 0) {pgd_rot100ang <- pgd_rot100ang + 360}

        Value <- c(npts1, npts2, tmax4penalty_in, tmin4penalty_in, damping,
                   GMRotI50_ang, pga_gmrotI50, pgv_gmrotI50, pgd_gmrotI50, pga_gmrot50,
                   pgv_gmrot50, pga_gmrot100ang, pga_gmrot100, pgv_gmrot100ang, pgv_gmrot100,
                   pga_1_ar, pga_2_ar,n_a_subset, pgv_1_ar, pgv_2_ar, n_v_subset, pgd_1_ar,
                   pgd_2_ar, n_d_subset, pga_rot00, pga_rot50, pga_rot100, pga_rot00ang,
                   pga_rot100ang, pgv_rot00, pgv_rot50, pgv_rot100, pgv_rot00ang, pgv_rot100ang,
                   pgd_rot00, pgd_rot50, pgd_rot100, pgd_rot00ang, pgd_rot100ang)
        frame.data2 <- data.frame(Name, Value)
        filename1 <- paste(temp_name1[[1]][1], "_dep.csv", sep = "")
        filename2 <- paste(temp_name1[[1]][1], "_indep.csv", sep = "")
        setwd(outputdatadir)
        write.table(frame.data1, file = filename1, col.names = FALSE, sep = ',')
        write.csv(frame.data2, file = filename2, row.names = FALSE)

        #################### Plot curves  ########################################
        setwd(outputplotdir)
        png(filename = paste(temp_name1[[1]][1], "_PSA.png", sep = ""), width = 1280, height = 768)
        par(oma=c(2,2,2,2))
        par(mar=c(5,5,4,2) + 0.1)
        plot(c(0.01,20), c(0,max(max(PSA1,PSA2))), log = "x", type = "n", xlab = "Period(s)", ylab = "PSA(g)", main = "PSA versus Period", cex.main=2, cex.lab=2, cex.axis=2)
        lines(period_t, PSA_1, lty = 2, lwd = 2, col = "red")
        lines(period_t, PSA_2, lty = 3, lwd = 2, col = "blue")
        legend(5,max(max(PSA1,PSA2)), c("PSA_H1","PSA_H2"), lty=c(2,3), lwd=c(2,2),col=c("red","blue"), cex = 2)
        dev.off()
        png(filename = paste(temp_name1[[1]][1], "_GMRotI50.png", sep = ""), width = 1280, height = 768)
        par(oma=c(2,2,2,2))
        par(mar=c(5,5,4,2) + 0.1)
        plot(period_t, GMRotI50, log = "x", type = "l", lwd = 2, xlab = "Period(s)", ylab = "GMRotI50(g)", main = "GMRotI50 versus Period", cex.lab=2, cex.axis=2, cex.main=2)
        dev.off()
        png(filename = paste(temp_name1[[1]][1], "_RotD.png", sep = ""), width = 1280, height = 768)
        par(oma=c(2,2,2,2))
        par(mar=c(5,5,4,2) + 0.1)
        plot(c(0.01,20), c(0,max(RotD100)), log = "x", type = "n", xlab = "Period(s)", ylab = "RotDxx(g)", main = "RotDxx versus Period", cex.lab=2, cex.axis=2, cex.main=2)
        lines(period_t, RotD00,  lty = 2, lwd = 2, col = 2)
        lines(period_t, RotD50,  lty = 3, lwd = 2, col = 3)
        lines(period_t, RotD100, lty = 4, lwd = 2, col = 4)
        legend(5, max(RotD100), c("RotD00","RotD50","RotD100"), lty = c(2,3,4),  lwd = c(2,2,2), col = c(2,3,4), cex = 2)
        dev.off()

      }

    } else{
      print("You should put the couple groups in the Inputdata folder, we detect you have odd number of time series")
      stop()
    }

    ###################### Create a combinated spreadsheet ########################
    ###############################################################################
    name_com <- c('Record_Sequence_Number')
    frame.data_com1 <- data.frame(Record_Sequence_Number)
    frame.data_com1 <- rbind(name_com, frame.data_com1)
    name_com <- c(paste('PGA_Rot', combine_index, sep = ''), paste('PGV_Rot', combine_index, sep = ''), paste('PGD_Rot', combine_index, sep = ''), Period)
    value_com <- cbind(pga_combine, pgv_combine, pgd_combine, rotd_combine)
    frame.data_com2 <- t(data.frame(name_com, t(value_com)))
    frame.data_com <- cbind(frame.data_com1, frame.data_com2)
    filename <- paste('summary', ".csv", sep = "")
    setwd(outputdatadir)
    write.table(frame.data_com, file = filename, col.names = FALSE, row.names = FALSE, sep = ',')


  }

  ############### process COSMOS formated data #################################
  else if(datatype == 'cosmos'){
    ## It only processes the corrected COSMOS data. The accelerations, velocity, and displacements
    ## of two horizontal and one vertical components of ground motions are all in the one file
    setwd(inputpath)
    temp = list.files()
    Station <- seq(1, length(temp))
    pga_combine <- seq(1, length(temp))*0
    pgv_combine <- seq(1, length(temp))*0
    pgd_combine <- seq(1, length(temp))*0
    rotd_combine <- matrix(seq(1, length(temp)*length(period_t)), nrow = length(temp), ncol = length(period_t))*0
    for(i in seq(1,length(temp))){
      setwd(inputpath)
      name <- basename(temp[i])
      temp_name <- strsplit(name, split = "_")
      station_name <- paste(temp_name[[1]][1], temp_name[[1]][2], sep = "_")
      Station[i] <- station_name
      assign(temp[i], head_dt1 <- scan(temp[i], sep = "", skip = 45, nlines = 1, what = character()))
      npts1 <- as.numeric(head_dt1[1])
      num_ind <- which(!is.na(as.numeric(head_dt1)))
      dt1 <- as.numeric(head_dt1[num_ind[2]])
      if(dt1 > 1 ){
        dt1 <- 1/dt1
      }
      len_h1 <- 45 + (ceiling(npts1/8)+1)*3 + 1
      assign(temp[i], head_dt2 <- scan(temp[i], sep = "", skip = len_h1 + 45, nlines = 1, what = character()))
      npts2 <- as.numeric(head_dt2[1])
      num_ind <- which(!is.na(as.numeric(head_dt2)))
      dt2 <- as.numeric(head_dt2[num_ind[2]])
      if(dt2 > 1 ){
        dt2 <- 1/dt2
      }
      if(dt1 != dt2){
          print("Samples per second are different! Quit!")
          stop()
        }
      line1 <- readLines(temp[i], n = ceiling(npts1/8)+46)
      data1 <- c()
      for(k in seq(47,ceiling(npts1/8)+46)){
        if(length(data1) == 0){
          data1 <- as.numeric(substring(line1[k], c(1,11,21,31,41,51,61,71),c(10,20,30,40,50,60,70,80)))
        } else{
          data1 <- c(data1,as.numeric(substring(line1[k], c(1,11,21,31,41,51,61,71),c(10,20,30,40,50,60,70,80))))
        }
      }
      data1 <- data1[1:npts1]

      line2 <- readLines(temp[i], n = len_h1 + ceiling(npts1/8)+46)
      data2 <- c()
      for(k in seq(len_h1 + 47, len_h1+ceiling(npts1/8)+46)){
        if(length(data2) == 0){
          data2 <- as.numeric(substring(line2[k], c(1,11,21,31,41,51,61,71),c(10,20,30,40,50,60,70,80)))
        } else{
          data2 <- c(data2,as.numeric(substring(line2[k], c(1,11,21,31,41,51,61,71),c(10,20,30,40,50,60,70,80))))
        }
      }
      data2 <- data2[1:npts2]

      ############################# Process data Begin ###################################
      ######### Sinc interpolate data ###############
      if(Interpolation_factor == 'auto'){
        # interp_factor <- ceiling(log(dt1/0.001)/log(2))

        # 1) the interpolation factor is based on f_nyq first -- we only care T_nyq,
        # because T less than T_nyq is not affected by interpolation,
        # therefore, dt/dt' = dt/(T_nyq/10) = dt/(dt/5) = 5, -> IF = 8
        # 2) if dt1 is less than the above calcualtion, then we use smaller interpolation
        interp_factor <- min(ceiling(log(dt1/0.001)/log(2)), 3)
        # if dt1 even smaller than the diresed time step (not common), we do not interpolate
        interp_factor <- max(interp_factor, 0)
        interp_factor <- 2^interp_factor
      }else if(log(Interpolation_factor)%%log(2) != 0){
        print("The Interpolation factor is not a number of power of 2, please try correct it!")
        stop()
      }else{
        interp_factor <- Interpolation_factor
      }

      dt1 <- dt1/interp_factor # compute new time step
      dt2 <- dt2/interp_factor
      data1 <- Interpft(data1, interp_factor) # new data after Sinc interpolation
      data2 <- Interpft(data2, interp_factor)
      number <- min(length(data1), length(data2)) # if two data sizes diff, we set both size as the smaller
      data1 <- data1[1:number]
      data2 <- data2[1:number]

      ############ Output PGA, PGV, PGD ###############
      ########## unit of original accel is cm/sec/sec
      ### change unit into g
      data1 <- data1/981 # unit, g
      data2 <- data2/981 # unit, g
      pga_1_ar <- max(abs(data1)) # unit, g
      pga_2_ar <- max(abs(data2)) # unit, g
      vel_1 <- cumtrapz(dt1*seq(1,length(data1)),data1) * 981 # unit, cm/s
      vel_2 <- cumtrapz(dt2*seq(1,length(data2)),data2) * 981 # unit, cm/s
      pgv_1_ar <- max(abs(vel_1))
      pgv_2_ar <- max(abs(vel_2))
      disp_1 <- cumtrapz(dt1*seq(1,length(vel_1)),vel_1) # unit, cm
      disp_2 <- cumtrapz(dt2*seq(1,length(vel_2)),vel_2) # unit, cm
      pgd_1_ar <- max(abs(disp_1))
      pgd_2_ar <- max(abs(disp_2))

      ############## Rotation values calculation ##############
      ######## subset slection #######
      length_min <- min(length(data1), length(data2))
      a_subset <- subset_select(data1, data2, fraction, length_min, dt1, 1)
      v_subset <- subset_select(vel_1, vel_2, fraction, length_min, dt1, 1)
      d_subset <- subset_select(disp_1, disp_2, fraction, length_min, dt1, 1)
      # define vectors
      pga_rot <- seq(1, 180)
      pga_gmrot <- seq(1,90)
      pgv_rot <- seq(1, 180)
      pgv_gmrot <- seq(1,90)
      pgd_rot <- seq(1, 180)
      pgd_gmrot <- seq(1,90)
      # calculation
      for(theta in seq(1,90)){
        # acceleration
        Rot_a1 <- a_subset[1,]*cos(theta/180*pi) + a_subset[2,]*sin(theta/180*pi)
        Rot_a2 <- -a_subset[1,]*sin(theta/180*pi) + a_subset[2,]*cos(theta/180*pi)
        pga_rot[theta] <- max(abs(Rot_a1))
        pga_rot[theta+90] <- max(abs(Rot_a2))
        pga_gmrot[theta] <- sqrt(max(abs(Rot_a1))*max(abs(Rot_a2)))
        # velocity
        Rot_v1 <- v_subset[1,]*cos(theta/180*pi) + v_subset[2,]*sin(theta/180*pi)
        Rot_v2 <- -v_subset[1,]*sin(theta/180*pi) + v_subset[2,]*cos(theta/180*pi)
        pgv_rot[theta] <- max(abs(Rot_v1))
        pgv_rot[theta+90] <- max(abs(Rot_v2))
        pgv_gmrot[theta] <- sqrt(max(abs(Rot_v1))*max(abs(Rot_v2)))
        # displacement
        Rot_d1 <- d_subset[1,]*cos(theta/180*pi) + d_subset[2,]*sin(theta/180*pi)
        Rot_d2 <- -d_subset[1,]*sin(theta/180*pi) + d_subset[2,]*cos(theta/180*pi)
        pgd_rot[theta] <- max(abs(Rot_d1))
        pgd_rot[theta+90] <- max(abs(Rot_d2))
        pgd_gmrot[theta] <- sqrt(max(abs(Rot_d1))*max(abs(Rot_d2)))
      }
      ## acceleration
      pga_rot00 <- min(pga_rot)
      pga_rot50 <- median(pga_rot)
      pga_rot100 <- max(pga_rot)
      pga_rot00ang <- which.min(pga_rot)
      pga_rot100ang <- which.max(pga_rot)
      pga_gmrot50 <- median(pga_gmrot)
      pga_gmrot100 <- max(pga_gmrot)
      pga_gmrot100ang <- which.max(pga_gmrot)
      ## velocity
      pgv_rot00 <- min(pgv_rot)
      pgv_rot50 <- median(pgv_rot)
      pgv_rot100 <- max(pgv_rot)
      pgv_rot00ang <- which.min(pgv_rot)
      pgv_rot100ang <- which.max(pgv_rot)
      pgv_gmrot50 <- median(pgv_gmrot)
      pgv_gmrot100 <- max(pgv_gmrot)
      pgv_gmrot100ang <- which.max(pgv_gmrot)
      ## displacement
      pgd_rot00 <- min(pgd_rot)
      pgd_rot50 <- median(pgd_rot)
      pgd_rot100 <- max(pgd_rot)
      pgd_rot00ang <- which.min(pgd_rot)
      pgd_rot100ang <- which.max(pgd_rot)
      pgd_gmrot50 <- median(pgd_gmrot)
      pgd_gmrot100 <- max(pgd_gmrot)
      pgd_gmrot100ang <- which.max(pgd_gmrot)

      ################# Applying PS_cal function to get PSA,PSV #############
      ########################## plug in the parameters ########################
      ##########################################
      # plug in the parameters for two components
      # using C++ to compute to speed up calculation
      # PSA calculation is from textbook of "Dynamics of Structure" 3rd edition by Chopra on page 168-169.
      ##### compute pseudo spectral acceleration
      PSA1 <- PS_cal_cpp(data1, period_t, damping, dt1, 1) # first row is absolute acceleration of oscillator, second row is psa
      PSA2 <- PS_cal_cpp(data2, period_t, damping, dt1, 1) # first row is absolute acceleration of oscillator, second row is psa
      PSA_gm_ar <- sqrt(PSA1[2,]*PSA2[2,])
      PSA_larger_ar <- seq(1, length(period_t))
      for(j in seq(1, length(PSA_larger_ar))){
        PSA_larger_ar[j] <- max(PSA1[2,j], PSA2[2,j])
      }

      ############## Applying GMRotD_cal function to get PSA for different periods for each rotated angle ####
      ########################## plug in the parameters ########################
      ##########################################################################
      results <- GM_RotD_cal(data1, data2, period_t, damping, dt1, fraction)  # the result is list, 1st is RotD180, 2nd is GMRotD
      RotD180 <- results[[1]]  # it is a list
      RotD180 <- matrix(unlist(RotD180), nrow = 180) # transform it into matrix
      GMRotD <- results[[2]]  # it is a list
      GMRotD <- matrix(unlist(GMRotD), nrow = 90) # transform it into matrix
      Num_points <- results[[3]]  # it is a list
      Num_points <- c(unlist(Num_points)) # transform it into matrix
      rd_alevel <- results[[4]]

      ##################### Applying GMRotI_cal function to get PSA for different periods for only one desired rotated angle
      ########################## plug in the parameters ########################
      ##########################################################################
      Result <- GMRotI50_cal(GMRotD, period_t, tmax = tmax4penalty_in, tmin = tmin4penalty_in)
      GMRotI50 <- Result[[1]]
      GMRotI50_ang <- Result[[2]]
      penalty <- Result[[3]]
      # acceleration
      pga_gmrotI50 <- pga_gmrot[GMRotI50_ang]
      # velocity
      pgv_gmrotI50 <- pgv_gmrot[GMRotI50_ang]
      # displacement
      pgd_gmrotI50 <- pgd_gmrot[GMRotI50_ang]

      ##################### Applying RotD_cal function to get PSA for different periods for each rotated angle
      ########################## plug in the parameters ########################
      ##########################################################################
      RotD00 <- apply(RotD180, 2, min)
      RotD100 <- apply(RotD180, 2, max)
      RotD50 <- apply(RotD180, 2, median)
      RotD_t <- t(RotD180) # transpose of RotD matrix, it is used to get the max index
      RotD100_ang <- max.col(RotD_t)
      RotD00_ang <- max.col(-RotD_t)

      ############################# Process data End ###################################
      ###############################################################################
      if(combine_index == 0){
        pga_combine[i] <- pga_rot00
        pgv_combine[i] <- pgv_rot00
        pgd_combine[i] <- pgd_rot00
        rotd_combine[i,] <- RotD00
      } else if(combine_index == 50){
        pga_combine[i] <- pga_rot50
        pgv_combine[i] <- pgv_rot50
        pgd_combine[i] <- pgd_rot50
        rotd_combine[i,] <- RotD50
      } else{
        pga_combine[i] <- pga_rot100
        pgv_combine[i] <- pgv_rot100
        pgd_combine[i] <- pgd_rot100
        rotd_combine[i,] <- RotD100
      }

      #################### Output calculation results ####################
      PSA_1 <- PSA1[2,]
      Absolute_acc_1 <- PSA1[1,] # absolute acceleration of oscillator
      PSA_2 <- PSA2[2,]
      Absolute_acc_2 <- PSA2[1,]
      npts1 <- length(data1)
      npts2 <- length(data2)
      n_a_subset <- length(a_subset[1,])
      n_v_subset <- length(v_subset[1,])
      n_d_subset <- length(d_subset[1,])
      Period <- c("T0.010S","T0.020S","T0.022S","T0.025S","T0.029S","T0.030S",
                  "T0.032S","T0.035S","T0.036S","T0.040S","T0.042S","T0.044S",
                  "T0.045S","T0.046S","T0.048S","T0.050S","T0.055S","T0.060S",
                  "T0.065S","T0.067S","T0.070S","T0.075S","T0.080S","T0.085S",
                  "T0.090S","T0.095S","T0.100S","T0.110S","T0.120S","T0.130S",
                  "T0.133S","T0.140S","T0.150S","T0.160S","T0.170S","T0.180S",
                  "T0.190S","T0.200S","T0.220S","T0.240S","T0.250S","T0.260S",
                  "T0.280S","T0.290S","T0.300S","T0.320S","T0.340S","T0.350S",
                  "T0.360S","T0.380S","T0.400S","T0.420S","T0.440S","T0.450S",
                  "T0.460S","T0.480S","T0.500S","T0.550S","T0.600S","T0.650S",
                  "T0.667S","T0.700S","T0.750S","T0.800S","T0.850S","T0.900S",
                  "T0.950S","T1.000S","T1.100S","T1.200S","T1.300S","T1.400S",
                  "T1.500S","T1.600S","T1.700S","T1.800S","T1.900S","T2.000S",
                  "T2.200S","T2.400S","T2.500S","T2.600S","T2.800S","T3.000S",
                  "T3.200S","T3.400S","T3.500S","T3.600S","T3.800S","T4.000S",
                  "T4.200S","T4.400S","T4.600S","T4.800S","T5.000S","T5.500S",
                  "T6.000S","T6.500S","T7.000S","T7.500S","T8.000S","T8.500S",
                  "T9.000S","T9.500S","T10.000S","T11.000S","T12.000S","T13.000S",
                  "T14.000S","T15.000S","T20.000S")

      # compute real azimuth and set angle in [0,360]
      RotD00_ang <- RotD00_ang + ang1
      RotD100_ang <- RotD100_ang + ang1
      RotD00_ang[which(RotD00_ang > 360)] <- RotD00_ang[which(RotD00_ang > 360)] - 360
      RotD00_ang[which(RotD00_ang < 0)] <- RotD00_ang[which(RotD00_ang < 0)] + 360
      RotD100_ang[which(RotD100_ang > 360)] <- RotD100_ang[which(RotD100_ang > 360)] - 360
      RotD100_ang[which(RotD100_ang < 0)] <- RotD100_ang[which(RotD100_ang < 0)] + 360

      frame.data1 <- data.frame(Period, PSA_1, Absolute_acc_1, PSA_2, Absolute_acc_2,
                                PSA_gm_ar, PSA_larger_ar, GMRotI50, RotD00, RotD50,
                                RotD100, RotD00_ang, RotD100_ang, Num_points)
      sub_frame.data1 <- t(frame.data1[,-1])
      frame.data1 <- rbind(t(frame.data1[1]), sub_frame.data1)
      Name <- c("npts1","npts2","tmax4penalty","tmin4penalty","damping", "GMRotI50angle",
                "PGA_GMRotI50", "PGV_GMRotI50",
                "PGD_GMRotI50", "PGA_GMRot50", "PGV_GMRot50", "PGA_GMRot100angle",
                "PGA_GMRot100", "PGV_GMRot100angle", "PGV_GMRot100", "PGA_1", "PGA_2",
                "n_a_subset","PGV_1", "PGV_2", "n_v_subset","PGD_1", "PGD_2", "n_d_subset",
                "PGA_Rot00", "PGA_Rot50", "PGA_Rot100", "PGA_Rot00angle", "PGA_Rot100angle",
                "PGV_Rot00", "PGV_Rot50", "PGV_Rot100", "PGV_Rot00angle", "PGV_Rot100angle",
                "PGD_Rot00", "PGD_Rot50", "PGD_Rot100", "PGD_Rot00angle", "PGD_Rot100angle")

      # compute real azimuth and set angel in [0,360]
      GMRotI50_ang <- GMRotI50_ang + ang1
      if(GMRotI50_ang > 360) {GMRotI50_ang <- GMRotI50_ang - 360}
      if(GMRotI50_ang < 0) {GMRotI50_ang <- GMRotI50_ang + 360}
      pga_gmrot100ang <- pga_gmrot100ang + ang1
      if(pga_gmrot100ang > 360) {pga_gmrot100ang <- pga_gmrot100ang - 360}
      if(pga_gmrot100ang < 0) {pga_gmrot100ang <- pga_gmrot100ang + 360}
      pgv_gmrot100ang <- pgv_gmrot100ang + ang1
      if(pgv_gmrot100ang > 360) {pgv_gmrot100ang <- pgv_gmrot100ang - 360}
      if(pgv_gmrot100ang < 0) {pgv_gmrot100ang <- pgv_gmrot100ang + 360}
      pga_rot00ang <- pga_rot00ang + ang1
      if(pga_rot00ang > 360) {pga_rot00ang <- pga_rot00ang - 360}
      if(pga_rot00ang < 0) {pga_rot00ang <- pga_rot00ang + 360}
      pga_rot100ang <- pga_rot100ang + ang1
      if(pga_rot100ang > 360) {pga_rot100ang <- pga_rot100ang - 360}
      if(pga_rot100ang < 0) {pga_rot100ang <- pga_rot100ang + 360}
      pgv_rot00ang <- pgv_rot00ang + ang1
      if(pgv_rot00ang > 360) {pgv_rot00ang <- pgv_rot00ang - 360}
      if(pgv_rot00ang < 0) {pgv_rot00ang <- pgv_rot00ang + 360}
      pgv_rot100ang <- pgv_rot100ang + ang1
      if(pgv_rot100ang > 360) {pgv_rot100ang <- pgv_rot100ang - 360}
      if(pgv_rot100ang < 0) {pgv_rot100ang <- pgv_rot100ang + 360}
      pgd_rot00ang <- pgd_rot00ang + ang1
      if(pgd_rot00ang > 360) {pgd_rot00ang <- pgd_rot00ang - 360}
      if(pgd_rot00ang < 0) {pgd_rot00ang <- pgd_rot00ang + 360}
      pgd_rot100ang <- pgd_rot100ang + ang1
      if(pgd_rot100ang > 360) {pgd_rot100ang <- pgd_rot100ang - 360}
      if(pgd_rot100ang < 0) {pgd_rot100ang <- pgd_rot100ang + 360}

      Value <- c(npts1, npts2, tmax4penalty_in, tmin4penalty_in, damping,
                 GMRotI50_ang, pga_gmrotI50, pgv_gmrotI50, pgd_gmrotI50, pga_gmrot50,
                 pgv_gmrot50, pga_gmrot100ang, pga_gmrot100, pgv_gmrot100ang, pgv_gmrot100,
                 pga_1_ar, pga_2_ar,n_a_subset, pgv_1_ar, pgv_2_ar, n_v_subset, pgd_1_ar,
                 pgd_2_ar, n_d_subset, pga_rot00, pga_rot50, pga_rot100, pga_rot00ang,
                 pga_rot100ang, pgv_rot00, pgv_rot50, pgv_rot100, pgv_rot00ang, pgv_rot100ang,
                 pgd_rot00, pgd_rot50, pgd_rot100, pgd_rot00ang, pgd_rot100ang)
      frame.data2 <- data.frame(Name, Value)
      filename1 <- paste(station_name, "_dep.csv", sep = "")
      filename2 <- paste(station_name, "_indep.csv", sep = "")
      setwd(outputdatadir)
      write.table(frame.data1, file = filename1, col.names = FALSE, sep = ',')
      write.csv(frame.data2, file = filename2, row.names = FALSE)

      #################### Plot curves  ########################################
      setwd(outputplotdir)
      png(filename = paste(station_name,"_PSA.png", sep = ""), width = 1280, height = 768)
      par(oma=c(2,2,2,2))
      par(mar=c(5,5,4,2) + 0.1)
      plot(c(0.01,20), c(0,max(max(PSA1,PSA2))), log = "x", type = "n", xlab = "Period(s)", ylab = "PSA(g)", main = "PSA versus Period", cex.main=2, cex.lab=2, cex.axis=2)
      lines(period_t, PSA_1, lty = 2, lwd = 2, col = "red")
      lines(period_t, PSA_2, lty = 3, lwd = 2, col = "blue")
      legend(5,max(max(PSA1,PSA2)), c("PSA_H1","PSA_H2"), lty=c(2,3), lwd=c(2,2),col=c("red","blue"), cex = 2)
      dev.off()
      png(filename = paste(station_name,"_GMRotI50.png", sep = ""), width = 1280, height = 768)
      par(oma=c(2,2,2,2))
      par(mar=c(5,5,4,2) + 0.1)
      plot(period_t, GMRotI50, log = "x", type = "l", lwd = 2, xlab = "Period(s)", ylab = "GMRotI50(g)", main = "GMRotI50 versus Period", cex.lab=2, cex.axis=2, cex.main=2)
      dev.off()
      png(filename = paste(station_name,"_RotD.png", sep = ""), width = 1280, height = 768)
      par(oma=c(2,2,2,2))
      par(mar=c(5,5,4,2) + 0.1)
      plot(c(0.01,20), c(0,max(RotD100)), log = "x", type = "n", xlab = "Period(s)", ylab = "RotDxx(g)", main = "RotDxx versus Period", cex.lab=2, cex.axis=2, cex.main=2)
      lines(period_t, RotD00,  lty = 2, lwd = 2, col = 2)
      lines(period_t, RotD50,  lty = 3, lwd = 2, col = 3)
      lines(period_t, RotD100, lty = 4, lwd = 2, col = 4)
      legend(5, max(RotD100), c("RotD00","RotD50","RotD100"), lty = c(2,3,4),  lwd = c(2,2,2), col = c(2,3,4), cex = 2)
      dev.off()

      }

    ###################### Create a combinated spreadsheet ########################
    ###############################################################################
    name_com <- c('Station_Name')
    # frame.data_com1 <- data.frame(Station)
    # frame.data_com1 <- rbind(name_com, frame.data_com1)
    frame.data_com1 <- rbind(name_com, t(t(Station)))
    name_com <- c(paste('PGA_Rot', combine_index, sep = ''), paste('PGV_Rot', combine_index, sep = ''), paste('PGD_Rot', combine_index, sep = ''), Period)
    value_com <- cbind(pga_combine, pgv_combine, pgd_combine, rotd_combine)
    frame.data_com2 <- t(data.frame(name_com, t(value_com)))
    frame.data_com <- cbind(frame.data_com1, frame.data_com2)
    filename <- paste('summary', ".csv", sep = "")
    setwd(outputdatadir)
    write.table(frame.data_com, file = filename, col.names = FALSE, row.names = FALSE, sep = ',')

  }


  ############### process USGS SMC formated data #################################
  else if(datatype == 'smc'){
    setwd(inputpath)
    temp = list.files()
    if(!length(temp)%%2){            # judge if it is even number of groups
      ## define combined spreadsheet variables ##
      Station <- seq(1, length(temp)/2)
      pga_combine <- seq(1, length(temp)/2)*0
      pgv_combine <- seq(1, length(temp)/2)*0
      pgd_combine <- seq(1, length(temp)/2)*0
      rotd_combine <- matrix(seq(1, (length(temp)/2)*length(period_t)) ,nrow = length(temp)/2, ncol = length(period_t))*0

      for(i in 1:(length(temp)/2)){
        setwd(inputpath)
        name1 <- basename(temp[2*i-1])
        temp_name1 <- strsplit(name1, split = "_")
        name2 <- basename(temp[2*i])
        temp_name2 <- strsplit(name2, split = "_")
        if(temp_name1[[1]][1] != temp_name2[[1]][1]){      # judge if there is a coupled group for name1
          print(c("You should put the coupled-group in the Inputdata folder!", name1, "is not coupled. Quit!"))
          stop()
        }
        station_name <- paste(temp_name1[[1]][1], temp_name1[[1]][2], sep = "_")
        Station[i] <- station_name
        commlines1 <- scan(temp[2*i-1], sep = "", skip = 12, nlines = 1, what = character())
        numcomm1 <- as.numeric(commlines1)[8]
        commlines1 <- scan(temp[2*i-1], sep = "", skip = 13, nlines = 1, what = character())
        npts1 <- as.numeric(commlines1)[1]
        commlines1 <- scan(temp[2*i-1], sep = "", skip = 17, nlines = 1, what = character())
        dt1 <- 1/as.numeric(commlines1)[2]

        commlines2 <- scan(temp[2*i], sep = "", skip = 12, nlines = 1, what = character())
        numcomm2 <- as.numeric(commlines2)[8]
        commlines2 <- scan(temp[2*i], sep = "", skip = 13, nlines = 1, what = character())
        npts2 <- as.numeric(commlines2)[1]
        commlines2 <- scan(temp[2*i], sep = "", skip = 17, nlines = 1, what = character())
        dt2 <- 1/as.numeric(commlines2)[2]

        if(dt1 != dt2){
          print("Samples per second are different! Quit!")
          stop()
        }

        line1 <- readLines(temp[2*i-1])
        data1 <- c()
        for(k in seq(28 + numcomm1, length(line1))){
          if(length(data1) == 0){
            data1 <- as.numeric(substring(line1[k], c(1,11,21,31,41,51,61,71),c(10,20,30,40,50,60,70,80)))
          } else{
            data1 <- c(data1,as.numeric(substring(line1[k], c(1,11,21,31,41,51,61,71),c(10,20,30,40,50,60,70,80))))
          }
        }

        line2 <- readLines(temp[2*i])
        data2 <- c()
        for(k in seq(28 + numcomm2, length(line2))){
          if(length(data2) == 0){
            data2 <- as.numeric(substring(line2[k], c(1,11,21,31,41,51,61,71),c(10,20,30,40,50,60,70,80)))
          } else{
            data2 <- c(data2,as.numeric(substring(line2[k], c(1,11,21,31,41,51,61,71),c(10,20,30,40,50,60,70,80))))
          }
        }

        # remove NAs
        data1 <- data1[!is.na(data1)]
        data2 <- data2[!is.na(data2)]

        ############################# Process data Begin ###################################
        ######### Sinc interpolate data ###############
        if(Interpolation_factor == 'auto'){
          # interp_factor <- ceiling(log(dt1/0.001)/log(2))
          # 1) the interpolation factor is based on f_nyq first -- we only care T_nyq,
          # because T less than T_nyq is not affected by interpolation,
          # therefore, dt/dt' = dt/(T_nyq/10) = dt/(dt/5) = 5, -> IF = 8
          # 2) if dt1 is less than the above calcualtion, then we use smaller interpolation
          interp_factor <- min(ceiling(log(dt1/0.001)/log(2)), 3)
          # if dt1 even smaller than the diresed time step (not common), we do not interpolate
          interp_factor <- max(interp_factor, 0)
          interp_factor <- 2^interp_factor
        }else if(log(Interpolation_factor)%%log(2) != 0){
          print("The Interpolation factor is not a number of power of 2, please try correct it!")
          stop()
        }else{
          interp_factor <- Interpolation_factor
        }

        dt1 <- dt1/interp_factor # compute new time step
        dt2 <- dt2/interp_factor
        data1 <- Interpft(data1, interp_factor) # new data after Sinc interpolation
        data2 <- Interpft(data2, interp_factor)
        number <- min(length(data1), length(data2)) # if two data sizes diff, we set both size as the smaller
        data1 <- data1[1:number]
        data2 <- data2[1:number]

        ############ Output PGA, PGV, PGD ###############
        ############ unit of orignial accel is cm/sec/sec
        ### change unit into g
        data1 <- data1/981 # unit, g
        data2 <- data2/981 # unit, g
        pga_1_ar <- max(abs(data1)) # unit, g
        pga_2_ar <- max(abs(data2)) # unit, g
        vel_1 <- cumtrapz(dt1*seq(1,length(data1)),data1) * 981 # unit, cm/s
        vel_2 <- cumtrapz(dt2*seq(1,length(data2)),data2) * 981 # unit, cm/s
        pgv_1_ar <- max(abs(vel_1))
        pgv_2_ar <- max(abs(vel_2))
        disp_1 <- cumtrapz(dt1*seq(1,length(vel_1)),vel_1) # unit, cm
        disp_2 <- cumtrapz(dt2*seq(1,length(vel_2)),vel_2) # unit, cm
        pgd_1_ar <- max(abs(disp_1))
        pgd_2_ar <- max(abs(disp_2))


        ############## Rotation values calculation ##############
        ######## subset slection #######
        length_min <- min(length(data1), length(data2))
        a_subset <- subset_select(data1, data2, fraction, length_min, dt1, 1)
        v_subset <- subset_select(vel_1, vel_2, fraction, length_min, dt1, 1)
        d_subset <- subset_select(disp_1, disp_2, fraction, length_min, dt1, 1)
        # define vectors
        pga_rot <- seq(1, 180)
        pga_gmrot <- seq(1,90)
        pgv_rot <- seq(1, 180)
        pgv_gmrot <- seq(1,90)
        pgd_rot <- seq(1, 180)
        pgd_gmrot <- seq(1,90)
        # calculation
        for(theta in seq(1,90)){
          # acceleration
          Rot_a1 <- a_subset[1,]*cos(theta/180*pi) + a_subset[2,]*sin(theta/180*pi)
          Rot_a2 <- -a_subset[1,]*sin(theta/180*pi) + a_subset[2,]*cos(theta/180*pi)
          pga_rot[theta] <- max(abs(Rot_a1))
          pga_rot[theta+90] <- max(abs(Rot_a2))
          pga_gmrot[theta] <- sqrt(max(abs(Rot_a1))*max(abs(Rot_a2)))
          # velocity
          Rot_v1 <- v_subset[1,]*cos(theta/180*pi) + v_subset[2,]*sin(theta/180*pi)
          Rot_v2 <- -v_subset[1,]*sin(theta/180*pi) + v_subset[2,]*cos(theta/180*pi)
          pgv_rot[theta] <- max(abs(Rot_v1))
          pgv_rot[theta+90] <- max(abs(Rot_v2))
          pgv_gmrot[theta] <- sqrt(max(abs(Rot_v1))*max(abs(Rot_v2)))
          # displacement
          Rot_d1 <- d_subset[1,]*cos(theta/180*pi) + d_subset[2,]*sin(theta/180*pi)
          Rot_d2 <- -d_subset[1,]*sin(theta/180*pi) + d_subset[2,]*cos(theta/180*pi)
          pgd_rot[theta] <- max(abs(Rot_d1))
          pgd_rot[theta+90] <- max(abs(Rot_d2))
          pgd_gmrot[theta] <- sqrt(max(abs(Rot_d1))*max(abs(Rot_d2)))
        }
        ## acceleration
        pga_rot00 <- min(pga_rot)
        pga_rot50 <- median(pga_rot)
        pga_rot100 <- max(pga_rot)
        pga_rot00ang <- which.min(pga_rot)
        pga_rot100ang <- which.max(pga_rot)
        pga_gmrot50 <- median(pga_gmrot)
        pga_gmrot100 <- max(pga_gmrot)
        pga_gmrot100ang <- which.max(pga_gmrot)
        ## velocity
        pgv_rot00 <- min(pgv_rot)
        pgv_rot50 <- median(pgv_rot)
        pgv_rot100 <- max(pgv_rot)
        pgv_rot00ang <- which.min(pgv_rot)
        pgv_rot100ang <- which.max(pgv_rot)
        pgv_gmrot50 <- median(pgv_gmrot)
        pgv_gmrot100 <- max(pgv_gmrot)
        pgv_gmrot100ang <- which.max(pgv_gmrot)
        ## displacement
        pgd_rot00 <- min(pgd_rot)
        pgd_rot50 <- median(pgd_rot)
        pgd_rot100 <- max(pgd_rot)
        pgd_rot00ang <- which.min(pgd_rot)
        pgd_rot100ang <- which.max(pgd_rot)
        pgd_gmrot50 <- median(pgd_gmrot)
        pgd_gmrot100 <- max(pgd_gmrot)
        pgd_gmrot100ang <- which.max(pgd_gmrot)

        ################# Applying PS_cal function to get PSA,PSV #############
        ########################## plug in the parameters ########################
        ##########################################
        # plug in the parameters for two components
        # using C++ to compute to speed up calculation
        # PSA calculation is from textbook of "Dynamics of Structure" 3rd edition by Chopra on page 168-169.
        ##### compute pseudo spectral acceleration
        PSA1 <- PS_cal_cpp(data1, period_t, damping, dt1, 1) # first row is absolute acceleration of oscillator, second row is psa
        PSA2 <- PS_cal_cpp(data2, period_t, damping, dt1, 1) # first row is absolute acceleration of oscillator, second row is psa
        PSA_gm_ar <- sqrt(PSA1[2,]*PSA2[2,])
        PSA_larger_ar <- seq(1, length(period_t))
        for(j in seq(1, length(PSA_larger_ar))){
          PSA_larger_ar[j] <- max(PSA1[2,j], PSA2[2,j])
        }

        ############## Applying GMRotD_cal function to get PSA for different periods for each rotated angle ####
        ########################## plug in the parameters ########################
        ##########################################################################
        results <- GM_RotD_cal(data1, data2, period_t, damping, dt1, fraction)  # the result is list, 1st is RotD180, 2nd is GMRotD
        RotD180 <- results[[1]]  # it is a list
        RotD180 <- matrix(unlist(RotD180), nrow = 180) # transform it into matrix
        GMRotD <- results[[2]]  # it is a list
        GMRotD <- matrix(unlist(GMRotD), nrow = 90) # transform it into matrix
        Num_points <- results[[3]]  # it is a list
        Num_points <- c(unlist(Num_points)) # transform it into matrix
        rd_alevel <- results[[4]]

        ##################### Applying GMRotI_cal function to get PSA for different periods for only one desired rotated angle
        ########################## plug in the parameters ########################
        ##########################################################################
        Result <- GMRotI50_cal(GMRotD, period_t, tmax = tmax4penalty_in, tmin = tmin4penalty_in)
        GMRotI50 <- Result[[1]]
        GMRotI50_ang <- Result[[2]]
        penalty <- Result[[3]]
        # acceleration
        pga_gmrotI50 <- pga_gmrot[GMRotI50_ang]
        # velocity
        pgv_gmrotI50 <- pgv_gmrot[GMRotI50_ang]
        # displacement
        pgd_gmrotI50 <- pgd_gmrot[GMRotI50_ang]

        ##################### Applying RotD_cal function to get PSA for different periods for each rotated angle
        ########################## plug in the parameters ########################
        ##########################################################################
        RotD00 <- apply(RotD180, 2, min)
        RotD100 <- apply(RotD180, 2, max)
        RotD50 <- apply(RotD180, 2, median)
        RotD_t <- t(RotD180) # transpose of RotD matrix, it is used to get the max index
        RotD100_ang <- max.col(RotD_t)
        RotD00_ang <- max.col(-RotD_t)

        ############################# Process data End ###################################
        ###############################################################################
        if(combine_index == 0){
          pga_combine[i] <- pga_rot00
          pgv_combine[i] <- pgv_rot00
          pgd_combine[i] <- pgd_rot00
          rotd_combine[i,] <- RotD00
        } else if(combine_index == 50){
          pga_combine[i] <- pga_rot50
          pgv_combine[i] <- pgv_rot50
          pgd_combine[i] <- pgd_rot50
          rotd_combine[i,] <- RotD50
        } else{
          pga_combine[i] <- pga_rot100
          pgv_combine[i] <- pgv_rot100
          pgd_combine[i] <- pgd_rot100
          rotd_combine[i,] <- RotD100
        }

        #################### Output calculation results ####################
        PSA_1 <- PSA1[2,]
        Absolute_acc_1 <- PSA1[1,] # absolute acceleration of oscillator
        PSA_2 <- PSA2[2,]
        Absolute_acc_2 <- PSA2[1,]
        npts1 <- length(data1)
        npts2 <- length(data2)
        n_a_subset <- length(a_subset[1,])
        n_v_subset <- length(v_subset[1,])
        n_d_subset <- length(d_subset[1,])
        Period <- c("T0.010S","T0.020S","T0.022S","T0.025S","T0.029S","T0.030S",
                    "T0.032S","T0.035S","T0.036S","T0.040S","T0.042S","T0.044S",
                    "T0.045S","T0.046S","T0.048S","T0.050S","T0.055S","T0.060S",
                    "T0.065S","T0.067S","T0.070S","T0.075S","T0.080S","T0.085S",
                    "T0.090S","T0.095S","T0.100S","T0.110S","T0.120S","T0.130S",
                    "T0.133S","T0.140S","T0.150S","T0.160S","T0.170S","T0.180S",
                    "T0.190S","T0.200S","T0.220S","T0.240S","T0.250S","T0.260S",
                    "T0.280S","T0.290S","T0.300S","T0.320S","T0.340S","T0.350S",
                    "T0.360S","T0.380S","T0.400S","T0.420S","T0.440S","T0.450S",
                    "T0.460S","T0.480S","T0.500S","T0.550S","T0.600S","T0.650S",
                    "T0.667S","T0.700S","T0.750S","T0.800S","T0.850S","T0.900S",
                    "T0.950S","T1.000S","T1.100S","T1.200S","T1.300S","T1.400S",
                    "T1.500S","T1.600S","T1.700S","T1.800S","T1.900S","T2.000S",
                    "T2.200S","T2.400S","T2.500S","T2.600S","T2.800S","T3.000S",
                    "T3.200S","T3.400S","T3.500S","T3.600S","T3.800S","T4.000S",
                    "T4.200S","T4.400S","T4.600S","T4.800S","T5.000S","T5.500S",
                    "T6.000S","T6.500S","T7.000S","T7.500S","T8.000S","T8.500S",
                    "T9.000S","T9.500S","T10.000S","T11.000S","T12.000S","T13.000S",
                    "T14.000S","T15.000S","T20.000S")

        # compute real azimuth and set angle in [0,360]
        RotD00_ang <- RotD00_ang + ang1
        RotD100_ang <- RotD100_ang + ang1
        RotD00_ang[which(RotD00_ang > 360)] <- RotD00_ang[which(RotD00_ang > 360)] - 360
        RotD00_ang[which(RotD00_ang < 0)] <- RotD00_ang[which(RotD00_ang < 0)] + 360
        RotD100_ang[which(RotD100_ang > 360)] <- RotD100_ang[which(RotD100_ang > 360)] - 360
        RotD100_ang[which(RotD100_ang < 0)] <- RotD100_ang[which(RotD100_ang < 0)] + 360

        frame.data1 <- data.frame(Period, PSA_1, Absolute_acc_1, PSA_2, Absolute_acc_2,
                                  PSA_gm_ar, PSA_larger_ar, GMRotI50, RotD00, RotD50,
                                  RotD100, RotD00_ang, RotD100_ang, Num_points)
        sub_frame.data1 <- t(frame.data1[,-1])
        frame.data1 <- rbind(t(frame.data1[1]), sub_frame.data1)
        Name <- c("npts1","npts2","tmax4penalty","tmin4penalty","damping", "GMRotI50angle",
                  "PGA_GMRotI50", "PGV_GMRotI50",
                  "PGD_GMRotI50", "PGA_GMRot50", "PGV_GMRot50", "PGA_GMRot100angle",
                  "PGA_GMRot100", "PGV_GMRot100angle", "PGV_GMRot100", "PGA_1", "PGA_2",
                  "n_a_subset","PGV_1", "PGV_2", "n_v_subset","PGD_1", "PGD_2", "n_d_subset",
                  "PGA_Rot00", "PGA_Rot50", "PGA_Rot100", "PGA_Rot00angle", "PGA_Rot100angle",
                  "PGV_Rot00", "PGV_Rot50", "PGV_Rot100", "PGV_Rot00angle", "PGV_Rot100angle",
                  "PGD_Rot00", "PGD_Rot50", "PGD_Rot100", "PGD_Rot00angle", "PGD_Rot100angle")

        # compute real azimuth and set angel in [0,360]
        GMRotI50_ang <- GMRotI50_ang + ang1
        if(GMRotI50_ang > 360) {GMRotI50_ang <- GMRotI50_ang - 360}
        if(GMRotI50_ang < 0) {GMRotI50_ang <- GMRotI50_ang + 360}
        pga_gmrot100ang <- pga_gmrot100ang + ang1
        if(pga_gmrot100ang > 360) {pga_gmrot100ang <- pga_gmrot100ang - 360}
        if(pga_gmrot100ang < 0) {pga_gmrot100ang <- pga_gmrot100ang + 360}
        pgv_gmrot100ang <- pgv_gmrot100ang + ang1
        if(pgv_gmrot100ang > 360) {pgv_gmrot100ang <- pgv_gmrot100ang - 360}
        if(pgv_gmrot100ang < 0) {pgv_gmrot100ang <- pgv_gmrot100ang + 360}
        pga_rot00ang <- pga_rot00ang + ang1
        if(pga_rot00ang > 360) {pga_rot00ang <- pga_rot00ang - 360}
        if(pga_rot00ang < 0) {pga_rot00ang <- pga_rot00ang + 360}
        pga_rot100ang <- pga_rot100ang + ang1
        if(pga_rot100ang > 360) {pga_rot100ang <- pga_rot100ang - 360}
        if(pga_rot100ang < 0) {pga_rot100ang <- pga_rot100ang + 360}
        pgv_rot00ang <- pgv_rot00ang + ang1
        if(pgv_rot00ang > 360) {pgv_rot00ang <- pgv_rot00ang - 360}
        if(pgv_rot00ang < 0) {pgv_rot00ang <- pgv_rot00ang + 360}
        pgv_rot100ang <- pgv_rot100ang + ang1
        if(pgv_rot100ang > 360) {pgv_rot100ang <- pgv_rot100ang - 360}
        if(pgv_rot100ang < 0) {pgv_rot100ang <- pgv_rot100ang + 360}
        pgd_rot00ang <- pgd_rot00ang + ang1
        if(pgd_rot00ang > 360) {pgd_rot00ang <- pgd_rot00ang - 360}
        if(pgd_rot00ang < 0) {pgd_rot00ang <- pgd_rot00ang + 360}
        pgd_rot100ang <- pgd_rot100ang + ang1
        if(pgd_rot100ang > 360) {pgd_rot100ang <- pgd_rot100ang - 360}
        if(pgd_rot100ang < 0) {pgd_rot100ang <- pgd_rot100ang + 360}

        Value <- c(npts1, npts2, tmax4penalty_in, tmin4penalty_in, damping,
                   GMRotI50_ang, pga_gmrotI50, pgv_gmrotI50, pgd_gmrotI50, pga_gmrot50,
                   pgv_gmrot50, pga_gmrot100ang, pga_gmrot100, pgv_gmrot100ang, pgv_gmrot100,
                   pga_1_ar, pga_2_ar,n_a_subset, pgv_1_ar, pgv_2_ar, n_v_subset, pgd_1_ar,
                   pgd_2_ar, n_d_subset, pga_rot00, pga_rot50, pga_rot100, pga_rot00ang,
                   pga_rot100ang, pgv_rot00, pgv_rot50, pgv_rot100, pgv_rot00ang, pgv_rot100ang,
                   pgd_rot00, pgd_rot50, pgd_rot100, pgd_rot00ang, pgd_rot100ang)
        frame.data2 <- data.frame(Name, Value)
        filename1 <- paste(station_name, "_dep.csv", sep = "")
        filename2 <- paste(station_name, "_indep.csv", sep = "")
        setwd(outputdatadir)
        write.table(frame.data1, file = filename1, col.names = FALSE, sep = ',')
        write.csv(frame.data2, file = filename2, row.names = FALSE)

        #################### Plot curves  ########################################
        setwd(outputplotdir)
        png(filename = paste(station_name,"_PSA.png", sep = ""), width = 1280, height = 768)
        par(oma=c(2,2,2,2))
        par(mar=c(5,5,4,2) + 0.1)
        plot(c(0.01,20), c(0,max(max(PSA1,PSA2))), log = "x", type = "n", xlab = "Period(s)", ylab = "PSA(g)", main = "PSA versus Period", cex.main=2, cex.lab=2, cex.axis=2)
        lines(period_t, PSA_1, lty = 2, lwd = 2, col = "red")
        lines(period_t, PSA_2, lty = 3, lwd = 2, col = "blue")
        legend(5,max(max(PSA1,PSA2)), c("PSA_H1","PSA_H2"), lty=c(2,3), lwd=c(2,2),col=c("red","blue"), cex = 2)
        dev.off()
        png(filename = paste(station_name,"_GMRotI50.png", sep = ""), width = 1280, height = 768)
        par(oma=c(2,2,2,2))
        par(mar=c(5,5,4,2) + 0.1)
        plot(period_t, GMRotI50, log = "x", type = "l", lwd = 2, xlab = "Period(s)", ylab = "GMRotI50(g)", main = "GMRotI50 versus Period", cex.lab=2, cex.axis=2, cex.main=2)
        dev.off()
        png(filename = paste(station_name,"_RotD.png", sep = ""), width = 1280, height = 768)
        par(oma=c(2,2,2,2))
        par(mar=c(5,5,4,2) + 0.1)
        plot(c(0.01,20), c(0,max(RotD100)), log = "x", type = "n", xlab = "Period(s)", ylab = "RotDxx(g)", main = "RotDxx versus Period", cex.lab=2, cex.axis=2, cex.main=2)
        lines(period_t, RotD00,  lty = 2, lwd = 2, col = 2)
        lines(period_t, RotD50,  lty = 3, lwd = 2, col = 3)
        lines(period_t, RotD100, lty = 4, lwd = 2, col = 4)
        legend(5, max(RotD100), c("RotD00","RotD50","RotD100"), lty = c(2,3,4),  lwd = c(2,2,2), col = c(2,3,4), cex = 2)
        dev.off()

      }

    } else{
      print("You should put the couple groups in the Inputdata folder, we detect you have odd number of time series")
      stop()
    }

    ###################### Create a combinated spreadsheet ########################
    ###############################################################################
    name_com <- c('Station_Name')
    # frame.data_com1 <- data.frame(Station)
    # frame.data_com1 <- rbind(name_com, frame.data_com1)
    frame.data_com1 <- rbind(name_com, t(t(Station)))
    name_com <- c(paste('PGA_Rot', combine_index, sep = ''), paste('PGV_Rot', combine_index, sep = ''), paste('PGD_Rot', combine_index, sep = ''), Period)
    value_com <- cbind(pga_combine, pgv_combine, pgd_combine, rotd_combine)
    frame.data_com2 <- t(data.frame(name_com, t(value_com)))
    frame.data_com <- cbind(frame.data_com1, frame.data_com2)
    filename <- paste('summary', ".csv", sep = "")
    setwd(outputdatadir)
    write.table(frame.data_com, file = filename, col.names = FALSE, row.names = FALSE, sep = ',')



  }


  ######## process only one column data (the first value is time step) #########
  else if(datatype == 'timeseries'){
    ### only one column data, the first vlaue of the column is time step
    setwd(inputpath)
    temp = list.files()
    if(!length(temp)%%2){            # judge if it is even number of groups
      ## define combined spreadsheet variables ##
      Station <- seq(1, length(temp)/2)
      pga_combine <- seq(1, length(temp)/2)*0
      pgv_combine <- seq(1, length(temp)/2)*0
      pgd_combine <- seq(1, length(temp)/2)*0
      rotd_combine <- matrix(seq(1, (length(temp)/2)*length(period_t)) ,nrow = length(temp)/2, ncol = length(period_t))*0

      for(i in 1:(length(temp)/2)){
        setwd(inputpath)
        name1 <- basename(temp[2*i-1])
        temp_name1 <- strsplit(name1, split = "_")
        name2 <- basename(temp[2*i])
        temp_name2 <- strsplit(name2, split = "_")
        if(temp_name1[[1]][1] != temp_name2[[1]][1]){      # judge if there is a coupled group for name1
          print(c("You should put the coupled-group in the Inputdata folder!", name1, "is not coupled. Quit!"))
          stop()
        }
        station_name <- paste(temp_name1[[1]][1], temp_name1[[1]][2], sep = "_")
        Station[i] <- station_name
        dt1 <- as.numeric(readLines(temp[2*i-1], n = 1))
        dt2 <- as.numeric(readLines(temp[2*i], n = 1))
        if(dt1 != dt2){
          print("Samples per second are different! Quit!")
          stop()
        }
        data1 <- as.numeric(readLines(temp[2*i-1]))
        data1 <- data1[2:length(data1)]
        data2 <- as.numeric(readLines(temp[2*i]))
        data2 <- data2[2:length(data2)]


        ############################# Process data Begin ###################################
        ######### Sinc interpolate data ###############
        if(Interpolation_factor == 'auto'){
          # interp_factor <- ceiling(log(dt1/0.001)/log(2))

          # 1) the interpolation factor is based on f_nyq first -- we only care T_nyq,
          # because T less than T_nyq is not affected by interpolation,
          # therefore, dt/dt' = dt/(T_nyq/10) = dt/(dt/5) = 5, -> IF = 8
          # 2) if dt1 is less than the above calcualtion, then we use smaller interpolation
          interp_factor <- min(ceiling(log(dt1/0.001)/log(2)), 3)
          # if dt1 even smaller than the diresed time step (not common), we do not interpolate
          interp_factor <- max(interp_factor, 0)
          interp_factor <- 2^interp_factor
        }else if(log(Interpolation_factor)%%log(2) != 0){
          print("The Interpolation factor is not a number of power of 2, please try correct it!")
          stop()
        }else{
          interp_factor <- Interpolation_factor
        }

        dt1 <- dt1/interp_factor # compute new time step
        dt2 <- dt2/interp_factor
        data1 <- Interpft(data1, interp_factor) # new data after Sinc interpolation
        data2 <- Interpft(data2, interp_factor)
        number <- min(length(data1), length(data2)) # if two data sizes diff, we set both size as the smaller
        data1 <- data1[1:number]
        data2 <- data2[1:number]

        ############ Output PGA, PGV, PGD ###############
        ########## unit of original accel is g
        pga_1_ar <- max(abs(data1)) # unit, g
        pga_2_ar <- max(abs(data2)) # unit, g
        vel_1 <- cumtrapz(dt1*seq(1,length(data1)),data1) * 981 # unit, cm/s
        vel_2 <- cumtrapz(dt2*seq(1,length(data2)),data2) * 981 # unit, cm/s
        pgv_1_ar <- max(abs(vel_1))
        pgv_2_ar <- max(abs(vel_2))
        disp_1 <- cumtrapz(dt1*seq(1,length(vel_1)),vel_1) # unit, cm
        disp_2 <- cumtrapz(dt2*seq(1,length(vel_2)),vel_2) # unit, cm
        pgd_1_ar <- max(abs(disp_1))
        pgd_2_ar <- max(abs(disp_2))

        ############## Rotation values calculation ##############
        ######## subset slection #######
        length_min <- min(length(data1), length(data2))
        a_subset <- subset_select(data1, data2, fraction, length_min, dt1, 1)
        v_subset <- subset_select(vel_1, vel_2, fraction, length_min, dt1, 1)
        d_subset <- subset_select(disp_1, disp_2, fraction, length_min, dt1, 1)
        # define vectors
        pga_rot <- seq(1, 180)
        pga_gmrot <- seq(1,90)
        pgv_rot <- seq(1, 180)
        pgv_gmrot <- seq(1,90)
        pgd_rot <- seq(1, 180)
        pgd_gmrot <- seq(1,90)
        # calculation
        for(theta in seq(1,90)){
          # acceleration
          Rot_a1 <- a_subset[1,]*cos(theta/180*pi) + a_subset[2,]*sin(theta/180*pi)
          Rot_a2 <- -a_subset[1,]*sin(theta/180*pi) + a_subset[2,]*cos(theta/180*pi)
          pga_rot[theta] <- max(abs(Rot_a1))
          pga_rot[theta+90] <- max(abs(Rot_a2))
          pga_gmrot[theta] <- sqrt(max(abs(Rot_a1))*max(abs(Rot_a2)))
          # velocity
          Rot_v1 <- v_subset[1,]*cos(theta/180*pi) + v_subset[2,]*sin(theta/180*pi)
          Rot_v2 <- -v_subset[1,]*sin(theta/180*pi) + v_subset[2,]*cos(theta/180*pi)
          pgv_rot[theta] <- max(abs(Rot_v1))
          pgv_rot[theta+90] <- max(abs(Rot_v2))
          pgv_gmrot[theta] <- sqrt(max(abs(Rot_v1))*max(abs(Rot_v2)))
          # displacement
          Rot_d1 <- d_subset[1,]*cos(theta/180*pi) + d_subset[2,]*sin(theta/180*pi)
          Rot_d2 <- -d_subset[1,]*sin(theta/180*pi) + d_subset[2,]*cos(theta/180*pi)
          pgd_rot[theta] <- max(abs(Rot_d1))
          pgd_rot[theta+90] <- max(abs(Rot_d2))
          pgd_gmrot[theta] <- sqrt(max(abs(Rot_d1))*max(abs(Rot_d2)))
        }
        ## acceleration
        pga_rot00 <- min(pga_rot)
        pga_rot50 <- median(pga_rot)
        pga_rot100 <- max(pga_rot)
        pga_rot00ang <- which.min(pga_rot)
        pga_rot100ang <- which.max(pga_rot)
        pga_gmrot50 <- median(pga_gmrot)
        pga_gmrot100 <- max(pga_gmrot)
        pga_gmrot100ang <- which.max(pga_gmrot)
        ## velocity
        pgv_rot00 <- min(pgv_rot)
        pgv_rot50 <- median(pgv_rot)
        pgv_rot100 <- max(pgv_rot)
        pgv_rot00ang <- which.min(pgv_rot)
        pgv_rot100ang <- which.max(pgv_rot)
        pgv_gmrot50 <- median(pgv_gmrot)
        pgv_gmrot100 <- max(pgv_gmrot)
        pgv_gmrot100ang <- which.max(pgv_gmrot)
        ## displacement
        pgd_rot00 <- min(pgd_rot)
        pgd_rot50 <- median(pgd_rot)
        pgd_rot100 <- max(pgd_rot)
        pgd_rot00ang <- which.min(pgd_rot)
        pgd_rot100ang <- which.max(pgd_rot)
        pgd_gmrot50 <- median(pgd_gmrot)
        pgd_gmrot100 <- max(pgd_gmrot)
        pgd_gmrot100ang <- which.max(pgd_gmrot)

        ################# Applying PS_cal function to get PSA,PSV #############
        ########################## plug in the parameters ########################
        ##########################################
        # plug in the parameters for two components
        # using C++ to compute to speed up calculation
        # PSA calculation is from textbook of "Dynamics of Structure" 3rd edition by Chopra on page 168-169.
        ##### compute pseudo spectral acceleration
        PSA1 <- PS_cal_cpp(data1, period_t, damping, dt1, 1) # first row is absolute acceleration of oscillator, second row is psa
        PSA2 <- PS_cal_cpp(data2, period_t, damping, dt1, 1) # first row is absolute acceleration of oscillator, second row is psa
        PSA_gm_ar <- sqrt(PSA1[2,]*PSA2[2,])
        PSA_larger_ar <- seq(1, length(period_t))
        for(j in seq(1, length(PSA_larger_ar))){
          PSA_larger_ar[j] <- max(PSA1[2,j], PSA2[2,j])
        }

        ############## Applying GMRotD_cal function to get PSA for different periods for each rotated angle ####
        ########################## plug in the parameters ########################
        ##########################################################################
        results <- GM_RotD_cal(data1, data2, period_t, damping, dt1, fraction)  # the result is list, 1st is RotD180, 2nd is GMRotD
        RotD180 <- results[[1]]  # it is a list
        RotD180 <- matrix(unlist(RotD180), nrow = 180) # transform it into matrix
        GMRotD <- results[[2]]  # it is a list
        GMRotD <- matrix(unlist(GMRotD), nrow = 90) # transform it into matrix
        Num_points <- results[[3]]  # it is a list
        Num_points <- c(unlist(Num_points)) # transform it into matrix
        rd_alevel <- results[[4]]

        ##################### Applying GMRotI_cal function to get PSA for different periods for only one desired rotated angle
        ########################## plug in the parameters ########################
        ##########################################################################
        Result <- GMRotI50_cal(GMRotD, period_t, tmax = tmax4penalty_in, tmin = tmin4penalty_in)
        GMRotI50 <- Result[[1]]
        GMRotI50_ang <- Result[[2]]
        penalty <- Result[[3]]
        # acceleration
        pga_gmrotI50 <- pga_gmrot[GMRotI50_ang]
        # velocity
        pgv_gmrotI50 <- pgv_gmrot[GMRotI50_ang]
        # displacement
        pgd_gmrotI50 <- pgd_gmrot[GMRotI50_ang]

        ##################### Applying RotD_cal function to get PSA for different periods for each rotated angle
        ########################## plug in the parameters ########################
        ##########################################################################
        RotD00 <- apply(RotD180, 2, min)
        RotD100 <- apply(RotD180, 2, max)
        RotD50 <- apply(RotD180, 2, median)
        RotD_t <- t(RotD180) # transpose of RotD matrix, it is used to get the max index
        RotD100_ang <- max.col(RotD_t)
        RotD00_ang <- max.col(-RotD_t)

        ############################# Process data End ###################################
        ###############################################################################
        if(combine_index == 0){
          pga_combine[i] <- pga_rot00
          pgv_combine[i] <- pgv_rot00
          pgd_combine[i] <- pgd_rot00
          rotd_combine[i,] <- RotD00
        } else if(combine_index == 50){
          pga_combine[i] <- pga_rot50
          pgv_combine[i] <- pgv_rot50
          pgd_combine[i] <- pgd_rot50
          rotd_combine[i,] <- RotD50
        } else{
          pga_combine[i] <- pga_rot100
          pgv_combine[i] <- pgv_rot100
          pgd_combine[i] <- pgd_rot100
          rotd_combine[i,] <- RotD100
        }

        #################### Output calculation results ####################
        PSA_1 <- PSA1[2,]
        Absolute_acc_1 <- PSA1[1,] # absolute acceleration of oscillator
        PSA_2 <- PSA2[2,]
        Absolute_acc_2 <- PSA2[1,]
        npts1 <- length(data1)
        npts2 <- length(data2)
        n_a_subset <- length(a_subset[1,])
        n_v_subset <- length(v_subset[1,])
        n_d_subset <- length(d_subset[1,])
        Period <- c("T0.010S","T0.020S","T0.022S","T0.025S","T0.029S","T0.030S",
                    "T0.032S","T0.035S","T0.036S","T0.040S","T0.042S","T0.044S",
                    "T0.045S","T0.046S","T0.048S","T0.050S","T0.055S","T0.060S",
                    "T0.065S","T0.067S","T0.070S","T0.075S","T0.080S","T0.085S",
                    "T0.090S","T0.095S","T0.100S","T0.110S","T0.120S","T0.130S",
                    "T0.133S","T0.140S","T0.150S","T0.160S","T0.170S","T0.180S",
                    "T0.190S","T0.200S","T0.220S","T0.240S","T0.250S","T0.260S",
                    "T0.280S","T0.290S","T0.300S","T0.320S","T0.340S","T0.350S",
                    "T0.360S","T0.380S","T0.400S","T0.420S","T0.440S","T0.450S",
                    "T0.460S","T0.480S","T0.500S","T0.550S","T0.600S","T0.650S",
                    "T0.667S","T0.700S","T0.750S","T0.800S","T0.850S","T0.900S",
                    "T0.950S","T1.000S","T1.100S","T1.200S","T1.300S","T1.400S",
                    "T1.500S","T1.600S","T1.700S","T1.800S","T1.900S","T2.000S",
                    "T2.200S","T2.400S","T2.500S","T2.600S","T2.800S","T3.000S",
                    "T3.200S","T3.400S","T3.500S","T3.600S","T3.800S","T4.000S",
                    "T4.200S","T4.400S","T4.600S","T4.800S","T5.000S","T5.500S",
                    "T6.000S","T6.500S","T7.000S","T7.500S","T8.000S","T8.500S",
                    "T9.000S","T9.500S","T10.000S","T11.000S","T12.000S","T13.000S",
                    "T14.000S","T15.000S","T20.000S")

        # compute real azimuth and set angle in [0,360]
        RotD00_ang <- RotD00_ang + ang1
        RotD100_ang <- RotD100_ang + ang1
        RotD00_ang[which(RotD00_ang > 360)] <- RotD00_ang[which(RotD00_ang > 360)] - 360
        RotD00_ang[which(RotD00_ang < 0)] <- RotD00_ang[which(RotD00_ang < 0)] + 360
        RotD100_ang[which(RotD100_ang > 360)] <- RotD100_ang[which(RotD100_ang > 360)] - 360
        RotD100_ang[which(RotD100_ang < 0)] <- RotD100_ang[which(RotD100_ang < 0)] + 360

        frame.data1 <- data.frame(Period, PSA_1, Absolute_acc_1, PSA_2, Absolute_acc_2,
                                  PSA_gm_ar, PSA_larger_ar, GMRotI50, RotD00, RotD50,
                                  RotD100, RotD00_ang, RotD100_ang, Num_points)
        sub_frame.data1 <- t(frame.data1[,-1])
        frame.data1 <- rbind(t(frame.data1[1]), sub_frame.data1)
        Name <- c("npts1","npts2","tmax4penalty","tmin4penalty","damping", "GMRotI50angle",
                  "PGA_GMRotI50", "PGV_GMRotI50",
                  "PGD_GMRotI50", "PGA_GMRot50", "PGV_GMRot50", "PGA_GMRot100angle",
                  "PGA_GMRot100", "PGV_GMRot100angle", "PGV_GMRot100", "PGA_1", "PGA_2",
                  "n_a_subset","PGV_1", "PGV_2", "n_v_subset","PGD_1", "PGD_2", "n_d_subset",
                  "PGA_Rot00", "PGA_Rot50", "PGA_Rot100", "PGA_Rot00angle", "PGA_Rot100angle",
                  "PGV_Rot00", "PGV_Rot50", "PGV_Rot100", "PGV_Rot00angle", "PGV_Rot100angle",
                  "PGD_Rot00", "PGD_Rot50", "PGD_Rot100", "PGD_Rot00angle", "PGD_Rot100angle")

        # compute real azimuth and set angel in [0,360]
        GMRotI50_ang <- GMRotI50_ang + ang1
        if(GMRotI50_ang > 360) {GMRotI50_ang <- GMRotI50_ang - 360}
        if(GMRotI50_ang < 0) {GMRotI50_ang <- GMRotI50_ang + 360}
        pga_gmrot100ang <- pga_gmrot100ang + ang1
        if(pga_gmrot100ang > 360) {pga_gmrot100ang <- pga_gmrot100ang - 360}
        if(pga_gmrot100ang < 0) {pga_gmrot100ang <- pga_gmrot100ang + 360}
        pgv_gmrot100ang <- pgv_gmrot100ang + ang1
        if(pgv_gmrot100ang > 360) {pgv_gmrot100ang <- pgv_gmrot100ang - 360}
        if(pgv_gmrot100ang < 0) {pgv_gmrot100ang <- pgv_gmrot100ang + 360}
        pga_rot00ang <- pga_rot00ang + ang1
        if(pga_rot00ang > 360) {pga_rot00ang <- pga_rot00ang - 360}
        if(pga_rot00ang < 0) {pga_rot00ang <- pga_rot00ang + 360}
        pga_rot100ang <- pga_rot100ang + ang1
        if(pga_rot100ang > 360) {pga_rot100ang <- pga_rot100ang - 360}
        if(pga_rot100ang < 0) {pga_rot100ang <- pga_rot100ang + 360}
        pgv_rot00ang <- pgv_rot00ang + ang1
        if(pgv_rot00ang > 360) {pgv_rot00ang <- pgv_rot00ang - 360}
        if(pgv_rot00ang < 0) {pgv_rot00ang <- pgv_rot00ang + 360}
        pgv_rot100ang <- pgv_rot100ang + ang1
        if(pgv_rot100ang > 360) {pgv_rot100ang <- pgv_rot100ang - 360}
        if(pgv_rot100ang < 0) {pgv_rot100ang <- pgv_rot100ang + 360}
        pgd_rot00ang <- pgd_rot00ang + ang1
        if(pgd_rot00ang > 360) {pgd_rot00ang <- pgd_rot00ang - 360}
        if(pgd_rot00ang < 0) {pgd_rot00ang <- pgd_rot00ang + 360}
        pgd_rot100ang <- pgd_rot100ang + ang1
        if(pgd_rot100ang > 360) {pgd_rot100ang <- pgd_rot100ang - 360}
        if(pgd_rot100ang < 0) {pgd_rot100ang <- pgd_rot100ang + 360}

        Value <- c(npts1, npts2, tmax4penalty_in, tmin4penalty_in, damping,
                   GMRotI50_ang, pga_gmrotI50, pgv_gmrotI50, pgd_gmrotI50, pga_gmrot50,
                   pgv_gmrot50, pga_gmrot100ang, pga_gmrot100, pgv_gmrot100ang, pgv_gmrot100,
                   pga_1_ar, pga_2_ar,n_a_subset, pgv_1_ar, pgv_2_ar, n_v_subset, pgd_1_ar,
                   pgd_2_ar, n_d_subset, pga_rot00, pga_rot50, pga_rot100, pga_rot00ang,
                   pga_rot100ang, pgv_rot00, pgv_rot50, pgv_rot100, pgv_rot00ang, pgv_rot100ang,
                   pgd_rot00, pgd_rot50, pgd_rot100, pgd_rot00ang, pgd_rot100ang)
        frame.data2 <- data.frame(Name, Value)
        filename1 <- paste(station_name, "_dep.csv", sep = "")
        filename2 <- paste(station_name, "_indep.csv", sep = "")
        setwd(outputdatadir)
        write.table(frame.data1, file = filename1, col.names = FALSE, sep = ',')
        write.csv(frame.data2, file = filename2, row.names = FALSE)

        #################### Plot curves  ########################################
        setwd(outputplotdir)
        png(filename = paste(station_name,"_PSA.png", sep = ""), width = 1280, height = 768)
        par(oma=c(2,2,2,2))
        par(mar=c(5,5,4,2) + 0.1)
        plot(c(0.01,20), c(0,max(max(PSA1,PSA2))), log = "x", type = "n", xlab = "Period(s)", ylab = "PSA(g)", main = "PSA versus Period", cex.main=2, cex.lab=2, cex.axis=2)
        lines(period_t, PSA_1, lty = 2, lwd = 2, col = "red")
        lines(period_t, PSA_2, lty = 3, lwd = 2, col = "blue")
        legend(5,max(max(PSA1,PSA2)), c("PSA_H1","PSA_H2"), lty=c(2,3), lwd=c(2,2),col=c("red","blue"), cex = 2)
        dev.off()
        png(filename = paste(station_name,"_GMRotI50.png", sep = ""), width = 1280, height = 768)
        par(oma=c(2,2,2,2))
        par(mar=c(5,5,4,2) + 0.1)
        plot(period_t, GMRotI50, log = "x", type = "l", lwd = 2, xlab = "Period(s)", ylab = "GMRotI50(g)", main = "GMRotI50 versus Period", cex.lab=2, cex.axis=2, cex.main=2)
        dev.off()
        png(filename = paste(station_name,"_RotD.png", sep = ""), width = 1280, height = 768)
        par(oma=c(2,2,2,2))
        par(mar=c(5,5,4,2) + 0.1)
        plot(c(0.01,20), c(0,max(RotD100)), log = "x", type = "n", xlab = "Period(s)", ylab = "RotDxx(g)", main = "RotDxx versus Period", cex.lab=2, cex.axis=2, cex.main=2)
        lines(period_t, RotD00,  lty = 2, lwd = 2, col = 2)
        lines(period_t, RotD50,  lty = 3, lwd = 2, col = 3)
        lines(period_t, RotD100, lty = 4, lwd = 2, col = 4)
        legend(5, max(RotD100), c("RotD00","RotD50","RotD100"), lty = c(2,3,4),  lwd = c(2,2,2), col = c(2,3,4), cex = 2)
        dev.off()

      }

    } else{
      print("You should put the couple groups in the Inputdata folder, we detect you have odd number of time series")
      stop()
    }

    ###################### Create a combinated spreadsheet ########################
    ###############################################################################
    name_com <- c('Station_Name')
    # frame.data_com1 <- data.frame(Station)
    # frame.data_com1 <- rbind(name_com, frame.data_com1)
    frame.data_com1 <- rbind(name_com, t(t(Station)))
    name_com <- c(paste('PGA_Rot', combine_index, sep = ''), paste('PGV_Rot', combine_index, sep = ''), paste('PGD_Rot', combine_index, sep = ''), Period)
    value_com <- cbind(pga_combine, pgv_combine, pgd_combine, rotd_combine)
    frame.data_com2 <- t(data.frame(name_com, t(value_com)))
    frame.data_com <- cbind(frame.data_com1, frame.data_com2)
    filename <- paste('summary', ".csv", sep = "")
    setwd(outputdatadir)
    write.table(frame.data_com, file = filename, col.names = FALSE, row.names = FALSE, sep = ',')

  }

}






