#' A function to calculate a suite of intensity measures and generate plots
#'
#' This is a comprehensive function that computes response spectra and generates spreadsheets and plots.
#' The function computes PSA of RotD00, RotD50, RotD100, GMRotI50, and other periods independent
#' variables (e.g., PGA, PGV, PGD, etc.). The outputs will be two cvs files and plots for each station. A summary csv file is also generated.
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
#' @importFrom utils tail write.table write.csv
#' @importFrom stats median
#' @importFrom grDevices png dev.off
#' @importFrom graphics par lines legend axis
#' @return It generates two csv files and plots for each input station recordings. Two folders are generated in the
#' \code{inputpath} directory, one is for "Outputdata" and another is for "Outputplot"
#' @keywords IMplot
#' @export

IMplot <- function(inputpath='/Users/PFW/Desktop/test/Inputdata', datatype = "ngaw2", tmax4penalty_in=10, tmin4penalty_in=0, combine_index=50,
                   ang1 = 0, damping=0.05, fraction=0.7, Interpolation_factor='auto'){

  # create directory
  maindir <- dirname(inputpath)
  dir.create(file.path(maindir, 'Outputdata'))
  outputdatadir <- file.path(maindir, 'Outputdata')
  dir.create(file.path(maindir, 'Outputplot'))
  outputplotdir <- file.path(maindir, 'Outputplot')

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
  period_names <- c("T0.010S","T0.020S","T0.022S","T0.025S","T0.029S","T0.030S",
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


  ## process the downloaded NGA-West2 .AT2 files
  if(datatype == 'ngaw2'){
    temp = list.files(inputpath, pattern = "*.AT2")
    if(!length(temp)%%2){
      # running only if the input files are paired NGA-West2 .AT2 files

      # define combined spreadsheet variables
      Record_Sequence_Number <- seq(1, length(temp)/2)
      EQID <- seq(1, length(temp)/2)
      pga_combine <- seq(1, length(temp)/2)*0
      pgv_combine <- seq(1, length(temp)/2)*0
      pgd_combine <- seq(1, length(temp)/2)*0
      rotd_combine <- matrix(seq(1, (length(temp)/2)*length(period_t)) ,nrow = length(temp)/2, ncol = length(period_t))*0
      highest_freq_combine <- seq(1, length(temp)/2)
      lowest_freq_combine <- seq(1, length(temp)/2)

      for(i in 1:(length(temp)/2)){
        name1 <- basename(temp[2*i-1])
        temp_name1 <- strsplit(name1, split = "[_]")
        name2 <- basename(temp[2*i])
        temp_name2 <- strsplit(name2, split = "[_]")
        if(temp_name1[[1]][1] != temp_name2[[1]][1]){
          # running only if the files are paired
          print(c("You should put the coupled-group in the Inputdata folder!", name1, "is not coupled. Quit!"))
          stop()
        }

        temp_rsn <- as.numeric(gsub("[^0-9\\.]","", temp_name1[[1]][1]))
        temp_ind <- which(NGAW2$Record.Sequence.Number == temp_rsn)
        if(grepl(utils::tail(temp_name1[[1]],1), as.character(NGAW2[temp_ind, 9]))){
          ang1 <- as.numeric(NGAW2[temp_ind, 11])
        } else if(grepl(utils::tail(temp_name1[[1]],1), as.character(NGAW2[temp_ind, 10]))){
          ang1 <- as.numeric(NGAW2[temp_ind, 12])
        } else{
          ang1 <- 0
        }

        assign(paste0(inputpath, '/', temp[2*i-1]), head_dt1 <- scan(paste0(inputpath, '/', temp[2*i-1]), sep = ",", skip = 3, nlines = 1, what = character()))
        dt1 <- as.numeric(gsub("[^0-9\\.]","", head_dt1[2]))
        assign(paste0(inputpath, '/', temp[2*i]), head_dt2 <- scan(paste0(inputpath, '/', temp[2*i]), sep = ",", skip = 3, nlines = 1, what = character()))
        dt2 <- as.numeric(gsub("[^0-9\\.]","", head_dt2[2]))
        if(dt1 != dt2){
          print("Samples per second are different! Quit!")
          stop()
        }
        assign(paste0(inputpath, '/', temp[2*i-1]), data1 <- scan(paste0(inputpath, '/', temp[2*i-1]), sep = "", skip = 4))
        assign(paste0(inputpath, '/', temp[2*i]), data2 <- scan(paste0(inputpath, '/', temp[2*i]), sep = "", skip = 4))

        m <- gregexpr('[0-9]+',temp_name1[[1]][1])
        number <- regmatches(temp_name1[[1]][1],m)
        Record_Sequence_Number[i] <- as.numeric(number[[1]][1])
        EQID[i] <- NGAW2[match(Record_Sequence_Number[i], NGAW2[,1]), 2]

        ## usable frequency: the combined lowest usable frequenct is the higher LUFreq of two components
        lowest_usable_freq <- pmax(NGAW2[match(Record_Sequence_Number[i], NGAW2[,1]), 6],
                                   NGAW2[match(Record_Sequence_Number[i], NGAW2[,1]), 7])
        lowest_freq_combine[i] <- lowest_usable_freq
        ## usable frequency: the combined highest usable frequenct is the lower HUFreq of two components
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

        ## find the new penalty periods based on NGA-West2 usable periods for GMRotI_50 calculation
        penalty_res <- penalty_fun(tmax4penalty_in, tmin4penalty_in, lowest_usable_freq, highest_usable_freq, period_t = period_t)
        tmax4penalty <- penalty_res[1]
        tmin4penalty <- penalty_res[2]
        tmax_index <- penalty_res[3]
        tmin_index <- penalty_res[4]


        ## Process data
        tmp_res <- main_proc(data1 = data1, data2 = data2, period_t = period_t, damping = damping, dt = dt1,
                             fraction = fraction, Interpolation_factor = Interpolation_factor, tmax4penalty = tmax4penalty,
                             tmin4penalty = tmin4penalty, outputdatadir = outputdatadir, outputplotdir = outputplotdir,
                             combine_index = combine_index, ang1 = ang1,
                             lowest_usable_freq = lowest_usable_freq, highest_usable_freq = highest_usable_freq,
                             flname = temp_name1[[1]][1], period_names = period_names)

        pga_combine[i] <- tmp_res$pga_rotxx
        pgv_combine[i] <- tmp_res$pgv_rotxx
        pgd_combine[i] <- tmp_res$pgd_rotxx
        rotd_combine[i,] <- tmp_res$rotxx
      }

    } else{
      print("You should put the couple groups in the Inputdata folder, we detect you have odd number of time series")
      stop()
    }

    ## Create a combined spreadsheet
    name_com <- c('Record_Sequence_Number', 'EQID', 'Lowest_Usable_Freq', 'Highest_Usable_Freq')
    frame.data_com1 <- data.frame(Record_Sequence_Number, EQID, lowest_freq_combine, highest_freq_combine)
    frame.data_com1 <- rbind(name_com, frame.data_com1)
    name_com <- c(paste('PGA_Rot', combine_index, sep = ''), paste('PGV_Rot', combine_index, sep = ''), paste('PGD_Rot', combine_index, sep = ''), period_names)
    value_com <- cbind(pga_combine, pgv_combine, pgd_combine, rotd_combine)
    frame.data_com2 <- t(data.frame(name_com, t(value_com)))
    frame.data_com <- cbind(frame.data_com1, frame.data_com2)
    filename <- paste('summary', ".csv", sep = "")
    write.table(frame.data_com, file = paste0(outputdatadir, '/', filename), col.names = FALSE, row.names = FALSE, sep = ',')
  }


  ## process PEER formatted data
  else if(datatype == 'peer_format'){
    temp = list.files(inputpath)
    if(!length(temp)%%2){
      # running only if the files are paired

      ## define combined spreadsheet variables
      Station <- seq(1, length(temp)/2)
      pga_combine <- seq(1, length(temp)/2)*0
      pgv_combine <- seq(1, length(temp)/2)*0
      pgd_combine <- seq(1, length(temp)/2)*0
      rotd_combine <- matrix(seq(1, (length(temp)/2)*length(period_t)) ,nrow = length(temp)/2, ncol = length(period_t))*0

      for(i in 1:(length(temp)/2)){
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
        assign(paste0(inputpath, '/', temp[2*i-1]), head_dt1 <- scan(paste0(inputpath, '/', temp[2*i-1]), sep = ",", skip = 3, nlines = 1, what = character()))
        dt1 <- as.numeric(gsub("[^0-9\\.]","", head_dt1[2]))
        assign(paste0(inputpath, '/', temp[2*i]), head_dt2 <- scan(paste0(inputpath, '/', temp[2*i]), sep = ",", skip = 3, nlines = 1, what = character()))
        dt2 <- as.numeric(gsub("[^0-9\\.]","", head_dt2[2]))
        if(dt1 != dt2){
          print("Samples per second are different! Quit!")
          stop()
        }
        assign(paste0(inputpath, '/', temp[2*i-1]), data1 <- scan(paste0(inputpath, '/', temp[2*i-1]), sep = "", skip = 4))
        assign(paste0(inputpath, '/', temp[2*i]), data2 <- scan(paste0(inputpath, '/', temp[2*i]), sep = "", skip = 4))


        ## Process data
        tmp_res <- main_proc(data1 = data1, data2 = data2, period_t = period_t, damping = damping, dt = dt1,
                             fraction = fraction, Interpolation_factor = Interpolation_factor, tmax4penalty = tmax4penalty_in,
                             tmin4penalty = tmin4penalty_in, outputdatadir = outputdatadir, outputplotdir = outputplotdir,
                             combine_index = combine_index, ang1 = ang1, flname = temp_name1[[1]][1], period_names = period_names)

        pga_combine[i] <- tmp_res$pga_rotxx
        pgv_combine[i] <- tmp_res$pgv_rotxx
        pgd_combine[i] <- tmp_res$pgd_rotxx
        rotd_combine[i,] <- tmp_res$rotxx
      }

    } else{
      print("You should put the couple groups in the Inputdata folder, we detect you have odd number of time series")
      stop()
    }

    ## Create a combined spreadsheet
    name_com <- c('Station_Name')
    frame.data_com1 <- data.frame(Station)
    frame.data_com1 <- rbind(name_com, frame.data_com1)
    name_com <- c(paste('PGA_Rot', combine_index, sep = ''), paste('PGV_Rot', combine_index, sep = ''), paste('PGD_Rot', combine_index, sep = ''), period_names)
    value_com <- cbind(pga_combine, pgv_combine, pgd_combine, rotd_combine)
    frame.data_com2 <- t(data.frame(name_com, t(value_com)))
    frame.data_com <- cbind(frame.data_com1, frame.data_com2)
    filename <- paste('summary', ".csv", sep = "")
    write.table(frame.data_com, file = paste0(outputdatadir, '/', filename), col.names = FALSE, row.names = FALSE, sep = ',')

  }


  # process NGA-East/Sub formatted data
  else if(datatype == 'nga'){
    temp = list.files(inputpath, pattern = "*.AT2")
    if(!length(temp)%%2){
      # running only if the file are paired

      ## define combined spreadsheet variables
      Record_Sequence_Number <- seq(1, length(temp)/2)
      pga_combine <- seq(1, length(temp)/2)*0
      pgv_combine <- seq(1, length(temp)/2)*0
      pgd_combine <- seq(1, length(temp)/2)*0
      rotd_combine <- matrix(seq(1, (length(temp)/2)*length(period_t)) ,nrow = length(temp)/2, ncol = length(period_t))*0

      for(i in 1:(length(temp)/2)){
        name1 <- basename(temp[2*i-1])
        temp_name1 <- strsplit(name1, split = "[_]")
        name2 <- basename(temp[2*i])
        temp_name2 <- strsplit(name2, split = "[_]")
        if(temp_name1[[1]][1] != temp_name2[[1]][1]){      # judge if there is a coupled group for name1
          print(c("You should put the coupled-group in the Inputdata folder!", name1, "is not coupled. Quit!"))
          stop()
        }

        assign(paste0(inputpath, '/', temp[2*i-1]), head_dt1 <- scan(paste0(inputpath, '/', temp[2*i-1]), sep = ",", skip = 3, nlines = 1, what = character()))
        dt1 <- as.numeric(gsub("[^0-9\\.]","", head_dt1[2]))
        assign(paste0(inputpath, '/', temp[2*i]), head_dt2 <- scan(paste0(inputpath, '/', temp[2*i]), sep = ",", skip = 3, nlines = 1, what = character()))
        dt2 <- as.numeric(gsub("[^0-9\\.]","", head_dt2[2]))
        if(dt1 != dt2){
          print("Samples per second are different! Quit!")
          stop()
        }

        m <- gregexpr('[0-9]+',temp_name1[[1]][1])
        number <- regmatches(temp_name1[[1]][1],m)
        Record_Sequence_Number[i] <- as.numeric(number[[1]][1])
        assign(paste0(inputpath, '/', temp[2*i-1]), data1 <- scan(paste0(inputpath, '/', temp[2*i-1]), sep = "", skip = 4))
        assign(paste0(inputpath, '/', temp[2*i]), data2 <- scan(paste0(inputpath, '/', temp[2*i]), sep = "", skip = 4))


        ## Process data
        tmp_res <- main_proc(data1 = data1, data2 = data2, period_t = period_t, damping = damping, dt = dt1,
                             fraction = fraction, Interpolation_factor = Interpolation_factor, tmax4penalty = tmax4penalty_in,
                             tmin4penalty = tmin4penalty_in, outputdatadir = outputdatadir, outputplotdir = outputplotdir,
                             combine_index = combine_index, ang1 = ang1,
                             flname = temp_name1[[1]][1], period_names = period_names)

        pga_combine[i] <- tmp_res$pga_rotxx
        pgv_combine[i] <- tmp_res$pgv_rotxx
        pgd_combine[i] <- tmp_res$pgd_rotxx
        rotd_combine[i,] <- tmp_res$rotxx
      }

    } else{
      print("You should put the couple groups in the Inputdata folder, we detect you have odd number of time series")
      stop()
    }

    ## Create a combined spreadsheet
    name_com <- c('Record_Sequence_Number')
    frame.data_com1 <- data.frame(Record_Sequence_Number)
    frame.data_com1 <- rbind(name_com, frame.data_com1)
    name_com <- c(paste('PGA_Rot', combine_index, sep = ''), paste('PGV_Rot', combine_index, sep = ''), paste('PGD_Rot', combine_index, sep = ''), period_names)
    value_com <- cbind(pga_combine, pgv_combine, pgd_combine, rotd_combine)
    frame.data_com2 <- t(data.frame(name_com, t(value_com)))
    frame.data_com <- cbind(frame.data_com1, frame.data_com2)
    filename <- paste('summary', ".csv", sep = "")
    write.table(frame.data_com, file = paste0(outputdatadir, '/', filename), col.names = FALSE, row.names = FALSE, sep = ',')

  }

  ## process COSMOS formatted data
  else if(datatype == 'cosmos'){
    ## It only processes the corrected COSMOS data. The accelerations, velocity, and displacements
    ## of two horizontal and one vertical components of ground motions are all in the one file
    temp = list.files(inputpath)
    Station <- seq(1, length(temp))
    pga_combine <- seq(1, length(temp))*0
    pgv_combine <- seq(1, length(temp))*0
    pgd_combine <- seq(1, length(temp))*0
    rotd_combine <- matrix(seq(1, length(temp)*length(period_t)), nrow = length(temp), ncol = length(period_t))*0
    for(i in seq(1,length(temp))){
      name <- basename(temp[i])
      temp_name <- strsplit(name, split = "_")
      station_name <- paste(temp_name[[1]][1], temp_name[[1]][2], sep = "_")
      Station[i] <- station_name
      assign(paste0(inputpath, '/', temp[i]), head_dt1 <- scan(paste0(inputpath, '/', temp[i]), sep = "", skip = 45, nlines = 1, what = character()))
      npts1 <- as.numeric(head_dt1[1])
      num_ind <- which(!is.na(as.numeric(head_dt1)))
      dt1 <- as.numeric(head_dt1[num_ind[2]])
      if(dt1 > 1 ){
        dt1 <- 1/dt1
      }
      len_h1 <- 45 + (ceiling(npts1/8)+1)*3 + 1
      assign(paste0(inputpath, '/', temp[i]), head_dt2 <- scan(paste0(inputpath, '/', temp[i]), sep = "", skip = len_h1 + 45, nlines = 1, what = character()))
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
      line1 <- readLines(paste0(inputpath, '/', temp[i]), n = ceiling(npts1/8)+46)
      data1 <- c()
      for(k in seq(47,ceiling(npts1/8)+46)){
        if(length(data1) == 0){
          data1 <- as.numeric(substring(line1[k], c(1,11,21,31,41,51,61,71),c(10,20,30,40,50,60,70,80)))
        } else{
          data1 <- c(data1,as.numeric(substring(line1[k], c(1,11,21,31,41,51,61,71),c(10,20,30,40,50,60,70,80))))
        }
      }
      data1 <- data1[1:npts1]

      line2 <- readLines(paste0(inputpath, '/', temp[i]), n = len_h1 + ceiling(npts2/8)+46)
      data2 <- c()
      for(k in seq(len_h1 + 47, len_h1+ceiling(npts2/8)+46)){
        if(length(data2) == 0){
          data2 <- as.numeric(substring(line2[k], c(1,11,21,31,41,51,61,71),c(10,20,30,40,50,60,70,80)))
        } else{
          data2 <- c(data2,as.numeric(substring(line2[k], c(1,11,21,31,41,51,61,71),c(10,20,30,40,50,60,70,80))))
        }
      }
      data2 <- data2[1:npts2]

      ## Process data
      tmp_res <- main_proc(data1 = data1, data2 = data2, period_t = period_t, damping = damping, dt = dt1,
                           fraction = fraction, Interpolation_factor = Interpolation_factor, tmax4penalty = tmax4penalty_in,
                           tmin4penalty = tmin4penalty_in, outputdatadir = outputdatadir, outputplotdir = outputplotdir,
                           combine_index = combine_index, ang1 = ang1,
                           flname = station_name, period_names = period_names)

      pga_combine[i] <- tmp_res$pga_rotxx
      pgv_combine[i] <- tmp_res$pgv_rotxx
      pgd_combine[i] <- tmp_res$pgd_rotxx
      rotd_combine[i,] <- tmp_res$rotxx

    }

    ## Create a combined spreadsheet
    name_com <- c('Station_Name')
    frame.data_com1 <- rbind(name_com, t(t(Station)))
    name_com <- c(paste('PGA_Rot', combine_index, sep = ''), paste('PGV_Rot', combine_index, sep = ''), paste('PGD_Rot', combine_index, sep = ''), period_names)
    value_com <- cbind(pga_combine, pgv_combine, pgd_combine, rotd_combine)
    frame.data_com2 <- t(data.frame(name_com, t(value_com)))
    frame.data_com <- cbind(frame.data_com1, frame.data_com2)
    filename <- paste('summary', ".csv", sep = "")
    write.table(frame.data_com, file = paste0(outputdatadir, '/', filename), col.names = FALSE, row.names = FALSE, sep = ',')

  }


  ## process USGS SMC formatted data
  else if(datatype == 'smc'){
    temp = list.files(inputpath)
    if(!length(temp)%%2){
      # running only if the files are paired

      ## define combined spreadsheet variables ##
      Station <- seq(1, length(temp)/2)
      pga_combine <- seq(1, length(temp)/2)*0
      pgv_combine <- seq(1, length(temp)/2)*0
      pgd_combine <- seq(1, length(temp)/2)*0
      rotd_combine <- matrix(seq(1, (length(temp)/2)*length(period_t)) ,nrow = length(temp)/2, ncol = length(period_t))*0

      for(i in 1:(length(temp)/2)){
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
        commlines1 <- scan(paste0(inputpath, '/', temp[2*i-1]), sep = "", skip = 12, nlines = 1, what = character())
        numcomm1 <- as.numeric(commlines1)[8]
        commlines1 <- scan(paste0(inputpath, '/', temp[2*i-1]), sep = "", skip = 13, nlines = 1, what = character())
        npts1 <- as.numeric(commlines1)[1]
        commlines1 <- scan(paste0(inputpath, '/', temp[2*i-1]), sep = "", skip = 17, nlines = 1, what = character())
        dt1 <- 1/as.numeric(commlines1)[2]

        commlines2 <- scan(paste0(inputpath, '/', temp[2*i]), sep = "", skip = 12, nlines = 1, what = character())
        numcomm2 <- as.numeric(commlines2)[8]
        commlines2 <- scan(paste0(inputpath, '/', temp[2*i]), sep = "", skip = 13, nlines = 1, what = character())
        npts2 <- as.numeric(commlines2)[1]
        commlines2 <- scan(paste0(inputpath, '/', temp[2*i]), sep = "", skip = 17, nlines = 1, what = character())
        dt2 <- 1/as.numeric(commlines2)[2]

        if(dt1 != dt2){
          print("Samples per second are different! Quit!")
          stop()
        }

        line1 <- readLines(paste0(inputpath, '/', temp[2*i-1]))
        data1 <- c()
        for(k in seq(28 + numcomm1, length(line1))){
          if(length(data1) == 0){
            data1 <- as.numeric(substring(line1[k], c(1,11,21,31,41,51,61,71),c(10,20,30,40,50,60,70,80)))
          } else{
            data1 <- c(data1,as.numeric(substring(line1[k], c(1,11,21,31,41,51,61,71),c(10,20,30,40,50,60,70,80))))
          }
        }

        line2 <- readLines(paste0(inputpath, '/', temp[2*i]))
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

        ## Process data
        tmp_res <- main_proc(data1 = data1, data2 = data2, period_t = period_t, damping = damping, dt = dt1,
                             fraction = fraction, Interpolation_factor = Interpolation_factor, tmax4penalty = tmax4penalty_in,
                             tmin4penalty = tmin4penalty_in, outputdatadir = outputdatadir, outputplotdir = outputplotdir,
                             combine_index = combine_index, ang1 = ang1,
                             flname = station_name, period_names = period_names)

        pga_combine[i] <- tmp_res$pga_rotxx
        pgv_combine[i] <- tmp_res$pgv_rotxx
        pgd_combine[i] <- tmp_res$pgd_rotxx
        rotd_combine[i,] <- tmp_res$rotxx
      }

    } else{
      print("You should put the couple groups in the Inputdata folder, we detect you have odd number of time series")
      stop()
    }

    ## Create a combined spreadsheet
    name_com <- c('Station_Name')
    frame.data_com1 <- rbind(name_com, t(t(Station)))
    name_com <- c(paste('PGA_Rot', combine_index, sep = ''), paste('PGV_Rot', combine_index, sep = ''), paste('PGD_Rot', combine_index, sep = ''), period_names)
    value_com <- cbind(pga_combine, pgv_combine, pgd_combine, rotd_combine)
    frame.data_com2 <- t(data.frame(name_com, t(value_com)))
    frame.data_com <- cbind(frame.data_com1, frame.data_com2)
    filename <- paste('summary', ".csv", sep = "")
    write.table(frame.data_com, file = paste0(outputdatadir, '/', filename), col.names = FALSE, row.names = FALSE, sep = ',')

  }


  ## process only one column data (the first value is time step)
  else if(datatype == 'timeseries'){
    ## only one column data, the first value of the column is time step
    temp = list.files(inputpath)
    if(!length(temp)%%2){
      # running only if the files are paired

      ## define combined spreadsheet variables ##
      Station <- seq(1, length(temp)/2)
      pga_combine <- seq(1, length(temp)/2)*0
      pgv_combine <- seq(1, length(temp)/2)*0
      pgd_combine <- seq(1, length(temp)/2)*0
      rotd_combine <- matrix(seq(1, (length(temp)/2)*length(period_t)) ,nrow = length(temp)/2, ncol = length(period_t))*0

      for(i in 1:(length(temp)/2)){
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
        dt1 <- as.numeric(readLines(paste0(inputpath, '/', temp[2*i-1]), n = 1))
        dt2 <- as.numeric(readLines(paste0(inputpath, '/', temp[2*i]), n = 1))
        if(dt1 != dt2){
          print("Samples per second are different! Quit!")
          stop()
        }
        data1 <- as.numeric(readLines(paste0(inputpath, '/', temp[2*i-1])))
        data1 <- data1[2:length(data1)]
        data2 <- as.numeric(readLines(paste0(inputpath, '/', temp[2*i])))
        data2 <- data2[2:length(data2)]


        ## Process data
        tmp_res <- main_proc(data1 = data1, data2 = data2, period_t = period_t, damping = damping, dt = dt1,
                             fraction = fraction, Interpolation_factor = Interpolation_factor, tmax4penalty = tmax4penalty_in,
                             tmin4penalty = tmin4penalty_in, outputdatadir = outputdatadir, outputplotdir = outputplotdir,
                             combine_index = combine_index, ang1 = ang1,
                             flname = station_name, period_names = period_names)

        pga_combine[i] <- tmp_res$pga_rotxx
        pgv_combine[i] <- tmp_res$pgv_rotxx
        pgd_combine[i] <- tmp_res$pgd_rotxx
        rotd_combine[i,] <- tmp_res$rotxx

      }

    } else{
      print("You should put the couple groups in the Inputdata folder, we detect you have odd number of time series")
      stop()
    }

    ## Create a combined spreadsheet
    name_com <- c('Station_Name')
    frame.data_com1 <- rbind(name_com, t(t(Station)))
    name_com <- c(paste('PGA_Rot', combine_index, sep = ''), paste('PGV_Rot', combine_index, sep = ''), paste('PGD_Rot', combine_index, sep = ''), period_names)
    value_com <- cbind(pga_combine, pgv_combine, pgd_combine, rotd_combine)
    frame.data_com2 <- t(data.frame(name_com, t(value_com)))
    frame.data_com <- cbind(frame.data_com1, frame.data_com2)
    filename <- paste('summary', ".csv", sep = "")
    write.table(frame.data_com, file = paste0(outputdatadir, '/', filename), col.names = FALSE, row.names = FALSE, sep = ',')

  }

}

## helper functions
penalty_fun <- function(tmax4penalty, tmin4penalty, f_lowest, f_highest, period_t){
  ## penalty calculation function
  # input parameters: 1, tmax4penalty: input maximum penalty period, 10s, only applied if lowest usable frequency from Flatfile is not available
  #                   2, tmin4penalty: input minimum penalty period, 0s, only applied if highest usable frequency from Flatfile is not available
  #                   3, f_lowest: the lowest usable frequency, obtaining from NGA Flatfile
  #                   4, f_highest: the highest usable frequency, obtaining from NGA Flatfile
  #                   5, period_t: the period array in NGA Flatfile
  # output: 1, tmax4penalty: the final maximum period for penalty;
  #         2, tmin4penalty: the final mimimum period for penalty;
  #         3, tmax_index: the index of tmax4penalty in period_t vector
  #         4, tmin_index: the index of tmin4penalty in period_t vector
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

########### main process function ############
main_proc <- function(data1, data2, period_t, damping, dt, fraction, Interpolation_factor, tmax4penalty, tmin4penalty, outputdatadir,
                      outputplotdir, combine_index, ang1, lowest_usable_freq = NA, highest_usable_freq = NA, flname, period_names) {
  ## main process function

  # Applying GMRotD_cal function to get PGA/PGV/PGD and PSA for different periods for each rotated angle,
  # the result is list, 1st is RotD180, 2nd is GMRotD
  results <- GM_RotD_cal(data1 = data1, data2 = data2, period_t = period_t, damping = damping, time_dt = dt, fraction = fraction,
                         Interpolation_factor = Interpolation_factor)
  RotD180 <- results[[1]]
  RotD180 <- matrix(unlist(RotD180), nrow = 180)
  GMRotD <- results[[2]]
  GMRotD <- matrix(unlist(GMRotD), nrow = 90)
  Num_points <- results[[3]]
  Num_points <- c(unlist(Num_points))
  rd_alevel <- results[[4]]

  ## acceleration
  pga_rot00 <- min(results$pga_rot)
  pga_rot50 <- stats::median(results$pga_rot)
  pga_rot100 <- max(results$pga_rot)
  pga_rot00ang <- which.min(results$pga_rot)
  pga_rot100ang <- which.max(results$pga_rot)
  pga_1_ar <- results$pga_rot[180]
  pga_2_ar <- results$pga_rot[90]
  pga_gmrot <- results$pga_gmrot
  pga_gmrot50 <- stats::median(pga_gmrot)
  pga_gmrot100 <- max(pga_gmrot)
  pga_gmrot100ang <- which.max(pga_gmrot)


  ## velocity
  pgv_rot00 <- min(results$pgv_rot)
  pgv_rot50 <- stats::median(results$pgv_rot)
  pgv_rot100 <- max(results$pgv_rot)
  pgv_rot00ang <- which.min(results$pgv_rot)
  pgv_rot100ang <- which.max(results$pgv_rot)
  pgv_1_ar <- results$pgv_rot[180]
  pgv_2_ar <- results$pgv_rot[90]
  pgv_gmrot <- results$pgv_gmrot
  pgv_gmrot50 <- stats::median(pgv_gmrot)
  pgv_gmrot100 <- max(pgv_gmrot)
  pgv_gmrot100ang <- which.max(pgv_gmrot)

  ## displacement
  pgd_rot00 <- min(results$pgd_rot)
  pgd_rot50 <- stats::median(results$pgd_rot)
  pgd_rot100 <- max(results$pgd_rot)
  pgd_rot00ang <- which.min(results$pgd_rot)
  pgd_rot100ang <- which.max(results$pgd_rot)
  pgd_1_ar <- results$pgd_rot[180]
  pgd_2_ar <- results$pgd_rot[90]
  pgd_gmrot <- results$pgd_gmrot
  pgd_gmrot50 <- stats::median(pgd_gmrot)
  pgd_gmrot100 <- max(pgd_gmrot)
  pgd_gmrot100ang <- which.max(pgd_gmrot)

  ## Applying GMRotI_cal function to get PSA for different periods for only one desired rotated angle
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


  ## Applying RotD_cal function to get PSA for different periods for each rotated angle
  RotD00 <- apply(RotD180, 2, min)
  RotD100 <- apply(RotD180, 2, max)
  RotD50 <- apply(RotD180, 2, stats::median)
  RotD_t <- t(RotD180) # transpose of RotD matrix, it is used to get the max index
  RotD100_ang <- max.col(RotD_t)
  RotD00_ang <- max.col(-RotD_t)


  ## Applying PS_cal function to get PSA,PSV
  if(Interpolation_factor == 'auto'){
    interp_factor <- min(ceiling(log(dt/0.001)/log(2)), 3)
    interp_factor <- max(interp_factor, 0)
    interp_factor <- 2^interp_factor
  }else if(Interpolation_factor %% 2 != 0){
    print("The Interpolation factor is not a number of power of 2, please try correct it!")
    stop()
  }else{
    interp_factor <- Interpolation_factor
  }
  dt <- dt/interp_factor # compute new time step
  data1 <- Interpft(data1, interp_factor) # new data after Sinc interpolation
  data2 <- Interpft(data2, interp_factor)
  number <- min(length(data1), length(data2)) # if two data sizes diff, we set both size as the smaller
  data1 <- data1[1:number]
  data2 <- data2[1:number]

  ## compute pseudo spectral acceleration, first row is absolute acceleration of oscillator, second row is psa
  PSA1 <- PS_cal_cpp(data1, period_t, damping, dt, 1)
  PSA2 <- PS_cal_cpp(data2, period_t, damping, dt, 1)
  PSA_gm_ar <- sqrt(PSA1[2,]*PSA2[2,])
  PSA_larger_ar <- seq(1, length(period_t))
  for(j in seq(1, length(PSA_larger_ar))){
    PSA_larger_ar[j] <- max(PSA1[2,j], PSA2[2,j])
  }
  PSA_1 <- PSA1[2,]
  Absolute_acc_1 <- PSA1[1,]
  PSA_2 <- PSA2[2,]
  Absolute_acc_2 <- PSA2[1,]
  npts1 <- length(data1)
  npts2 <- length(data2)


  ## compute real azimuth and set angle in [0,360]
  RotD00_ang <- RotD00_ang + ang1
  RotD100_ang <- RotD100_ang + ang1
  RotD00_ang[which(RotD00_ang > 360)] <- RotD00_ang[which(RotD00_ang > 360)] - 360
  RotD00_ang[which(RotD00_ang < 0)] <- RotD00_ang[which(RotD00_ang < 0)] + 360
  RotD100_ang[which(RotD100_ang > 360)] <- RotD100_ang[which(RotD100_ang > 360)] - 360
  RotD100_ang[which(RotD100_ang < 0)] <- RotD100_ang[which(RotD100_ang < 0)] + 360

  frame.data1 <- data.frame(period_names, PSA_1, Absolute_acc_1, PSA_2, Absolute_acc_2,
                            PSA_gm_ar, PSA_larger_ar, GMRotI50, RotD00, RotD50,
                            RotD100, RotD00_ang, RotD100_ang, Num_points)
  sub_frame.data1 <- t(frame.data1[,-1])
  frame.data1 <- rbind(t(frame.data1[1]), sub_frame.data1)
  Name <- c("npts1","npts2","tmax4penalty","tmin4penalty","damping","Lowest usable freq",
            "Highest usable freq", "GMRotI50angle", "PGA_GMRotI50", "PGV_GMRotI50",
            "PGD_GMRotI50", "PGA_GMRot50", "PGV_GMRot50", "PGA_GMRot100angle",
            "PGA_GMRot100", "PGV_GMRot100angle", "PGV_GMRot100", "PGA_1", "PGA_2", "PGV_1", "PGV_2","PGD_1", "PGD_2",
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

  Value <- c(npts1, npts2, tmax4penalty, tmin4penalty, damping)

  ifelse(is.na(lowest_usable_freq),
         Value <- c(Value, 1/tmax4penalty),
         Value <- c(Value, lowest_usable_freq))

  ifelse(is.na(highest_usable_freq),
         Value <- c(Value, 1/tmin4penalty),
         Value <- c(Value, highest_usable_freq))

  Value <- c(Value,
             GMRotI50_ang, pga_gmrotI50, pgv_gmrotI50, pgd_gmrotI50, pga_gmrot50,
             pgv_gmrot50, pga_gmrot100ang, pga_gmrot100, pgv_gmrot100ang, pgv_gmrot100,
             pga_1_ar, pga_2_ar, pgv_1_ar, pgv_2_ar, pgd_1_ar,
             pgd_2_ar, pga_rot00, pga_rot50, pga_rot100, pga_rot00ang,
             pga_rot100ang, pgv_rot00, pgv_rot50, pgv_rot100, pgv_rot00ang, pgv_rot100ang,
             pgd_rot00, pgd_rot50, pgd_rot100, pgd_rot00ang, pgd_rot100ang)
  frame.data2 <- data.frame(Name, Value)
  filename1 <- paste(flname, "_dep.csv", sep = "")
  filename2 <- paste(flname, "_indep.csv", sep = "")
  write.table(frame.data1, file = paste0(outputdatadir, '/', filename1), col.names = FALSE, sep = ',')
  write.csv(frame.data2, file = paste0(outputdatadir, '/', filename2), row.names = FALSE)


  ## generate plots
  png(filename = paste0(outputplotdir, '/', paste(flname,"_PSA.png", sep = "")), width = 1280, height = 768)
  par(oma=c(2,2,2,2))
  par(mar=c(5,5,4,2) + 0.1)
  plot(c(0.01,20), c(0,max(c(PSA_1,PSA_2))), log = "x", type = "n", xlab = "Period(s)", ylab = "PSA(g)",
       main = "PSA versus Period", cex.main=2, cex.lab=2, cex.axis=2, xaxt = 'n')
  x_range <- c(0.01,20)
  x_range_log10 <- c(ceiling(log10(x_range[1])), floor(log10(x_range[2])))
  # minor ticks
  x_range_minor_log10 <- c(floor(log10(x_range[1])), ceiling(log10(x_range[2])))
  atx_minor <- outer(1:9, 10^(x_range_minor_log10[1]:x_range_minor_log10[2]))
  axis(1, at = atx_minor, labels = FALSE, tcl = par("tcl") * 0.5)
  # major ticks
  atx <- 10^(x_range_log10[1]:x_range_log10[2])
  axis(1, at = atx, labels = atx, cex.axis = 2, cex.lab = 2)
  lines(period_t, PSA_1, lty = 2, lwd = 2, col = "red")
  lines(period_t, PSA_2, lty = 3, lwd = 2, col = "blue")
  legend(5,max(c(PSA_1,PSA_2)), c("PSA_H1","PSA_H2"), lty=c(2,3), lwd=c(2,2),col=c("red","blue"), cex = 2)
  dev.off()
  png(filename = paste0(outputplotdir, '/', paste(flname,"_GMRotI50.png", sep = "")), width = 1280, height = 768)
  par(oma=c(2,2,2,2))
  par(mar=c(5,5,4,2) + 0.1)
  plot(period_t, GMRotI50, log = "x", type = "l", lwd = 2, xlab = "Period(s)", ylab = "GMRotI50(g)", main = "GMRotI50 versus Period",
       cex.lab=2, cex.axis=2, cex.main=2, xaxt = 'n')
  x_range <- c(0.01,20)
  x_range_log10 <- c(ceiling(log10(x_range[1])), floor(log10(x_range[2])))
  # minor ticks
  x_range_minor_log10 <- c(floor(log10(x_range[1])), ceiling(log10(x_range[2])))
  atx_minor <- outer(1:9, 10^(x_range_minor_log10[1]:x_range_minor_log10[2]))
  axis(1, at = atx_minor, labels = FALSE, tcl = par("tcl") * 0.5)
  # major ticks
  atx <- 10^(x_range_log10[1]:x_range_log10[2])
  axis(1, at = atx, labels = atx, cex.axis = 2, cex.lab = 2)
  dev.off()
  png(filename = paste0(outputplotdir, '/', paste(flname,"_RotD.png", sep = "")), width = 1280, height = 768)
  par(oma=c(2,2,2,2))
  par(mar=c(5,5,4,2) + 0.1)
  plot(c(0.01,20), c(0,max(RotD100)), log = "x", type = "n", xlab = "Period(s)", ylab = "RotDxx(g)", main = "RotDxx versus Period", cex.lab=2, cex.axis=2, cex.main=2)
  x_range <- c(0.01,20)
  x_range_log10 <- c(ceiling(log10(x_range[1])), floor(log10(x_range[2])))
  # minor ticks
  x_range_minor_log10 <- c(floor(log10(x_range[1])), ceiling(log10(x_range[2])))
  atx_minor <- outer(1:9, 10^(x_range_minor_log10[1]:x_range_minor_log10[2]))
  axis(1, at = atx_minor, labels = FALSE, tcl = par("tcl") * 0.5)
  # major ticks
  atx <- 10^(x_range_log10[1]:x_range_log10[2])
  axis(1, at = atx, labels = atx, cex.axis = 2, cex.lab = 2)
  lines(period_t, RotD00,  lty = 2, lwd = 2, col = 2)
  lines(period_t, RotD50,  lty = 3, lwd = 2, col = 3)
  lines(period_t, RotD100, lty = 4, lwd = 2, col = 4)
  legend(5, max(RotD100), c("RotD00","RotD50","RotD100"), lty = c(2,3,4),  lwd = c(2,2,2), col = c(2,3,4), cex = 2)
  dev.off()

  ## return
  res <- list()
  if(combine_index == 0){
    res$pga_rotxx <- pga_rot00
    res$pgv_rotxx <- pgv_rot00
    res$pgd_rotxx <- pgd_rot00
    res$rotxx <- RotD00
  } else if(combine_index == 50){
    res$pga_rotxx <- pga_rot50
    res$pgv_rotxx <- pgv_rot50
    res$pgd_rotxx <- pgd_rot50
    res$rotxx <- RotD50
  } else{
    res$pga_rotxx <- pga_rot100
    res$pgv_rotxx <- pgv_rot100
    res$pgd_rotxx <- pgd_rot100
    res$rotxx <- RotD100
  }

  return(res)
}



