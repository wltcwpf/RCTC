#' A Function for Renaming Paired Ground Motions Files and Placed into A New Directory
#'
#' This function changes filenames into standard forms that can be read by IMplot.
#' This operation is not required for the "ngaw2" and "nga" data type (data downloaded from PEER)
#' or if users manually rename data files as instructed in the report.
#' This function transfers data file names into standard names: Station name_x_H1 and Station name_x_H2,
#' while the number x represents sensor number to identify ground motions from different sensors
#' in case of multiple sensors at the same station. For example
#' there are two pairs of ground motions recorded at station LA-BH, then the standardized filenames of
#' the first pair of ground motions are LA-BH_1_H1 and LA-BH_1_H2, the standardized filenames of
#' the second pair of ground motions are LA-BH_2_H1 and LA-BH_2_H2.
#' @param filedir1 The file path (directory + filename) of the first horizontal component of ground motions.
#' @param filedir2 The file path (directory + filename) of the second horizontal component of ground motions. Set it
#' as 'null' if two horizontal components of ground motions are in one file
#' @param stationname The station name, a string
#' @param sn Sensor number, default value is 1. If there is more than one sensor at the station,
#' please enter 2,3,... and so on to identify
#' @param outputdir The directory of output folder where the renamed files will be placed.
#' It can be the \code{Inputdata} folder.
#' @keywords nametransfer
#' @export
#' @examples
#' nametransfer(filedir1 = '/Users/PFW/Desktop/RCTC/20160101231556A1', filedir2 = '/Users/PFW/Desktop/RCTC/20160101231556A2', stationname = 'LA-BH', sn = 1, outputdir = '/Users/PFW/Desktop/RCTC/Inputdata')



nametransfer <- function(filedir1, filedir2, stationname, sn=1, outputdir){
  oldname1 <- basename(filedir1)
  oldname2 <- basename(filedir2)
  if(filedir2 != 'null'){
    newname1 <- paste(stationname, sn, "H1", sep = "_")
    newname2 <- paste(stationname, sn, "H2", sep = "_")
    file.copy(filedir1, outputdir)
    file.rename(file.path(outputdir, oldname1), file.path(outputdir, newname1))
    file.copy(filedir2, outputdir)
    file.rename(file.path(outputdir, oldname2), file.path(outputdir, newname2))
  } else{
    newname1 <- paste(stationname, sn, "H12", sep = "_")
    file.copy(filedir1, outputdir)
    file.rename(file.path(outputdir, oldname1), file.path(outputdir, newname1))
  }
}
