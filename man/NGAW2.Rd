% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{NGAW2}
\alias{NGAW2}
\title{Database from NGA-West2 Flatfile}
\format{
A data frame with 21540 rows and 12 variables
\describe{
  \item{Record.Sequence.Number}{record sequence number, identifies a specfic recording}
  \item{EQID}{earthquake ID, identifies a specific earthquake event}
  \item{LP.H1..Hz.}{low pass corner frequency of the frist horizontal component of ground motions}
  \item{LP.H2..Hz.}{low pass corner frequency of the second horizontal component of ground motions}
  \item{Factor}{the ratio for usbale frequency from corner frequency}
  \item{Lowest.Usable.Freq...H1..Hz.}{lowest usable frequency of the frist horizontal component of ground motions}
  \item{Lowest.Usable.Freq...H2..H2.}{lowest usable frequency of the second horizontal component of ground motions}
  \item{Lowest.Usable.Freq...Ave..Component..Hz.}{the average of lowest usable frequency of two horizontal components}
  \item{File.Name..Horizontal.1.}{the file name for the first horziontal component of ground motion in NGA-West2 database}
  \item{File.Name..Horizontal.2.}{the file name for the second horziontal component of ground motion in NGA-West2 database}
  \item{H1.azimth..degrees.}{the azimuth degree of the first horizontal component of ground motion}
  \item{H2.azimth..degrees.}{the azimuth degree of the second horizontal component of ground motion}
}
}
\source{
\url{http://peer.berkeley.edu/ngawest2/databases/}
}
\usage{
NGAW2
}
\description{
A dataset containing 21540 recordings from NGA-West2 Flatfile, which is a reference for metadata.
It helps to link stations and earthquakes, obtain usable frequencies and real azimuth degree of
the sensor.
}
\keyword{datasets}
