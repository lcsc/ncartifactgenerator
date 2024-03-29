% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/chunk.R
\name{write_nc_xy_chunk_dir}
\alias{write_nc_xy_chunk_dir}
\title{Create the binary chunks directory for the nc oriented to the maps of each date.}
\usage{
write_nc_xy_chunk_dir(in_file, out_file)
}
\arguments{
\item{in_file}{netCDF file with chunking oriented to the maps of each date.}

\item{out_file}{Chunks directory file.}
}
\description{
Create the binary chunks directory for the nc oriented to the maps of each date.
}
\examples{
write_nc_xy_chunk_dir(in_file = "/path/ETo-xy.nc", out_file = "/path/ETo-xy.bin")
}
