#' @name Chunk
#' @author
#' Borja Latorre-Garcés \url{https://www.eead.csic.es/home/staffinfo?Id=215}; Soil and Water, EEAD, CSIC \url{https://www.eead.csic.es}
#' Fergus Reig-Gracia \url{http://fergusreig.es}; Environmental Hydrology, Climate and Human Activity Interactions, Geoenvironmental Processes, IPE, CSIC \url{https://www.ipe.csic.es/hidrologia-ambiental}
#' Eduardo Moreno-Lamana \url{https://apuntes.eduardofilo.es}; Environmental Hydrology, Climate and Human Activity Interactions, Geoenvironmental Processes, IPE, CSIC \url{https://www.ipe.csic.es/hidrologia-ambiental}
#' Daniel Vilas-Perulán \url{https://www.eead.csic.es/home/staffinfo?Id=754}; Soil and Water, EEAD, CSIC \url{https://www.eead.csic.es}
#' Manuel Arretxea-Iriarte; Physics of climate and climate change, IGEO, CSIC \url{https://igeo.ucm-csic.es/}
#' @title Chunk functions
#' @details
#' \tabular{ll}{
#'   Version: \tab 1.0.0\cr
#'   License: \tab GPL version 3 or newer\cr
#' }
#'
#' @description
#' From a netCDF file, generate two versions with the same data but with different chunk
#' configurations. In one case, favor the retrieval of temporal series of the main
#' variable in each pixel, and in the other, favor the retrieval of planes for each date.

#####################################################################
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/> <http://www.gnu.org/licenses/gpl.txt/>.
#####################################################################

#' @import ncdf4
#' @import hdf5r
#' @import bit64

library(ncdf4)
library(hdf5r)
library(bit64)

# Constants
OFFSET_TYPE_SIZE <- 8
SIZE_TYPE_SIZE <- 4

#' Create the new netCDF file with chunk dimensions that favor obtaining temporal series
#' for each pixel.
#' @param in_file Original netCDF file
#' @param out_file netCDF file with the same information as the original but with new chunk structure.
#' @param lon_by Number of pixels horizontally that will be read as a block during the read/write loop. -1 to read all at once.
#' @param lat_by Number of pixels vertically that will be read as a block during the read/write loop. -1 to read all at once.
#' @param lon_name Name of longitude dimension.
#' @param lat_name Name of latitude dimension.
#' @param signif_digits Number of significant digits to keep in the data.
#' @export
#' @examples
#' write_nc_chunk_t(in_file = "/path/ETo.nc", out_file = "/path/ETo-t.nc", lon_by = 100, lat_by = 100, lon_name = "lon", lat_name = "lat")
write_nc_chunk_t <- function(in_file, out_file, lon_by = -1, lat_by = -1, lon_name = "lon", lat_name = "lat", signif_digits) {
  # Open the original netCDF file
  nc_in_file <- nc_open(in_file)

  var_name <- getVarName(nc_in_file)

  # Reads global attributes
  global_att <- ncatt_get(nc_in_file, 0)

  # Read attributes of dimensions and variable
  lon_longname_att <- ncatt_get(nc_in_file, lon_name, "long_name")
  lon_longname <- if (lon_longname_att$hasatt) lon_longname_att$value else "longitude"
  lon_units_att <- ncatt_get(nc_in_file, lon_name, "units")
  lon_units <- if (lon_units_att$hasatt) lon_units_att$value else "m"

  lat_longname_att <- ncatt_get(nc_in_file, lat_name, "long_name")
  lat_longname <- if (lat_longname_att$hasatt) lat_longname_att$value else "latitude"
  lat_units_att <- ncatt_get(nc_in_file, lat_name, "units")
  lat_units <- if (lat_units_att$hasatt) lat_units_att$value else "m"

  time_longname_att <- ncatt_get(nc_in_file, "time", "long_name")
  time_longname <- if (time_longname_att$hasatt) time_longname_att$value else "time"
  time_units_att <- ncatt_get(nc_in_file, "time", "units")
  time_units <- if (time_units_att$hasatt) time_units_att$value else "days since 1970-01-01"
  time_calendar_att <- ncatt_get(nc_in_file, "time", "calendar")
  time_calendar <- if (time_calendar_att$hasatt) time_calendar_att$value else "gregorian"
  time_unlim <- nc_in_file$dim$time$unlim

  var_longname_att <- ncatt_get(nc_in_file, var_name, "long_name")
  var_longname <- if (var_longname_att$hasatt) var_longname_att$value else NULL
  var_units_att <- ncatt_get(nc_in_file, var_name, "units")
  var_units <- if (var_units_att$hasatt) var_units_att$value else ""
  # On rechunked netCDFs we force missvall to NaN
  var_missval <- NaN # nc_in_file$var[[var_name]]$missval

  # Reads the dimensions of the original file
  lon_data <- ncvar_get(nc_in_file, lon_name)
  lat_data <- ncvar_get(nc_in_file, lat_name)
  time_data <- ncvar_get(nc_in_file, "time")

  # Sizes of dimensions
  lon_num <- length(lon_data)
  lat_num <- length(lat_data)
  time_num <- length(time_data)

  # Checks the size of the read/write batches used for processing large files
  lon_by <- if (lon_by < 1 || lon_by >= lon_num) lon_num else lon_by
  lat_by <- if (lat_by < 1 || lat_by >= lat_num) lat_num else lat_by

  # Define the dimensions for the final file
  lon <- ncdim_def(lon_name, lon_units, lon_data, longname = lon_longname)
  lat <- ncdim_def(lat_name, lat_units, lat_data, longname = lat_longname)
  time <- ncdim_def("time", time_units, time_data,
    unlim = time_unlim,
    longname = time_longname, calendar = time_calendar
  )
  args <- list(
    name = var_name, units = var_units, dim = list(lon, lat, time),
    chunksizes = c(1, 1, time_num), compression = 9
  )
  if (!is.null(var_longname)) {
    args$longname <- var_longname
  }
  if (!is.null(var_missval)) {
    args$missval <- var_missval
  }
  var <- do.call(ncvar_def, args)

  # Final file creation
  nc_out_file <- nc_create(out_file, list(var), force_v4 = TRUE)

  # Select signif function or pass
  if (missing(signif_digits)) {
    adjust_prec <- function(x) x
  } else {
    adjust_prec <- function(x) signif(x, signif_digits)
  }

  if (lon_by == lon_num && lat_by == lat_num) {
    # Read/write data at once
    var_data <- adjust_prec(ncvar_get(nc_in_file, var))
    ncvar_put(nc_out_file, var, var_data)
  } else {
    # Read/write data in batches
    for (x in seq(1, lon_num, by = lon_by)) {
      x_rest <- lon_num - x + 1
      x_count <- if (x_rest >= lon_by) lon_by else x_rest
      for (y in seq(1, lat_num, by = lat_by)) {
        y_rest <- lat_num - y + 1
        y_count <- if (y_rest >= lat_by) lat_by else y_rest
        var_data <- adjust_prec(ncvar_get(nc_in_file, var, start = c(x, y, 1), count = c(x_count, y_count, time_num)))
        ncvar_put(nc_out_file, var, var_data, start = c(x, y, 1), count = c(x_count, y_count, time_num))
      }
    }
  }

  # Writes the global attributes to the final file
  for (name in names(global_att)) {
    ncatt_put(nc_out_file, 0, name, global_att[[name]])
  }

  nc_close(nc_out_file)
  nc_close(nc_in_file)
}


#' Create the new netCDF file with favorable chunk dimensions to obtain plans for each date.
#' @param in_file Original netCDF file
#' @param out_file netCDF file with the same information as the original but with new chunk structure.
#' @param time_by Number of dates that will be read as a block during the read/write loop. -1 to read all at once.
#' @param lon_name Name of longitude dimension.
#' @param lat_name Name of latitude dimension.
#' @param signif_digits Number of significant digits to keep in the data.
#' @export
#' @examples
#' write_nc_chunk_xy(in_file = "/path/ETo.nc", out_file = "/path/ETo-xy.nc", time_by = 100, lon_name = "lon", lat_name = "lat")
write_nc_chunk_xy <- function(in_file, out_file, time_by = -1, lon_name = "lon", lat_name = "lat", signif_digits) {
  # Open the original netCDF file
  nc_in_file <- nc_open(in_file)

  var_name <- getVarName(nc_in_file)

  # Reads global attributes
  global_att <- ncatt_get(nc_in_file, 0)

  # Read attributes of dimensions and variable
  lon_longname_att <- ncatt_get(nc_in_file, lon_name, "long_name")
  lon_longname <- if (lon_longname_att$hasatt) lon_longname_att$value else "longitude"
  lon_units_att <- ncatt_get(nc_in_file, lon_name, "units")
  lon_units <- if (lon_units_att$hasatt) lon_units_att$value else "m"

  lat_longname_att <- ncatt_get(nc_in_file, lat_name, "long_name")
  lat_longname <- if (lat_longname_att$hasatt) lat_longname_att$value else "latitude"
  lat_units_att <- ncatt_get(nc_in_file, lat_name, "units")
  lat_units <- if (lat_units_att$hasatt) lat_units_att$value else "m"

  time_longname_att <- ncatt_get(nc_in_file, "time", "long_name")
  time_longname <- if (time_longname_att$hasatt) time_longname_att$value else "time"
  time_units_att <- ncatt_get(nc_in_file, "time", "units")
  time_units <- if (time_units_att$hasatt) time_units_att$value else "days since 1970-01-01"
  time_calendar_att <- ncatt_get(nc_in_file, "time", "calendar")
  time_calendar <- if (time_calendar_att$hasatt) time_calendar_att$value else "gregorian"
  time_unlim <- nc_in_file$dim$time$unlim

  var_longname_att <- ncatt_get(nc_in_file, var_name, "long_name")
  var_longname <- if (var_longname_att$hasatt) var_longname_att$value else NULL
  var_units_att <- ncatt_get(nc_in_file, var_name, "units")
  var_units <- if (var_units_att$hasatt) var_units_att$value else ""
  # On rechunked netCDFs we force missvall to NaN
  var_missval <- NaN # nc_in_file$var[[var_name]]$missval

  # Reads the dimensions of the original file
  lon_data <- ncvar_get(nc_in_file, lon_name)
  lat_data <- ncvar_get(nc_in_file, lat_name)
  time_data <- ncvar_get(nc_in_file, "time")

  # Sizes of dimensions
  lon_num <- length(lon_data)
  lat_num <- length(lat_data)
  time_num <- length(time_data)

  # Checks the size of the read/write batches used for processing large files
  time_by <- if (time_by < 1 || time_by >= time_num) time_num else time_by

  # Define the dimensions for the final file
  lon <- ncdim_def(lon_name, lon_units, lon_data, longname = lon_longname)
  lat <- ncdim_def(lat_name, lat_units, lat_data, longname = lat_longname)
  time <- ncdim_def("time", time_units, time_data,
    unlim = time_unlim,
    longname = time_longname, calendar = time_calendar
  )
  args <- list(
    name = var_name, units = var_units, dim = list(lon, lat, time),
    chunksizes = c(lon_num, lat_num, 1), compression = 9
  )
  if (!is.null(var_longname)) {
    args$longname <- var_longname
  }
  if (!is.null(var_missval)) {
    args$missval <- var_missval
  }
  var <- do.call(ncvar_def, args)

  # Final file creation
  nc_out_file <- nc_create(out_file, list(var), force_v4 = TRUE)

  # Select signif function or pass
  if (missing(signif_digits)) {
    adjust_prec <- function(x) x
  } else {
    adjust_prec <- function(x) signif(x, signif_digits)
  }

  if (time_by == time_num) {
    # Read/write data at once
    var_data <- adjust_prec(ncvar_get(nc_in_file, var))
    ncvar_put(nc_out_file, var, var_data)
  } else {
    # Read/write data in batches
    for (t in seq(1, time_num, by = time_by)) {
      t_rest <- time_num - t + 1
      t_count <- if (t_rest >= time_by) time_by else t_rest
      var_data <- adjust_prec(ncvar_get(nc_in_file, var_name, start = c(1, 1, t), count = c(lon_num, lat_num, t_count)))
      ncvar_put(nc_out_file, var, var_data, start = c(1, 1, t), count = c(lon_num, lat_num, t_count))
    }
  }

  # Writes the global attributes to the final file
  for (name in names(global_att)) {
    ncatt_put(nc_out_file, 0, name, global_att[[name]])
  }

  nc_close(nc_out_file)
  nc_close(nc_in_file)
}


#' Filename generator for rechunked nc and chunk directories.
#' @param file_name Original netCDF filename
#' @param sufix Optional suffix to be added to the end of the file name before the extension
#' @param ext File extension. If given the empty value, the original extension is maintained
#' @return New filename
#' @export
#' @examples
#' create_nc_name(file_name = "ETo.nc", sufix = "-t", ext = ".bin")
create_nc_name <- function(file_name, sufix = "-t", ext = "") {
  pos <- unlist(gregexpr(".nc", file_name))
  ext_pos <- pos[length(pos)]
  if (ext == "") {
    ext <- substr(file_name, ext_pos, nchar(file_name))
  }
  return(paste0(substr(file_name, 1, ext_pos - 1), sufix, ext))
}


#' Create the binary chunks directory for the nc oriented to the time series of each pixel.
#' @param in_file netCDF file with chunking oriented to the time series of each pixel.
#' @param out_file Chunks directory file.
#' @export
#' @examples
#' write_nc_t_chunk_dir(in_file = "/path/ETo-t.nc", out_file = "/path/ETo-t.bin")
write_nc_t_chunk_dir <- function(in_file, out_file) {
  nc_in_file <- nc_open(in_file)
  var_name <- getVarName(nc_in_file)
  nc_close(nc_in_file)

  # Open the netCDF file with hdf5r
  nc_in_file <- h5file(in_file, mode = "r")

  # Open a binary file in write mode
  bin_out_file <- file(out_file, "wb")

  lon_num <- nc_in_file[[var_name]]$dims[[1]]
  lat_num <- nc_in_file[[var_name]]$dims[[2]]
  time_num <- nc_in_file[[var_name]]$dims[[3]]

  for (y in seq(1, lat_num)) {
    for (x in seq(1, lon_num)) {
      chunk_info <- nc_in_file[[var_name]]$get_chunk_info_by_coord(c(0, y - 1, x - 1))
      # Write the number pairs to the binary file
      writeBin(as.vector(as.integer64(chunk_info$addr)), bin_out_file, endian = "little")
      writeBin(as.integer(chunk_info$size), bin_out_file, size = SIZE_TYPE_SIZE, endian = "little")
    }
  }

  # Close files
  close(bin_out_file)
  nc_in_file$close()
}


#' Create the binary chunks directory for the nc oriented to the time series of each pixel
#' using the new H5Dchunk_iter function (https://docs.hdfgroup.org/hdf5/v1_14/group___h5_d.html#title6).
#' @param in_file netCDF file with chunking oriented to the time series of each pixel.
#' @param out_file Chunks directory file.
#' @export
#' @examples
#' write_nc_t_chunk_dir_iter(in_file = "/path/ETo-t.nc", out_file = "/path/ETo-t.bin")
write_nc_t_chunk_dir_iter <- function(in_file, out_file) {
  nc_in_file <- nc_open(in_file)
  var_name <- getVarName(nc_in_file)
  nc_close(nc_in_file)

  # Open the netCDF file with hdf5r
  nc_in_file <- h5file(in_file, mode = "r")

  # Open a binary file in write mode
  bin_out_file <- file(out_file, "wb")

  lon_num <- nc_in_file[[var_name]]$dims[[1]]
  lat_num <- nc_in_file[[var_name]]$dims[[2]]
  time_num <- nc_in_file[[var_name]]$dims[[3]]

  nc_in_file[[var_name]]$chunk_iter(function(chunk_info) {
    lat_index <- chunk_info$offset[[2]]
    lon_index <- chunk_info$offset[[3]]
    dir_pos <- (OFFSET_TYPE_SIZE + SIZE_TYPE_SIZE) * (lon_index + lat_index * lon_num)
    seek(bin_out_file, dir_pos)
    writeBin(as.vector(as.integer64(chunk_info$addr)), bin_out_file, endian = "little")
    writeBin(as.integer(chunk_info$size), bin_out_file, size = SIZE_TYPE_SIZE, endian = "little")
  })

  # Close files
  close(bin_out_file)
  nc_in_file$close()
}


#' Create the binary chunks directory for the nc oriented to the maps of each date.
#' @param in_file netCDF file with chunking oriented to the maps of each date.
#' @param out_file Chunks directory file.
#' @export
#' @examples
#' write_nc_xy_chunk_dir(in_file = "/path/ETo-xy.nc", out_file = "/path/ETo-xy.bin")
write_nc_xy_chunk_dir <- function(in_file, out_file) {
  nc_in_file <- nc_open(in_file)
  var_name <- getVarName(nc_in_file)
  nc_close(nc_in_file)

  # Open the netCDF file with hdf5r
  nc_in_file <- h5file(in_file, mode = "r")

  # Open a binary file in write mode
  bin_out_file <- file(out_file, "wb")

  lon_num <- nc_in_file[[var_name]]$dims[[1]]
  lat_num <- nc_in_file[[var_name]]$dims[[2]]
  time_num <- nc_in_file[[var_name]]$dims[[3]]

  for (t in seq(1, time_num)) {
    chunk_info <- nc_in_file[[var_name]]$get_chunk_info_by_coord(c(t - 1, 0, 0))
    # Write the number pairs to the binary file
    writeBin(as.vector(as.integer64(chunk_info$addr)), bin_out_file, endian = "little")
    writeBin(as.integer(chunk_info$size), bin_out_file, size = SIZE_TYPE_SIZE, endian = "little")
  }

  # Close files
  close(bin_out_file)
  nc_in_file$close()
}


#' Create the binary chunks directory for the nc oriented to the maps of each date
#' using the new H5Dchunk_iter function (https://docs.hdfgroup.org/hdf5/v1_14/group___h5_d.html#title6).
#' @param in_file netCDF file with chunking oriented to the maps of each date.
#' @param out_file Chunks directory file.
#' @export
#' @examples
#' write_nc_xy_chunk_dir_iter(in_file = "/path/ETo-xy.nc", out_file = "/path/ETo-xy.bin")
write_nc_xy_chunk_dir_iter <- function(in_file, out_file) {
  nc_in_file <- nc_open(in_file)
  var_name <- getVarName(nc_in_file)
  nc_close(nc_in_file)

  # Open the netCDF file with hdf5r
  nc_in_file <- h5file(in_file, mode = "r")

  # Open a binary file in write mode
  bin_out_file <- file(out_file, "wb")

  nc_in_file[[var_name]]$chunk_iter(function(chunk_info) {
    time_index <- chunk_info$offset[[1]]
    dir_pos <- (OFFSET_TYPE_SIZE + SIZE_TYPE_SIZE) * time_index
    seek(bin_out_file, dir_pos)
    writeBin(as.vector(as.integer64(chunk_info$addr)), bin_out_file, endian = "little")
    writeBin(as.integer(chunk_info$size), bin_out_file, size = SIZE_TYPE_SIZE, endian = "little")
  })

  # Close files
  close(bin_out_file)
  nc_in_file$close()
}


#' Union of peninsular and Canary Islands netCDF files.
#' @param can_filename netCDF file with data for the Canary Islands.
#' @param pen_filename netCDF file with data for the Iberian Peninsula.
#' @param fusion_filename netCDF file with the union of the data from the two files.
#' @param nc_chunk_num Number of chunks to use in the new netCDF file in each dimension.
#' @param signif_digits Number of significant digits to keep in the data.
#' @param lon_name Name of longitude dimension.
#' @param lat_name Name of latitude dimension.
#' @param time_name Name of time dimension.
#' @param grid_size Size of the grid in degrees.
#' @export
#' @examples
#' fusion_pen_can(can_filename = "/path/ETo_can.nc", pen_filename = "/path/ETo_pen.nc", fusion_filename = "/path/ETo_all.nc", nc_chunk_num = 16, signif_digits = 4)
fusion_pen_can <- function(can_filename,
                           pen_filename,
                           fusion_filename,
                           nc_chunk_num = 16,
                           signif_digits,
                           lon_name = "lon",
                           lat_name = "lat",
                           time_name = "time",
                           grid_size) {
  # Select signif function or pass
  if (missing(signif_digits)) {
    adjust_prec <- function(x) x
  } else {
    adjust_prec <- function(x) signif(x, signif_digits)
  }

  # Open the original netCDF files
  can <- nc_open(can_filename, write = FALSE)
  pen <- nc_open(pen_filename, write = FALSE)

  min_lon <- min(ncvar_get(pen, lon_name), ncvar_get(can, lon_name))
  max_lon <- max(ncvar_get(pen, lon_name), ncvar_get(can, lon_name))
  min_lat <- min(ncvar_get(pen, lat_name), ncvar_get(can, lat_name))
  max_lat <- max(ncvar_get(pen, lat_name), ncvar_get(can, lat_name))

  # Define dimensions
  if (missing(grid_size)) {
    grid_size <- min(abs(diff(ncvar_get(pen, lon_name))), abs(diff(ncvar_get(pen, lat_name))))
  }

  # Define dimensions (modified to exclude the first row and column of pixels)
  dimLon <- ncdf4::ncdim_def(
    name = lon_name,
    units = "degrees_east",
    vals = seq(min_lon,
      max_lon,
      by = grid_size
    ),
    longname = "longitude"
  )

  dimLat <- ncdf4::ncdim_def(
    name = lat_name,
    units = "degrees_north",
    vals = seq(min_lat,
      max_lat,
      by = grid_size
    ),
    longname = "latitude"
  )

  time_longname_att <- ncatt_get(pen, time_name, "long_name")
  time_longname <- if (time_longname_att$hasatt) time_longname_att$value else "time"
  time_units_att <- ncatt_get(pen, time_name, "units")
  time_units <- if (time_units_att$hasatt) time_units_att$value else "days since 1961-01-01"
  time_calendar_att <- ncatt_get(pen, time_name, "calendar")
  time_calendar <- if (time_calendar_att$hasatt) time_calendar_att$value else "gregorian"
  time_unlim <- pen$dim$time$unlim
  time_data <- ncvar_get(pen, time_name)
  dimTime <- ncdf4::ncdim_def(
    name = time_name,
    units = time_units,
    vals = time_data,
    unlim = time_unlim,
    calendar = time_calendar,
    longname = time_longname
  )

  # Define the CRS variable
  varCRS <- ncdf4::ncvar_def(
    name = "crs",
    units = "",
    dim = list(),
    longname = "CRS definition",
    prec = "integer"
  )

  # Calc chunk sizes
  chunk_lon <- ceiling(dimLon$len / nc_chunk_num)
  chunk_lat <- ceiling(dimLat$len / nc_chunk_num)
  if (nc_chunk_num > dimTime$len) nc_chunk_num <- 1
  chunk_time <- ceiling(dimTime$len / nc_chunk_num)

  # Define the variable
  var_name <- getVarName(pen)
  var_longname_att <- ncatt_get(pen, var_name, "long_name")
  var_longname <- if (var_longname_att$hasatt) var_longname_att$value else NULL
  var_units_att <- ncatt_get(pen, var_name, "units")
  var_units <- if (var_units_att$hasatt) var_units_att$value else ""
  # On fusioned netCDF we force missvall to NaN
  var_missval <- NA # pen$var[[var_name]]$missval
  args <- list(
    name = var_name, units = var_units, dim = list(dimLon, dimLat, dimTime),
    chunksizes = c(chunk_lon, chunk_lat, chunk_time),
    prec = "float", compression = 9
  )
  if (!is.null(var_longname)) {
    args$longname <- var_longname
  }
  if (!is.null(var_missval)) {
    args$missval <- var_missval
  }
  var <- do.call(ncvar_def, args)
  # Create a new netCDF file and add the variables
  nc <- ncdf4::nc_create(fusion_filename, vars = list(var, varCRS), force_v4 = TRUE, verbose = TRUE)

  # Link the variable to the projection variable using the
  # grid_mapping attribute
  ncdf4::ncatt_put(nc, var_name, "grid_mapping", "crs")

  # CRS attributes
  ncdf4::ncatt_put(nc, "crs", "grid_mapping_name", "latitude_longitude")
  ncdf4::ncatt_put(nc, "crs", "longitude_of_prime_meridian", 0.0)
  ncdf4::ncatt_put(nc, "crs", "semi_major_axis", 6378137.0)
  ncdf4::ncatt_put(nc, "crs", "inverse_flattening", 298.257223563)
  ncdf4::ncatt_put(nc, "crs", "crs_wkt", 'GEOGCRS["WGS 84", ENSEMBLE["World Geodetic System 1984 ensemble",
                 MEMBER["World Geodetic System 1984 (Transit)"],
                 MEMBER["World Geodetic System 1984 (G730)"],
                 MEMBER["World Geodetic System 1984 (G873)"],
                 MEMBER["World Geodetic System 1984 (G1150)"],
                 MEMBER["World Geodetic System 1984 (G1674)"],
                 MEMBER["World Geodetic System 1984 (G1762)"],
                 MEMBER["World Geodetic System 1984 (G2139)"],
                 ELLIPSOID["WGS 84",6378137,298.257223563,
                           LENGTHUNIT["metre",1]],
                 ENSEMBLEACCURACY[2.0]],
        PRIMEM["Greenwich",0,
               ANGLEUNIT["degree",0.0174532925199433]],
        CS[ellipsoidal,2],
        AXIS["geodetic latitude (Lat)",north,
             ORDER[1],
             ANGLEUNIT["degree",0.0174532925199433]],
        AXIS["geodetic longitude (Lon)",east,
             ORDER[2],
             ANGLEUNIT["degree",0.0174532925199433]],
        USAGE[
          SCOPE["Horizontal component of 3D system."],
          AREA["World."],
          BBOX[-90,-180,90,180]],
        ID["EPSG",4326]]')

  # Add longitude attributes
  ncdf4::ncatt_put(nc, lon_name, "long_name", "longitude")
  ncdf4::ncatt_put(nc, lon_name, "standard_name", "longitude")
  ncdf4::ncatt_put(nc, lon_name, "axis", "X")
  ncdf4::ncatt_put(nc, lon_name, "comment", "Longitude geographical coordinates, WGS84 projection")
  ncdf4::ncatt_put(nc, lon_name, "reference_datum", "geographical coordinates, WGS84 projection")

  # Add latitude attributes
  ncdf4::ncatt_put(nc, lat_name, "long_name", "latitude")
  ncdf4::ncatt_put(nc, lat_name, "standard_name", "latitude")
  ncdf4::ncatt_put(nc, lat_name, "axis", "Y")
  ncdf4::ncatt_put(nc, lat_name, "comment", "Latitude geographical coordinates, WGS84 projection")
  ncdf4::ncatt_put(nc, lat_name, "reference_datum", "geographical coordinates, WGS84 projection")

  # Add time attributes
  ncdf4::ncatt_put(nc, time_name, "axis", "T")

  # Add attributes
  ncdf4::ncatt_put(nc, var_name, "standard_name", var_longname)

  # Add global attributes
  ncdf4::ncatt_put(nc, 0, "Conventions", "CF-1.8")

  ######################################

  # Round the longitude and latitude values to 4 decimal places
  nc_lon_rounded <- round(ncvar_get(nc, lon_name), 4)
  nc_lat_rounded <- round(ncvar_get(nc, lat_name), 4)
  can_lon_rounded <- round(ncvar_get(can, lon_name), 4)
  can_lat_rounded <- round(ncvar_get(can, lat_name), 4)
  pen_lon_rounded <- round(ncvar_get(pen, lon_name), 4)
  pen_lat_rounded <- round(ncvar_get(pen, lat_name), 4)

  # Determine the corresponding extent in the new NetCDF file
  can_lon_start <- which(nc_lon_rounded == min(can_lon_rounded))
  can_lat_start <- which(nc_lat_rounded == min(can_lat_rounded))
  pen_lon_start <- which(nc_lon_rounded == min(pen_lon_rounded))
  pen_lat_start <- which(nc_lat_rounded == min(pen_lat_rounded))
  print("Start of fusion")
  # Read all data for the Canary Islands from the existing file
  print("Canary read")
  var_data_can <- adjust_prec(ncvar_get(can, var_name))
  if (length(dim(var_data_can)) == 2) {
    var_data_can <- array(var_data_can, dim = c(dim(var_data_can)[1], dim(var_data_can)[2], 1))
  }
  # Write all data for the Canary Islands to the new file
  ncvar_put(nc, var, var_data_can,
    start = c(can_lon_start, can_lat_start, 1),
    count = dim(var_data_can)
  )
  print("Canary written")

  # Read all data for the Iberian Peninsula from the existing file in batches
  print("Peninsula read")
  # Reads the dimensions of the original file
  lon_data <- ncvar_get(pen, lon_name)
  lat_data <- ncvar_get(pen, lat_name)
  time_data <- ncvar_get(pen, time_name)
  # Sizes of dimensions
  lon_num <- length(lon_data)
  lat_num <- length(lat_data)
  time_num <- length(time_data)

  # Read/write data in batches
  for (t in seq(1, time_num, by = chunk_time)) {
    t_rest <- time_num - t + 1
    t_count <- if (t_rest >= chunk_time) chunk_time else t_rest
    var_data <- adjust_prec(ncvar_get(pen, var_name, start = c(1, 1, t), count = c(lon_num, lat_num, t_count)))
    if (length(dim(var_data)) == 2) {
      var_data <- array(var_data, dim = c(dim(var_data)[1], dim(var_data)[2], 1))
    }
    ncvar_put(nc, var, var_data, start = c(pen_lon_start, pen_lat_start, t), count = dim(var_data))
  }
  print("Peninsula written")

  # Close all the NetCDF files to save the changes
  nc_close(nc)
  nc_close(can)
  nc_close(pen)
}


## Usage

# nc_route = "../viewer/nc"
# ncFile = "ETo.nc"
# file = file.path(nc_route, ncFile)
# t_file = file.path(nc_route, create_nc_name(ncFile))
# lon_by = 100
# lat_by = 100
# write_nc_chunk_t(in_file=file, out_file=t_file, lon_by=lon_by, lat_by=lat_by)


# nc_route = "../viewer/nc"
# ncFile = "ETo.nc"
# file = file.path(nc_route, ncFile)
# xy_file = file.path(nc_route, create_nc_name(ncFile, sufix="-xy"))
# time_by = 100
# write_nc_chunk_xy(in_file=file, out_file=xy_file, time_by=time_by)


# nc_route = "../viewer/nc"
# ncFile = "ETo.nc"
# file = file.path(nc_route, ncFile)
# ncEnv_route = "../viewer"
# ncEnvFile = "ncEnv.js"
# envFile = file.path(ncEnv_route, ncEnvFile)
# write_nc_env(in_file=file, out_file=envFile)


# nc_route = "/home/docker/workdir/proto_eto/viewer/nc"
# ncFile = "ETo-t.nc"
# file = file.path(nc_route, ncFile)
# bin_file = file.path(nc_route, create_nc_name(ncFile, ext="bin"))
# write_nc_t_chunk_dir(in_file = file, out_file = bin_file)
#  or
# write_nc_t_chunk_dir_iter(in_file = file, out_file = bin_file)


# nc_route = "/home/docker/workdir/proto_eto/viewer/nc"
# ncFile = "ETo-xy.nc"
# file = file.path(nc_route, ncFile)
# bin_file = file.path(nc_route, create_nc_name(ncFile, ext="bin"))
# write_nc_xy_chunk_dir(in_file = file, out_file = bin_file)
#  or
# write_nc_xy_chunk_dir_iter(in_file = file, out_file = bin_file)
