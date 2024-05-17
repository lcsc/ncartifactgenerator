#' @name Functions
#' @author
#' Borja Latorre-Garcés \url{https://www.eead.csic.es/home/staffinfo?Id=215}; Soil and Water, EEAD, CSIC \url{https://www.eead.csic.es}
#' Fergus Reig-Gracia \url{http://fergusreig.es}; Environmental Hydrology, Climate and Human Activity Interactions, Geoenvironmental Processes, IPE, CSIC \url{https://www.ipe.csic.es/hidrologia-ambiental}
#' Eduardo Moreno-Lamana \url{https://apuntes.eduardofilo.es}; Environmental Hydrology, Climate and Human Activity Interactions, Geoenvironmental Processes, IPE, CSIC \url{https://www.ipe.csic.es/hidrologia-ambiental}
#' Daniel Vilas-Perulán \url{https://www.eead.csic.es/home/staffinfo?Id=754}; Soil and Water, EEAD, CSIC \url{https://www.eead.csic.es}
#' Manuel Arretxea-Iriarte; Physics of climate and climate change, IGEO, CSIC \url{https://igeo.ucm-csic.es/}
#' @title Auxiliary functions
#' @details
#' \tabular{ll}{
#'   Version: \tab 1.0.0\cr
#'   License: \tab GPL version 3 or newer\cr
#' }
#'
#' @description
#' aux functions

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

library(ncdf4)


#' Read EPSG code from NetCDF file
#'
#' This function reads the EPSG code from a NetCDF file.
#'
#' @param nc The NetCDF file object.
#'
#' @return The EPSG code as a numeric value.
read_epsg <- function(nc) {
  epsg <- 0

  for (att in ncatt_get(nc, 0)) {
    x <- gsub(".*epsg:([0-9]+).*", "\\1", att, ignore.case = TRUE)
    if (x != att) {
      epsg <- x
    }
  }

  return(epsg)
}


#' Return the names of the X and Y dimensions in a NetCDF object
#'
#' This function takes a NetCDF object as input and returns the names of the X and Y dimensions.
#' If the dimensions are named "Y" and "X", the function will return these names along with their positions.
#' If the dimensions are named "latitude" and "longitude", the function will return "latitude" as the Y dimension name,
#' "longitude" as the X dimension name, and their positions.
#' If the dimensions are named "lat" and "lon", the function will return "lat" as the Y dimension name,
#' "lon" as the X dimension name, and their positions.
#'
#' @param nc A NetCDF object.
#'
#' @return A list containing the names and positions of the X and Y dimensions.
returnXYNames <- function(nc) {
  dim <- array(NA, dim = nc$ndims)
  for (i in 1:nc$ndims) {
    dim[i] <- nc$dim[[i]]$name
  }
  if (sum("Y" %in% dim) > 0) {
    return(list(Y = "Y", X = "X", YPos = which("Y" == dim), XPos = which("X" == dim)))
  } else if (sum("y" %in% dim) > 0) {
    return(list(Y = "y", X = "x", YPos = which("y" == dim), XPos = which("x" == dim)))
  } else if (sum("longitude" %in% dim) > 0) {
    return(list(Y = "latitude", X = "longitude", YPos = which("latitude" == dim), XPos = which("longitude" == dim)))
  } else {
    return(list(Y = "lat", X = "lon", YPos = which("lat" == dim), XPos = which("lon" == dim)))
  }
}


#' read_times Function
#'
#' This function reads the time values from a NetCDF file and converts them to Date objects.
#'
#' @param nc The NetCDF file object.
#' @param formatdates Optional parameter to specify the format of the output dates.
#'
#' @return A vector of Date objects representing the time values.
read_times <- function(nc, formatdates) {
  times <- nc$dim[[timePosition(nc)]]$vals
  units <- strsplit(nc$dim[[timePosition(nc)]]$units, " ")[[1]]
  origin <- as.Date(units[which(units == "since") + 1]) # length(units)
  if (length(origin) > 0) {
    times <- as.Date(times, origin = origin)
  }
  if (!missing(formatdates)) {
    times <- format(times, formatdates)
  }
  return(times)
}


#' Read coordinates from a NetCDF file
#'
#' This function reads the longitude and latitude coordinates from a NetCDF file
#' and returns them as a data frame.
#'
#' @param nc The NetCDF file object.
#' @param epsg The EPSG code specifying the coordinate reference system (CRS).
#'
#' @return A data frame containing the longitude and latitude coordinates.
read_coords <- function(nc, epsg) {
  # read spatial dims
  dimNames <- returnXYNames(nc)
  lon <- nc$dim[[dimNames$X]]$vals
  if (lon[1] > lon[length(lon)]) {
    lon <- rev(lon)
  }
  lat <- nc$dim[[dimNames$Y]]$vals
  if (lat[1] > lat[length(lat)]) {
    lat <- rev(lat)
  }

  coords <- expand.grid(lon, lat)
  colnames(coords) <- c("lon", "lat")
  coords <- SpatialPoints(coords, proj4string = CRS(paste0("+init=epsg:", epsg)))
  coords_sf <- st_as_sf(coords)
  coords <- st_transform(coords_sf, CRS("+init=epsg:4326"))
  coords <- st_coordinates(coords)
  return(coords)
}


#' readMinMax Function
#'
#' This function reads the minimum and maximum values of a variable from a NetCDF file.
#'
#' @param nc The NetCDF file object.
#'
#' @return A list containing the minimum and maximum values of the variable.
readMinMax <- function(nc) {
  time_position <- timePosition(nc)
  var_name <- getVarName(nc)
  length_var <- nc$var[[var_name]]$dim[[time_position]]$len
  var_range <- length_var %/% 100
  minimum <- array(NA, length_var)
  maximum <- array(NA, length_var)
  i_var <- 1
  for(i_var in seq(1, length_var, var_range)){
    r_var <- min(var_range, length_var - i_var + 1)
    data <- ncvar_get(nc, nc$var[[var_name]]$name, c(1, 1, i_var), c(-1, -1, r_var))
    minimum[c(i_var:(i_var + r_var - 1))] <- suppressWarnings(apply(data, time_position, min, na.rm = TRUE))
    maximum[c(i_var:(i_var + r_var - 1))] <- suppressWarnings(apply(data, time_position, max, na.rm = TRUE))
  }
  minMax <- list(minimum = minimum, maximum = maximum)
  minMax$minimum[is.na(minMax$min)] <- 0
  minMax$maximum[is.na(minMax$max)] <- 0
  return(minMax)
}

#' merge arrays using aggregation function
#' @param arr1 array 1
#' @param arr2 array 2
#' @param positions positions
#' @param aggregation_func aggregation function to use (e.g., min, max, sum)
#' @return merged array
mergeArrays <- function(arr1, arr2, positions, aggregation_func) {
  if (length(arr1) < length(arr2)) {
    merged_arr <- arr2
  } else {
    merged_arr <- arr1
  }
  i <- 1
  for (i in positions) {
    merged_arr[i] <- aggregation_func(arr1[i], arr2[i], na.rm = TRUE)
  }
  return(merged_arr)
}


#' timePosition
#'
#' This function returns the position of the "time" dimension in a netCDF object.
#'
#' @param nc A netCDF object.
#'
#' @return The position of the "time" dimension.
timePosition <- function(nc) {
  return(grep("time", tolower(names(nc$dim))))
}


#' Get var name of main variable in nc file
#' @param nc open nc file
#'
#' @return var name
getVarName <- function(nc) {
  var_names <- names(nc$var)
  var_name <- var_names[!var_names %in% c("crs")][1]
  return(var_name)
}


#' Convert R array to JavaScript object
#'
#' This function converts an R array into a JavaScript object. It supports different data types such as character, date, and numeric.
#'
#' @param name The name of the JavaScript object.
#' @param value The R array to be converted.
#' @param type The data type of the values in the R array. Default is "character".
#' @param digits The number of digits to round the numeric values. Default is 3.
#' @param value_array Logical value indicating whether the values should be converted into an array in JavaScript. Default is TRUE.
#'
#' @return A string representing the JavaScript object.
arrayRtojs <- function(name, value, type = "character", digits = 3, value_array = TRUE) {
  times <- ""
  for (t in names(value)) {
    if (times != "") {
      times <- paste0(times, ", ")
    }
    if (type == "character") {
      sep <- "'"
      values <- value[[t]]
    } else if (type == "date") {
      sep <- "'"
      values <- as.Date(value[[t]])
    } else {
      sep <- ""
      values <- round(value[[t]], digits = digits)
    }
    if (value_array) {
      times <- paste0(times, "'", t, "': [", sep, paste(values, collapse = paste0(sep, ",", sep)), sep, "]")
    } else {
      times <- paste0(times, "'", t, "': ", sep, values[1], sep)
    }
  }
  times.write <- paste0("var ", name, " = {", times, "};\n")
  return(times.write)
}


#' Converts the data type of a variable or dimension of a netCDF to the data type codes used
#' by the struct library.
#' @param nc_type netCDF data type
#'
#' @return struct data type
get_struct_typecode <- function(nc_type) {
  result <- switch(nc_type,
    "float" = "f",
    "double" = "d",
    "int" = "i"
  )
  return(result)
}


#' Generate artifacts from NetCDF file
#'
#' This function generates artifacts from a NetCDF file, including chunked NetCDF files and corresponding binary files.
#'
#' @param nc_root The root directory of the NetCDF file.
#' @param out_root The root directory where the artifacts will be generated.
#' @param nc_filename The name of the NetCDF file.
#' @param portion The suffix to be added to the generated artifact files.
#' @param var_id The variable ID.
#' @param epsg The EPSG code.
#' @param info_js The information for the JavaScript configuration.
#' @param lon_name The name of the longitude dimension.
#' @param lat_name The name of the latitude dimension.
#' @param write A logical value indicating whether to write the artifacts to disk. Default is FALSE.
#'
#' @return The information for the JavaScript configuration.
#'
#' @export
generate_artifacts <- function(nc_root, out_root,
                               nc_filename, portion, var_id,
                               epsg, info_js, lon_name = NA, lat_name = NA,
                               write = FALSE) {
  print(paste0("Processing: ", var_id, portion, " from file ", nc_filename, portion, ".nc"))

  var_nc_out_folder <- file.path(out_root, "nc")
  dir.create(var_nc_out_folder, showWarnings = FALSE, recursive = TRUE)
  nc_file <- file.path(nc_root, paste0(nc_filename, portion, ".nc"))
  infoNc <- file.info(nc_file)

  # Previous analysis/preparation of NC file
  nc <- nc_open(nc_file, write = TRUE)
  ## read spatial dims
  if (missing(lon_name) || missing(lat_name)) {
    dimNames <- returnXYNames(nc)
    lon_name <- dimNames$X
    lat_name <- dimNames$Y
  }
  var_name <- getVarName(nc)

  print(" Step 0: Order dimensions")
  ## Order dimensions

  var_range <- 100
  lon <- ncvar_get(nc, lon_name)
  lat <- ncvar_get(nc, lat_name)
  if (lon[1] > lon[length(lon)]) {
    for(i_lat in seq(1, length(lat), var_range)){
      r_lat <- min(var_range, length(lat) - i_lat + 1)
      var <- ncvar_get(nc, var_name, c(1, i_lat, 1), c(-1, r_lat, -1))
      var <- var[order(lon), , ]
      ncvar_put(nc, var_name, var, c(1, i_lat, 1), c(-1, r_lat, -1))
      nc_sync(nc)
    }
    lon <- rev(lon)
    ncvar_put(nc, lon_name, lon)
    print(paste0("  ", lon_name, " reversed"))
  }
  if (lat[1] > lat[length(lat)]) {
    for(i_lon in seq(1, length(lon), var_range)){
      r_lon <- min(var_range, length(lon) - i_lon + 1)
      var <- ncvar_get(nc, var_name, c(i_lon, 1, 1), c(r_lon, -1, -1))
      var <- var[, order(lat), ]
      ncvar_put(nc, var_name, var, c(i_lon, 1, 1), c(r_lon, -1, -1))
      nc_sync(nc)
    }
    lat <- rev(lat)
    ncvar_put(nc, lat_name, lat)
    print(paste0("  ", lat_name, " reversed"))
  }
  nc_close(nc)

  # NC t chunks
  print(" Step 1: Chunk_t")
  t_nc_filename <- paste0(var_id, portion, "-t", ".nc")
  t_nc_file <- file.path(var_nc_out_folder, t_nc_filename)
  lon_by <- 100
  lat_by <- 100
  infoT <- file.info(t_nc_file)
  if (is.na(infoT$mtime) || infoT$mtime < infoNc$mtime) {
    write_nc_chunk_t(in_file = nc_file, out_file = t_nc_file, lon_by = lon_by, lat_by = lat_by, lon_name = lon_name, lat_name = lat_name)
  } else {
    print("  Skipped (already newer)")
  }
  # BIN t chunks directory
  print(" Step 2: Bin_t")
  t_bin_filename <- paste0(var_id, portion, "-t", ".bin")
  t_bin_file <- file.path(var_nc_out_folder, t_bin_filename)
  infoT <- file.info(t_bin_file)
  if (is.na(infoT$mtime) || infoT$mtime < infoNc$mtime) {
    write_nc_t_chunk_dir_iter(in_file = t_nc_file, out_file = t_bin_file)
  } else {
    print("  Skipped (already newer)")
  }

  # NC xy chunks
  print(" Step 3: Chunk_xy")
  xy_nc_filename <- paste0(var_id, portion, "-xy", ".nc")
  xy_nc_file <- file.path(var_nc_out_folder, xy_nc_filename)
  time_by <- 100
  infoT <- file.info(xy_nc_file)
  if (is.na(infoT$mtime) || infoT$mtime < infoNc$mtime) {
    write_nc_chunk_xy(in_file = nc_file, out_file = xy_nc_file, time_by = time_by, lon_name = lon_name, lat_name = lat_name)
  } else {
    print("  Skipped (already newer)")
  }

  # BIN xy chunks directory
  print(" Step 4: Bin_xy")
  xy_bin_filename <- paste0(var_id, portion, "-xy", ".bin")
  xy_bin_file <- file.path(var_nc_out_folder, xy_bin_filename)
  infoT <- file.info(xy_bin_file)
  if (is.na(infoT$mtime) || infoT$mtime < infoNc$mtime) {
    write_nc_xy_chunk_dir_iter(in_file = xy_nc_file, out_file = xy_bin_file)
  } else {
    print("  Skipped (already newer)")
  }

  # times.json
  print(" Step 5: times.json")
  nc <- nc_open(nc_file)
  dim_time <- dim(ncvar_get(nc, "time"))
  nc_close(nc)
  args <- list(
    file = nc_file,
    folder = out_root,
    epsg = epsg,
    varName = var_id,
    infoJs = info_js,
    lon_name = lon_name,
    lat_name = lat_name,
    portion = portion,
    varmax = -1,
    varmin = -1
  )
  if (!missing(epsg)) {
    args$epsg <- epsg
  }
  info_js <- do.call(config_web, args)
  return(info_js)
}
