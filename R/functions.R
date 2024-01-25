#' @name Functions
#' @author
#' Borja Latorre Garcés \url{http://eead.csic.es/home/staffinfo?Id=215}; Soil and Water, EEAD, CSIC \url{http://www.eead.csic.es}
#' Fergus Reig Gracia \url{http://fergusreig.es}; Environmental Hydrology, Climate and Human Activity Interactions, Geoenvironmental Processes, IPE, CSIC \url{http://www.ipe.csic.es/hidrologia-ambiental}
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

#' read epsg
#' @param nc nc
#' @return epsg
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

#' raster 3857
#' @param nc nc
#' @return project Raster
raster_3857 <- function(nc, epsg) {
  # read spatial dims
  dimNames <- returnXYNames(nc)
  nrow <- nc$dim[[dimNames$Y]]$len
  ncol <- nc$dim[[dimNames$X]]$len
  lon <- nc$dim[[dimNames$X]]$vals
  if (lon[1] > lon[length(lon)]) {
    lon <- rev(lon)
  }
  lat <- nc$dim[[dimNames$Y]]$vals
  if (lat[1] > lat[length(lat)]) {
    lat <- rev(lat)
  }
  dx <- lon[2] - lon[1]
  dy <- lat[2] - lat[1]

  # create empty raster
  data <- array(numeric(), c(nrow, ncol))
  r <- raster(data[c(nrow:1), ])
  extent(r) <- c(lon[1] - dx / 2, rev(lon)[1] + dx / 2, lat[1] - dy / 2, rev(lat)[1] + dy / 2)
  crs(r) <- crs(paste0("+init=epsg:", epsg))

  # warp to mercator
  r.crs <- projectRaster(from = r, method = "ngb", crs = CRS(paste0("+init=epsg:", "3857")))
  return(r.crs)
}

#' read dims names
#' @param nc nc
#' @return list Y X
returnXYNames <- function(nc) {
  dim <- array(NA, dim = nc$ndims)
  for (i in 1:nc$ndims) {
    dim[i] <- nc$dim[[i]]$name
  }
  if (sum("Y" %in% dim) > 0) {
    return(list(Y = "Y", X = "X", YPos = which("Y" == dim), XPos = which("X" == dim)))
  } else {
    if (sum("longitude" %in% dim) > 0) {
      return(list(Y = "latitude", X = "longitude", YPos = which("latitude" == dim), XPos = which("longitude" == dim)))
    } else {
      # if(sum("lon"%in%dim)>0){
      return(list(Y = "lat", X = "lon", YPos = which("lat" == dim), XPos = which("lon" == dim)))
      # }
    }
  }
}

#' read times
#' @param nc nc
#' @param formatdates formatdates
#' @return times
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

#' read max zoom
#' @param r.crs r.crs
#' @return max zoom
readMaxZoom <- function(r.crs) {
  pixelSize <- min(res(r.crs)) # /1000
  maxzoom <- as.integer(log((20037508.343 * 2) / (256 * pixelSize)) / log(2) + 0.4999)
  return(maxzoom)
}

#' read coords
#' @param nc nc
#' @param epsg epsg
#' @return coords
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

#' read min max
#' @param nc nc
#' @return day min max
readMinMax <- function(nc) {
  var_name <- getVarName(nc)
  data <- ncvar_get(nc, nc$var[[var_name]]$name)
  minMax <- list(min = apply(data, 3, min, na.rm = TRUE), max = apply(data, 3, max, na.rm = TRUE))
  minMax$min[is.na(minMax$min)] <- 0
  minMax$max[is.na(minMax$max)] <- 0
  return(minMax)
}

#' merge min arrays
#' @param min1 min1
#' @param min2 min2
#' @param positions positions
#' @return min
minFusion <- function(min1, min2, positions) {
  if (length(min1) < length(min2)) {
    min.new <- min2
  } else {
    min.new <- min1
  }
  i <- 1
  for (i in positions) {
    min.new[i] <- min(min1[i], min2[i], na.rm = TRUE)
  }
  return(min.new)
}

#' merge max arrays
#' @param max1 max1
#' @param max2 max2
#' @param positions positions
#' @return max
maxFusion <- function(max1, max2, positions) {
  if (length(max1) < length(max2)) {
    max.new <- max2
  } else {
    max.new <- max1
  }
  i <- 1
  for (i in positions) {
    max.new[i] <- max(max1[i], max2[i], na.rm = TRUE)
  }
  return(max.new)
}

#' dim time position
#' @param nc open nc file
#' @return max
timePosition <- function(nc) {
  return(grep("time", tolower(names(nc$dim))))
}


#' Get var name of main variable in nc file
#' @param nc open nc file
#' @return var name
#' @export
getVarName <- function(nc) {
  var_names <- names(nc$var)
  var_name <- var_names[!var_names %in% c("crs")][1]
  return(var_name)
}
