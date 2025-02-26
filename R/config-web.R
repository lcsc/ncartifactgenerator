#' @name Config-web
#' @author
#' Borja Latorre-Garcés \url{https://www.eead.csic.es/home/staffinfo?Id=215}; Soil and Water, EEAD, CSIC \url{https://www.eead.csic.es}
#' Fergus Reig-Gracia \url{http://fergusreig.es}; Environmental Hydrology, Climate and Human Activity Interactions, Geoenvironmental Processes, IPE, CSIC \url{https://www.ipe.csic.es/hidrologia-ambiental}
#' Eduardo Moreno-Lamana \url{https://apuntes.eduardofilo.es}; Environmental Hydrology, Climate and Human Activity Interactions, Geoenvironmental Processes, IPE, CSIC \url{https://www.ipe.csic.es/hidrologia-ambiental}
#' Daniel Vilas-Perulán \url{https://www.eead.csic.es/home/staffinfo?Id=754}; Soil and Water, EEAD, CSIC \url{https://www.eead.csic.es}
#' Manuel Arretxea-Iriarte; Physics of climate and climate change, IGEO, CSIC \url{https://igeo.ucm-csic.es/}
#' @title Config-web functions
#' @details
#' \tabular{ll}{
#'   Version: \tab 1.0.0\cr
#'   License: \tab GPL version 3 or newer\cr
#' }
#'
#' @description
#' write web configuration

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

#' @import js
#' @import sf
#' @import ncdf4
#' @import raster
#' @import jsonlite

library(js)
library(sf)
library(ncdf4)
library(raster)
library(jsonlite)


#' Web configuration
#'
#' This function generates the web configuration for a specific netCDF file.
#'
#' @param file Path to the netCDF file.
#' @param folder Path to the output folder.
#' @param epsg EPSG code for spatial projection.
#' @param formatdates Date format for time variables.
#' @param varmin Minimum value for the main variable.
#' @param varmax Maximum value for the main variable.
#' @param varName Name of the netCDF file without the .nc extension.
#' @param infoJs List with additional information in JavaScript format.
#' @param lon_name Name of the longitude dimension.
#' @param lat_name Name of the latitude dimension.
#' @param time_name Name of the time dimension.
#' @param time_by Number of time steps that will be read as a block during the read/write loop. -1 to read all at once. Default is 100.
#' @param portion Portion of the netCDF file to process.
#'
#' @return List with the generated web configuration.
#'
#' @export
config_web <- function(file, folder, epsg, formatdates, varmin, varmax, varName, infoJs = NA, lon_name, lat_name, time_name = "time", time_by = 100, portion) {
  if (missing(infoJs) || sum(!is.na(infoJs)) == 0) {
    infoJs <- list(
      times = list(),
      latMin = list(), latMax = list(), latNum = list(),
      lonMin = list(), lonMax = list(), lonNum = list(),
      minVal = list(), maxVal = list(),
      varType = "f", espg = epsg,
      portions = list()
    )
  }

  # CLARIFICATION!
  # var_name is the name of the primary variable within the ncfile
  # varName is the name of the nc file without the .nc extension
  if (!missing(file)) {
    # open nc
    nc <- nc_open(file)
    var_name <- getVarName(nc)
  }

  if (missing(varName)) {
    if (!missing(file)) {
      # nc name
      varName <- basename(gsub(".nc", "", file))
    }
    if (missing(file) | write) {
      varName <- "NaN"
    }
  }

  # read epsg
  if (missing(epsg)) {
    epsg <- read_epsg(nc)
    if (epsg == 0) {
      cat("Error: EPSG not found\n")
      return()
    }
  }
  infoJs$epsg <- epsg

  dir.create(folder, showWarnings = FALSE, recursive = TRUE)

  # Calculates the minimum and maximum of the netCDF dimensions corresponding to the spatial coordinates
  lon <- nc$dim[[lon_name]]$vals
  lat <- nc$dim[[lat_name]]$vals
  infoJs$lonMin[[paste0(varName, portion)]] <- min(lon)
  infoJs$lonMax[[paste0(varName, portion)]] <- max(lon)
  infoJs$lonNum[[paste0(varName, portion)]] <- length(lon)
  infoJs$latMin[[paste0(varName, portion)]] <- min(lat)
  infoJs$latMax[[paste0(varName, portion)]] <- max(lat)
  infoJs$latNum[[paste0(varName, portion)]] <- length(lat)
  
  varMinMax <- NULL;
  if (missing(varmin) | missing(varmax)) {
    print ("Reading Min Max...")
    varMinMax <- readMinMax(nc, time_by)
  } else {
    varMinMax <- list(minimum = varmin, maximum = varmax)
  }
  infoJs$minVal[[paste0(varName)]] <- min(varMinMax$minimum)
  infoJs$maxVal[[paste0(varName)]] <- max(varMinMax$maximum)
  
  
  infoJs$portions[[varName]] <- I(c(infoJs$portions[[varName]], portion))

  infoJs$varType <- get_struct_typecode(nc$var[[var_name]]$prec)

  # Available times
  if (!missing(formatdates)) {
    times.write <- read_times(nc, formatdates)
  } else {
    times.write <- read_times(nc)
  }

  if (missing(varmin) | missing(varmax)) {
    varmin <- varMinMax$minimum
    varmax <- varMinMax$maximum
  } else {
    if (!missing(file)) {
      varmin <- array(varmin, dim = length(nc$dim[[timePosition(nc)]]$vals))
      varmax <- array(varmax, dim = length(nc$dim[[timePosition(nc)]]$vals))
    }
  }

  # Delete Inf
  aux <- varmin == Inf | varmin == -Inf | varmax == Inf | varmax == -Inf
  varmin[aux] <- 0
  varmax[aux] <- 100

  positions <- 1:length(varmin)
  if (!is.null(infoJs$times[[varName]])) {
    if (length(times.write) == sum(times.write %in% infoJs$times[[varName]])) {
      positions <- match(times.write, infoJs$times[[varName]])
      times.write <- infoJs$times[[varName]]
    }
  }

  # Checks if there is varmin and varmax information for varName in infoJs.
  # If there is, the varmin and varmax values are merged with the infoJs values
  # using the mergeArrays() function.
  if (!is.null(infoJs$varmin[[varName]]) & !is.null(infoJs$varmax[[varName]])) {
    varmin <- mergeArrays(arr1 = varmin, arr2 = infoJs$varmin[[varName]], positions = positions, aggregation_func = min)
    varmax <- mergeArrays(arr1 = varmax, arr2 = infoJs$varmax[[varName]], positions = positions, aggregation_func = max)
  }

  infoJs$varmin[[varName]] <- varmin
  infoJs$varmax[[varName]] <- varmax
  infoJs$times[[varName]] <- times.write

  return(infoJs)
}


#' Write JSON string with configuration data
#'
#' This function generates a JSON string with configuration data similar to the JavaScript variables in writeJs.
#'
#' @param infoJs A list containing the configuration data, including times, variable minimum and maximum values, latitude and longitude ranges, and EPSG code.
#' @param varTitle The title of the variable.
#' @param legendTitle The title of the legend (optional, default is "Legend").
#' @param offsetType The type of offset (optional, default is "Q").
#' @param sizeType The type of size (optional, default is "I").
#' @param minify Whether the JSON string should be minified (optional, default is TRUE).
#'
#' @return The generated JSON string.
#'
#' @export
writeJson <- function(folder, infoJs, varTitle, legendTitle = "Legend", offsetType = "Q", sizeType = "I", minify = TRUE) {
  file <- file.path(folder, "times.json")

  json <- list()

  # Center of all the maps
  lonM <- mean(c(min(unlist(infoJs$lonMin)), max(unlist(infoJs$lonMax))))
  latM <- mean(c(min(unlist(infoJs$latMin)), max(unlist(infoJs$latMax))))
  center_point <- st_sfc(st_point(c(lonM, latM)))

  # Transform the point to the geographic coordinate system (EPSG:4326)
  st_crs(center_point) <- as.numeric(infoJs$epsg)
  point_transformed <- st_transform(center_point, crs = 4326)
  center_transformed <- st_coordinates(point_transformed)

  json$center <- list(lat = center_transformed[1, ][["Y"]], lng = center_transformed[1, ][["X"]])
  json$times <- infoJs$times
  json$varMin <- infoJs$varmin
  json$varMax <- infoJs$varmax
  json$varTitle <- varTitle
  if (length(legendTitle) > 1) {
    json$legendTitle <- legendTitle
  } else {
    json$legendTitle <- list("NaN" = list(legendTitle))
  }
  json$portions <- infoJs$portions
  json$lonMin <- infoJs$lonMin
  json$lonMax <- infoJs$lonMax
  json$lonNum <- infoJs$lonNum
  json$latMin <- infoJs$latMin
  json$latMax <- infoJs$latMax
  json$latNum <- infoJs$latNum
  json$minVal <- infoJs$minVal
  json$maxVal <- infoJs$maxVal
  json$varType <- infoJs$varType
  json$offsetType <- offsetType
  json$sizeType <- sizeType
  json$projection <- paste0("EPSG:", infoJs$epsg)

  jsonString <- toJSON(json, auto_unbox = TRUE, pretty = !minify)

  write(jsonString, file = file, append = FALSE)
  return(jsonString)
}
