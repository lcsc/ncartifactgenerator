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
#' @param portion Portion of the netCDF file to process.
#'
#' @return List with the generated web configuration.
#'
#' @export
config_web <- function(file, folder, epsg, formatdates, varmin, varmax, varName, infoJs = NA, lon_name, lat_name, time_name = "time", portion) {
  if (missing(infoJs) || sum(!is.na(infoJs)) == 0) {
    infoJs <- list(
      times = list(),
      latMin = list(), latMax = list(), latNum = list(),
      lonMin = list(), lonMax = list(), lonNum = list(),
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

  # read spatial dims
  if (missing(lon_name) || missing(lat_name)) {
    dimNames <- returnXYNames(nc)
    lon_name <- dimNames$X
    lat_name <- dimNames$Y
  }

  coords <- read_coords(nc, epsg)
  infoJs$lonMin[[paste0(varName, portion)]] <- min(coords[, 1])
  infoJs$lonMax[[paste0(varName, portion)]] <- max(coords[, 1])
  infoJs$lonNum[[paste0(varName, portion)]] <- length(unique(coords[, 1]))
  infoJs$latMin[[paste0(varName, portion)]] <- min(coords[, 2])
  infoJs$latMax[[paste0(varName, portion)]] <- max(coords[, 2])
  infoJs$latNum[[paste0(varName, portion)]] <- length(unique(coords[, 2]))

  infoJs$portions[[varName]] <- I(c(infoJs$portions[[varName]], portion))

  infoJs$varType <- get_struct_typecode(nc$var[[var_name]]$prec)

  # Available times
  if (!missing(formatdates)) {
    times.write <- read_times(nc, formatdates)
  } else {
    times.write <- read_times(nc)
  }

  if (missing(varmin) | missing(varmax)) {
    varMinMax <- readMinMax(nc)
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


#' Write JavaScript file with configuration data
#'
#' This function writes a JavaScript file with configuration data for VisorServiciosClimaticos.
#'
#' @param folder The folder where the JavaScript file will be saved.
#' @param infoJs A list containing the configuration data, including times, variable minimum and maximum values, latitude and longitude ranges, and EPSG code.
#' @param varTitle The title of the variable.
#' @param legendTitle The title of the legend (optional, default is "Legend").
#' @param offsetType The type of offset (optional, default is "Q").
#' @param sizeType The type of size (optional, default is "I").
#'
#' @return The generated JavaScript code as a character string.
#'
#' @export
writeJs <- function(folder, infoJs, varTitle, legendTitle = "Legend", offsetType = "Q", sizeType = "I", minify = TRUE) {
  file <- file.path(folder, "times.js")

  if (missing(varTitle)) {
    if (length(infoJs$varmin) > 1) {
      varTitle <- names(infoJs$varmin)
      names(varTitle) <- names(infoJs$varmin)
    } else {
      varTitle <- names(infoJs$varmin)
    }
  }

  text.js <- ""

  text.js <- paste(text.js, paste0("var center = {'lat': ", infoJs$latM, ", 'lng': ", infoJs$lonM, "};\n"))

  text.js <- paste(text.js, arrayRtojs(name = "times", value = infoJs$times, type = "date"))
  text.js <- paste(text.js, arrayRtojs(name = "varMin", value = infoJs$varmin, type = "numeric"))
  text.js <- paste(text.js, arrayRtojs(name = "varMax", value = infoJs$varmax, type = "numeric"))
  text.js <- paste(text.js, arrayRtojs(name = "varTitle", value = varTitle))
  if (length(legendTitle) > 1) {
    text.js <- paste(text.js, arrayRtojs(name = "legendTitle", value = legendTitle))
  } else {
    text.js <- paste(text.js, paste0("var legendTitle = {NaN:['", legendTitle, "']};\n"))
  }
  text.js <- paste(text.js, arrayRtojs(name = "portions", value = infoJs$portions, type = "character"))
  text.js <- paste(text.js, arrayRtojs(name = "lonMin", value = infoJs$lonMin, type = "numeric", value_array = FALSE))
  text.js <- paste(text.js, arrayRtojs(name = "lonMax", value = infoJs$lonMax, type = "numeric", value_array = FALSE))
  text.js <- paste(text.js, arrayRtojs(name = "lonNum", value = infoJs$lonNum, type = "numeric", value_array = FALSE))
  text.js <- paste(text.js, arrayRtojs(name = "latMin", value = infoJs$latMin, type = "numeric", value_array = FALSE))
  text.js <- paste(text.js, arrayRtojs(name = "latMax", value = infoJs$latMax, type = "numeric", value_array = FALSE))
  text.js <- paste(text.js, arrayRtojs(name = "latNum", value = infoJs$latNum, type = "numeric", value_array = FALSE))
  text.js <- paste(text.js, paste0("var varType = '", infoJs$varType, "';\n"))
  text.js <- paste(text.js, paste0("var offsetType = '", offsetType, "';\n"))
  text.js <- paste(text.js, paste0("var sizeType = '", sizeType, "';\n"))
  text.js <- paste(text.js, paste0("var projection = 'EPSG:", infoJs$epsg, "';\n"))

  if (minify) {
    text.js <- uglify_optimize(text.js)
  }

  write(text.js, file = file, append = FALSE)
  return(text.js)
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

  json$center <- list(lat = latM, lng = lonM)
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
  json$varType <- infoJs$varType
  json$offsetType <- offsetType
  json$sizeType <- sizeType
  json$projection <- paste0("EPSG:", infoJs$epsg)

  jsonString <- toJSON(json, auto_unbox = TRUE, pretty = !minify)

  write(jsonString, file = file, append = FALSE)
  return(jsonString)
}
