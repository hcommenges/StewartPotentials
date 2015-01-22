
#' Stewart's Potentials
#'
#' Compute and plot Stewart's Potentials
#'
#' This function compute the potentials of spatial interaction as defined by J.Q. Stewart (1950) and plot the results.
#' It provides two kinds of outputs: a data.frame or a raster map (ggplot) representing potentials of interaction. 
#'
#' @param knownpts, object of class sp (SpatialPointsDataFrame or SpatialPolygonsDataFrame).
#' @param mask, object of class sp (SpatialPolygonsDataFrame) to clip the final map of potentials.
#' @param opportgrid, output of the PointOpport() function. This object is a list containing a grid of unknown points, 
#' a matrix of opportunities and a projection.
#' @param varname CHARACTER, name of the variable (in the attribute table) from which potentials are computed.
#' @param typedist CHARACTER, type of distance. Options are "euclidean" (default) or "orthodromic". 
#' If the distance is euclidean, X and Y coordinates are expected. If orthodromic longitud and latitud are expected 
#' in decimal degrees, in WGS84 reference system.)
#' @param typefct CHARACTER, spatial interaction function. Options are "pareto" (default) or "exponential".
#' If "pareto" the interaction is defined as: (1 + alpha * mDistance) ^ (-beta).
#' If "exponential" the interaction is defined as: exp(- alpha * mDistance ^ beta).
#' The alpha parameter is computed from parameters given by the user (beta and reach).
#' @param beta NUMERIC, impedance factor for the spatial interaction function.
#' @param span NUMERIC, distance where the density of probability of the spatial interaction function equals 0.5.
#' @param nbclass INTEGER, number of classes for the cartographic representation.
#' @param resolution INTEGER, resolution of the grid for raster representation.
#' @references Stewart J. Q., Demographic gravitation: evidence and applications, Sociometry, 11(1-2), 31-58.


# load packages ----

library(sp)
library(rgdal)
library(rgeos)
library(raster)
library(fields)
library(RColorBrewer)


# load functions ----

source("Potentials_functions.R")


# load data ----

spatPts <- readOGR(dsn = "DATA/parispc_pts.shp", layer = "parispc_pts", stringsAsFactors = FALSE)
spatUnits <- readOGR(dsn = "DATA/parispc_com.shp", layer = "parispc_com", stringsAsFactors = FALSE)
spatMask <- readOGR(dsn = "DATA/parispc_mask.shp", layer = "parispc_mask", stringsAsFactors = FALSE)

plot(spatMask)
plot(spatUnits)
plot(spatPts, add = T)


# compute opportunities ----

opportunities <- PointsOpport(knownpts = spatPts, 
                              varname = "POPULATION",
                              span = 1500, 
                              beta = 2, 
                              mask = spatMask, 
                              resolution = 400)


# compute potentials ----

potentials <- OpportPotentials(opportgrid = opportunities, 
                               nbclass = 8, 
                               mask = spatMask)


# plot potentials ----

plot(potentials, col = plyr::mapvalues(potentials$layer, 
                                       from = potentials$layer, 
                                       to = brewer.pal(n = 8, name = "Reds")))


# generalize polygons geometry ----

library(spgrass6)
initGRASS(gisBase = "/usr/lib/grass64", home = getwd(), override = TRUE)
writeVECT6(SDF = potentials, vname = "inpol", v.in.ogr_flags = c("o", "overwrite"))
execGRASS("v.generalize", flags = c("overwrite"), input = "inpol", output = "outpol", threshold = 1, method = "snakes", alpha = 1, beta = 1)
genPol <- readVECT6(vname = "outpol")

plot(genPol, col = plyr::mapvalues(genPol$cat, 
                                   from = genPol$cat, 
                                   to = brewer.pal(n = 8, name = "Reds")))




