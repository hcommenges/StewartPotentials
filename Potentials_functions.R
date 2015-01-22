


######################################
########### MAIN FUNCTIONS ###########
######################################


PointsOpport <- function(knownpts,
                         varname,
                         span,
                         beta,
                         mask = FALSE, 
                         typedist = "euclidean", 
                         typefct = "exponential", 
                         resolution = 200)
{
  if(is.projected(knownpts) == FALSE) {
    stop("Your input (knownpts) isn't a spatial object or isn't projected")
  }
  
  getproj <- proj4string(knownpts)
  unknownpts <- CreateGrid(knownpts = knownpts, mask = mask, resolution = resolution)
  matdist <- GetDistMat(knownpts = knownpts, unknownpts = unknownpts, typedist = typedist)
  interact <- ComputeInteractDensity(matdist = matdist, unknownpts = unknownpts, typefct = typefct, beta = beta, span = span)
  unknownpts <- interact$PTS
  matdens <- interact$MAT
  matopport <- ComputeOpportunity(knownpts = knownpts, matdens = matdens, varname = varname)
  return(list(UNKWPTS = unknownpts, OPPORT = matopport, PROJ = getproj))
}

OpportPotentials <- function(opportgrid, nbclass = 8, mask = FALSE)
{
  polyPot <- ComputePotentials(unknownpts = opportgrid$UNKWPTS, 
                               matopport = opportgrid$OPPORT, 
                               myproj = opportgrid$PROJ,
                               nbclass = nbclass, 
                               mask = mask)
  return(polyPot)
}


# Discretize and map ----




######################################
######## BACKGROUND FUNCTIONS ########
######################################

# Create grid of unknown points ----

CreateGrid <- function (knownpts, mask, resolution)
{
  coordPts <- data.frame(COORDX = coordinates(knownpts)[, 1], COORDY = coordinates(knownpts)[, 2], stringsAsFactors = FALSE)
  spatSteps <- ifelse(max(coordPts[ , "COORDX"]) - min(coordPts[ , "COORDX"]) >= max(coordPts[ , "COORDY"]) - min(coordPts[ , "COORDY"]),
                      round((max(coordPts[ , "COORDX"]) - min(coordPts[ , "COORDX"])) / resolution, digits = 4),
                      round((max(coordPts[ , "COORDY"]) - min(coordPts[ , "COORDY"])) / resolution, digits = 4))
  
  if(is.null(dim(mask))){
    boundingBox <- bbox(knownpts)
  } else {
    boundingBox <- bbox(mask)
  }
  
  boxCoordX <- round(seq(boundingBox[1, 1] - spatSteps, boundingBox[1, 2] + spatSteps, spatSteps), digits = 4)
  boxCoordY <- round(seq(boundingBox[2, 1] - spatSteps, boundingBox[2, 2] + spatSteps, spatSteps), digits = 4)
  spatGrid <- expand.grid(boxCoordX, boxCoordY)
  idSeq <- seq(1, nrow(spatGrid), 1)
  spatGrid <- data.frame(ID = idSeq, COORDX = spatGrid[, 1], COORDY = spatGrid[, 2])
  return(spatGrid)
}


# Compute matrix of distances ----

GetDistMat <- function(knownpts, unknownpts, typedist)
{
  knownPts <- data.frame(COORDX = coordinates(knownpts)[, 1], COORDY = coordinates(knownpts)[, 2], stringsAsFactors = FALSE)
  
  if(typedist == "euclidean") {
    matDist <- rdist(knownPts[ , c("COORDX", "COORDY")], unknownpts[ , c("COORDX", "COORDY")])
  } else if(typedist == "orthodromic") {
    # x = longitud, y = latitud (decimal, WGS84)
    matDist <- rdist.earth(knownPts[ , c("COORDX", "COORDY")], unknownpts[ , c("COORDX", "COORDY")], miles = FALSE, R = NULL)
  } else {
    stop("Please choose a valid distance argument (typedist)")
  }
  
  diag(matDist) <- 0 
  return(round(matDist, digits = 8))
}


# Compute density of interaction ----

ComputeInteractDensity <- function(matdist, unknownpts, typefct, beta, span)
{
  if(typefct == "pareto") {
    alpha  <- (2 ^ (1 / beta) - 1) / span
    matDens <- (1 + alpha * matdist) ^ (-beta)
  } else if(typefct == "exponential") {
    alpha  <- log(2) / span ^ beta
    matDens <- exp(- alpha * matdist ^ beta)
  } else {
    stop("Please choose a valid interaction function argument (typefct)")
  }
  
  diag(matDens) <- 1
  matDens <- round(matDens, digits = 8)
  unknownpts$MAXDENSITY <- apply(matDens, 2, max, na.rm = TRUE)
  
  return(list(PTS = unknownpts, MAT = matDens))
}


# Compute opportunities ----

ComputeOpportunity <- function(knownpts, matdens, varname = varname)
{
  matOpport <- knownpts@data[, varname] * matdens
  return(round(matOpport, digits = 8))
}


# Compute potentials ----

ComputePotentials <- function(unknownpts, matopport, nbclass, mask, myproj)
{
  # compute potentials and create raster
  unknownpts$POTENTIALS <- apply(matopport, 2, sum, na.rm = TRUE)
  equalInterval <- (max(unknownpts$POTENTIALS) - min(unknownpts$POTENTIALS)) / nbclass
  brksVal <- seq(min(unknownpts$POTENTIALS), max(unknownpts$POTENTIALS), equalInterval)
  
  unknownpts$POTDISCRET <- as.integer(cut(unknownpts$POTENTIALS, 
                                          breaks = brksVal,
                                          labels = seq(1, nbclass, 1), 
                                          right = FALSE, 
                                          include.lowest = TRUE))
  
  spatUnknownPts <- SpatialPointsDataFrame(coords = unknownpts[ , c(2, 3)], 
                                           data = unknownpts, 
                                           proj4string = CRS(myproj))

  rastGrid <- raster(xmn = min(unknownpts$COORDX), 
                     xmx = max(unknownpts$COORDX),
                     ymn = min(unknownpts$COORDY),
                     ymx = max(unknownpts$COORDY),
                     nrows = length(unique(unknownpts$COORDY)),
                     ncols = length(unique(unknownpts$COORDX)),
                     crs = CRS(myproj))
  
  rastFilled <- rasterize(spatUnknownPts, rastGrid, field = "POTDISCRET")
  
  if(!is.null(dim(mask))){
    rastFilled <- mask(rastFilled, mask = mask)
  }
  
  getPoly <- rasterToPolygons(rastFilled, n = 8, digits = 16, dissolve = TRUE)
  return(getPoly)
}




