#####################################################
#
# @author: Clement Hardy (clem.hardy@pm.me)
#
# This script generates a Minimal Planar Graph for the
#for the rest of the corridor analysis.
#####################################################



#### LOADING PACKAGES ####

# install.packages("grainscape")
library(grainscape)
library(igraph)
library(raster)
library(sf)
library(Dict)
library(rlang)


#### LOADING DATA ####

# Resistance map (see CreatingResistanceMap_PitherEtAl2023.py for details)

directoryWithFiles = "./CorridorAnalysis"

resistanceMap <- raster(paste(directoryWithFiles, "../ResistanceCostMapCreation/CostMaps/MEDIUM (100x100m)/resistanceCostMapPitherEtAl2023.tif", sep = ""))

# Patches raster map 

patchs_SacredMountainsAndProtectedAreas <- raster(paste(directoryWithFiles, "raster_ToutesAiresProtégées_MontSacr_ExtendedGouin.tif", sep = ""))

# Combining the two for the grainscape functions (else, always bug because of extent) :
# Patches will have the value 99999. 

resistanceMapAndPatches = resistanceMap
resistanceMapAndPatches[patchs_SacredMountainsAndProtectedAreas == 1] <- -999

# To check if all is OK

# plot(resistanceMapAndPatches)

#### CREATING AND DISPLAYING THE MPG ####

patchyMPG <- MPG(resistanceMapAndPatches, patch = resistanceMapAndPatches == -999)

# plot(patchyMPG, quick = "mpgPlot", theme = FALSE)

# Exports the data as shapefiles for visualisation
pathOfExport = "./CorridorAnalysis/MinimalPlanarGraph"
directoryOfExport = "CreateMinimalPlanarGraph"
fullPathOfExport = paste(pathOfExport, directoryOfExport, sep = "")
export(patchyMPG, dirname = directoryOfExport,
       path = pathOfExport,
       vorBound = FALSE, overwrite = TRUE)