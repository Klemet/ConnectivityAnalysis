# -*- coding: utf-8 -*-
"""
@author: Clement Hardy (clem.hardy@pm.me)

The goal of this script is to run a Circuitscape Omnidirectional analysis 
(with current going from North to South, South to North, East to West and
 West to East), with the resistance map having been generated with the script
CreatingResistanceMap_PitherEtAl2023.py (see file
resistanceCostMapPitherEtAl2023.tif).

To run circuitscape, you will need to install the juliacall Python module.
See https://juliapy.github.io/PythonCall.jl/stable/juliacall/ . I don't think
the module needs Julia to be installed on your computer, but you might want to
do it just in case (see https://julialang.org/downloads/).
"""

#%% LOADING PACKAGES

import subprocess
import os, glob
# import importlib
# import geopandas as gpd
import rasterio
# from rasterio.features import rasterize
# from rasterio.enums import Resampling
# from rasterio.mask import mask
import numpy as np
# import matplotlib.pyplot as plt
# from rasterio.enums import Compression

#%% LOADING FUNCTIONS

def getRasterData(path):
    """
    Get the data from a raster file with the rasterio package (avoids using GDAL),
    and put it in a numpy array for easy operations.
    """
    with rasterio.open(path) as src:
        # Read the raster data into a numpy array
        raster_data = src.read(1)  # Read the first band
    return(raster_data)


def writeNewRasterData(pathToSave, dataArray, dtypeRaster = rasterio.int32):
    """
    Writes a numpy array to a raster file.
    The extent of the raster file is set inside the function, and correspond to the area of the
    territory in Mauricie that we studied.
    
    /!\ YOU SHOULD CHANGE THE EXTENT HERE TO THE ONE OF YOUR STUDY AREA
    
    The resolution is set to 100x100m.
    """
    extent = (-551184.1237665037624538,
              -320384.1237665037624538,
              289450.5385142607265152,
              462250.5385142607265152)

    resolution = (100 , 100)

    cols = int((extent[1] - extent[0]) / resolution[0])
    rows = int((extent[3] - extent[2]) / resolution[1])
    transform = rasterio.Affine(resolution[0], 0, extent[0], 0, -resolution[1], extent[3])
    
    rasterToSave = rasterio.open(
        pathToSave,
        'w',
        driver='GTiff',
        height=rows,
        width=cols,
        count=1,
        dtype=dtypeRaster,
        crs='EPSG:32198',  # Coordinate Reference System (CRS)
        transform=transform
    )
    
    rasterToSave.write(dataArray, indexes=1)
    
    rasterToSave.close()

def convertTIFToASCRaster(pathToTifRaster):
    """
    Converts a .tif raster to a .asc raster.
    Circuitscape works better with .asc rasters, so it's needed for that.'
    """
    with rasterio.open(pathToTifRaster) as src:
        # Read the raster data
        data = src.read(1)  # Read the first band
        
        # Get the raster metadata
        meta = src.meta
        meta.update(driver='AAIGrid')  # Update the driver to AAIGrid for ASCII raster
        # Apparently, NoData value for Circuitscape needs to be -9999..??
        # https://www.researchgate.net/post/Why_do_I_keep_getting_an_error_when_trying_to_run_Circuitscape
        meta.update(nodata=-9999)
        

        # Create a new ASCII raster file
        with rasterio.open(pathToTifRaster[0: -3] + "asc", 'w', **meta) as dst:
            # Write the data to the new ASCII raster file
            dst.write(data, 1)

#%% READING RESISTANCE MAP
# I assume that it has been created with the previous script. I use the same folder structure.

resistanceMap = getRasterData(".\ResistanceCostMapCreation\CostMaps\MEDIUM (100x100m)\resistanceCostMapPitherEtAl2023.tif")

#%% MAKING THE SOURCE/GROUND RASTERS FOR CIRCUITSCAPE

# Circuitscape needs raster maps that indicate where the current is coming from
# (source) and where it is going (ground).
# Here, these maps will simply contain a band of 1 pixel equal to 1 on one of the side
# of the map (north, south, east or west). These maps can be used as ground or source
# in circuitscape, allowing us to do the omnidirectional analysis.

folderOfOutputsMaps_MEDIUM = r"./CircuitscapeAnalysis/Maps/MEDIUM (100x100m)"

# North Band
# We create a template array based on the resistance map
north_band = np.full(resistanceMap.shape, -9999)
north_band[0, :] = 1 # Set the first row to 1s
# We output
writeNewRasterData(folderOfOutputsMaps_MEDIUM + '/northBandWallToWallCircuitscape_MEDIUM.tif',
                   north_band, size = "MEDIUM")
convertTIFToASCRaster(folderOfOutputsMaps_MEDIUM + '/northBandWallToWallCircuitscape_MEDIUM.tif')
# East Band
east_band = np.full(resistanceMap.shape, -9999)
east_band[:, -1] = 1 # Set the last column to 1s
writeNewRasterData(folderOfOutputsMaps_MEDIUM + '/eastBandWallToWallCircuitscape_MEDIUM.tif',
                   east_band, size = "MEDIUM")
convertTIFToASCRaster(folderOfOutputsMaps_MEDIUM + '/eastBandWallToWallCircuitscape_MEDIUM.tif')

# West Band
west_band = np.full(resistanceMap.shape, -9999)
west_band[:, 0] = 1 # Set the first column to 1s
writeNewRasterData(folderOfOutputsMaps_MEDIUM + '/westBandWallToWallCircuitscape_MEDIUM.tif',
                   west_band, size = "MEDIUM")
convertTIFToASCRaster(folderOfOutputsMaps_MEDIUM + '/westBandWallToWallCircuitscape_MEDIUM.tif')

# South Band
south_band = np.full(resistanceMap.shape, -9999)
south_band[-1, :] = 1 # Set the last row to 1s
writeNewRasterData(folderOfOutputsMaps_MEDIUM + '/southBandWallToWallCircuitscape_MEDIUM.tif',
                   south_band, size = "MEDIUM")
convertTIFToASCRaster(folderOfOutputsMaps_MEDIUM + '/southBandWallToWallCircuitscape_MEDIUM.tif')


#%% USING CIRCUITSCAPE IN JULIA

# The latest version of Circuitscape (v5) is used directly in the Julia
# proramming langage.
# As such, you'll have to install Julia on your computer to use the
# following section. 
# - Install Julia using https://julialang.org/downloads/
#   WARNING : Go to "OFFICIAL BINARIES FOR MANUAL DOWNLOAD" below in the page;
#   first method didn't work for me (see https://julialang.org/downloads/#official_binaries_for_manual_download)
# - Install circuitscape in Julia (launch a julia prompt with "julia" in any terminal)
#   using Pkg
#   Pkg.add("Circuitscape")

# We then run Circuitscape by running a Julia script that contains the instructions.
# See JuliaScript_RunCircuitscapeComputations.jl.
# The script is very simple; but running circuitscape requires "ini files" which
# define what you want to do in your simulation.
# These files are in the folder "CircuitScape_IniFiles"; make sure that they point
# to the right location for the input and output files of the circuitscape runs.

# You can also run the commands in the Julia script by yourself in a command prompt after calling
# Julia. This makes it easier to debug if something goes wrong.

cmd = 'julia JuliaScript_RunCircuitscapeComputations.jl'
subprocess.run(cmd, shell=True)

#%% AVERAGING THE OUTPUT RASTERS

# Circuitscape has been run into four directions
# Now, we will create a map that will sum the currents.
folderOfOutputs_MEDIUM = r"./CircuitscapeAnalysis/Outputs/MEDIUM (100x100m)"
NorthToSouthRaster = getRasterData(folderOfOutputs_MEDIUM + "/NorthToSouth_MEDIUM_Outputs_curmap.asc")
SouthToNorthRaster = getRasterData(folderOfOutputs_MEDIUM + "/SouthToNorth_MEDIUM_Outputs_curmap.asc")
EastToWestRaster = getRasterData(folderOfOutputs_MEDIUM + "/EastToWest_MEDIUM_Outputs_curmap.asc")
WestToEastRaster = getRasterData(folderOfOutputs_MEDIUM + "/WestToEast_MEDIUM_Outputs_curmap.asc")

stacked_arrays = np.stack([NorthToSouthRaster, SouthToNorthRaster, EastToWestRaster, WestToEastRaster], axis=0)
avg_array = np.mean(stacked_arrays, axis=0)

writeNewRasterData(folderOfOutputs_MEDIUM + '/averageCurrent4WallToWall_MEDIUM.tif',
                   avg_array, size = "MEDIUM", dtypeRaster=rasterio.float32)

# You can now go and check the map in QGIS.
# One advice : tweak the maximum value of your pseudocolor coloring to better
# visualize the pinch points in your landscape !