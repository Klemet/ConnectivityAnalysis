# -*- coding: utf-8 -*-
"""
@author: Clement Hardy (clem.hardy@pm.me)

This script is used to define ecological corridors of mouvement in a landscape
based on a resistance map (created with previous scripts in this repository),
and the methodology of the model Linkage Mapper. I couldn't use Linkage Mapper
directly because of the cost of a ArcGIS licence.

WARNING : The script will require a temporary folder where it will generate
a large number of rasters; one for every individual protected area to connect
with corridors in your landscape. This folder can easily reach 10-20GB of disk space.
This is due to avoid keeping the raster in the RAM, which would be complicated.
There might be a better way to do this, but I don't have the time to improve the script.
Be sure to have such space ready ! 


The methodology of LinkageMapper used is described in detail https://github.com/linkagescape/linkage-mapper/blob/main/toolbox/doc/Linkage%20Pathways%20Linkage%20Mapper%20User%20Guide.docx
Summary :
    
1. Load pairs of patches on a MPG defined with grainscape or Graphlab
2. For each patch, create a map of cost-weighted distance to the patch
3. For each pair of patches, compute a corridor raster
4. Fuse all corridor rasters together by keeping the smallest value in each

The script needs to save all of the cost-weighted distance raster to avoid saving
them in the RAM (around 10GB for 500 patches and a landscape of the size of our study area).
"""

#%% LOADING PACKAGES

import os
import subprocess
import geopandas as gpd
import rasterio
import numpy as np
from tqdm import tqdm
import heapq

import whitebox
wbt = whitebox.WhiteboxTools()

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

def write_raster_from_array(output_path, data_array, template_raster_path):
    """
    Write a new raster file based on a numpy array and a template raster.

    Parameters:
    - output_path (str): Path where the new raster will be saved.
    - data_array (numpy.ndarray): 2D numpy array containing the values for the new raster.
    - template_raster_path (str): Path to the template raster file.

    Returns:
    None
    """
    # Open the template raster to get metadata
    with rasterio.open(template_raster_path) as template:
        # Check if the data_array shape matches the template raster shape
        if data_array.shape != template.shape:
            raise ValueError("The input array shape does not match the template raster shape.")

        # Copy the metadata from the template
        metadata = template.meta.copy()

        # Update the metadata with the new data type if necessary
        metadata.update(dtype=data_array.dtype)

        # Write the new raster
        with rasterio.open(output_path, 'w', **metadata) as dst:
            dst.write(data_array, 1)  # Write to the first (and only) band

    # print(f"New raster has been written to: {output_path}")
    

# The following functions are used to find the least-cost path between two cells
# in a numpy array.
# It's needed to compute the corridors properly, there is an issue in the Grainscape
# package that results in the wrong least-cost path being computed, making for
# negative values (see https://github.com/achubaty/grainscape/issues/72)
# Courtesy of perplexity :
# https://www.perplexity.ai/search/in-python-implement-the-follow-dhPPzjk.T0W1GDNa8dEiEw
def neighbors(rows, cols, r, c):
    for dr, dc in [(-1, 0), (1, 0), (0, -1), (0, 1)]:
        rr, cc = r + dr, c + dc
        if 0 <= rr < rows and 0 <= cc < cols:
            yield rr, cc

def dijkstra(resistanceRaster, patchesRaster, patchID):
    
    rows, cols = resistanceRaster.shape
    unique_patch_ids = np.unique(patchesRaster)

    # Initialize the cost array with infinity
    cost = np.full((rows, cols), np.inf)
    pq = []

    # Initialize the priority queue with the starting points (patch_id locations)
    for r in range(rows):
        for c in range(cols):
            if patchesRaster[r, c] == patchID:
                cost[r, c] = 0
                heapq.heappush(pq, (0, r, c))

    while pq:
        current_cost, r, c = heapq.heappop(pq)

        if current_cost > cost[r, c]:
            continue

        for rr, cc in neighbors(rows, cols, r, c):
            new_cost = current_cost + resistanceRaster[rr, cc]
            if new_cost < cost[rr, cc]:
                cost[rr, cc] = new_cost
                heapq.heappush(pq, (new_cost, rr, cc))

    # At this point, `cost` contains the cumulative cost distance for the current patch_id
    # print(f"Cumulative cost distance for patch ID {patch_id}:")
    # print(cost)
    return(cost)

def find_least_cost_path(cost_raster, patch_raster, start_id, end_id):
    # Find all possible start coordinates
    start_coords = np.argwhere(patch_raster == start_id)
    
    # Edit cost raster so that end and start = 0
    cost_raster_function = copy.deepcopy(cost_raster)
    # cost_raster_function = np.where(patch_raster == start_id, 0, cost_raster_function)
    # cost_raster_function = np.where(patch_raster == end_id, 0, cost_raster_function)
    
    # Initialize data structures
    rows, cols = cost_raster.shape
    distances = np.full((rows, cols), np.inf)
    predecessors = np.full((rows, cols, 2), -1, dtype=int)
    
    # Priority queue for Dijkstra's algorithm
    pq = []
    
    # Initialize all start points
    for start in start_coords:
        distances[tuple(start)] = 0
        heapq.heappush(pq, (0, start[0], start[1]))
    
    # Possible movements
    directions = [(-1, 0), (1, 0), (0, -1), (0, 1), (1, 1), (-1, 1), (1, -1), (-1, -1)]
    
    while pq:
        current_cost, r, c = heapq.heappop(pq)
        
        # Check if the current cell is an end cell
        if patch_raster[r, c] == end_id:
            end = [r, c]
            break
        
        for dr, dc in directions:
            nr, nc = r + dr, c + dc
            if 0 <= nr < rows and 0 <= nc < cols:
                movement_cost = ((cost_raster_function[r, c] + cost_raster_function[nr, nc])/2) * (math.sqrt(dr * dr + dc * dc))
                new_cost = current_cost + movement_cost
                
                if new_cost < distances[nr, nc]:
                    distances[nr, nc] = new_cost
                    predecessors[nr, nc] = [r, c]
                    heapq.heappush(pq, (new_cost, nr, nc))
    
    # Reconstruct path
    path = []
    current = end
    while not np.array_equal(current, predecessors[tuple(current)]):  # Stop when we reach a start point
        path.append(current)
        current = predecessors[tuple(current)]
    path.append(current)  # Add the start point
    path.reverse()
    
    return distances[tuple(end)]


#%% LOADING FILES

folderWithData = r"./CorridorAnalysis/"

# Need two raster maps : 
    
# Resistance map
resistanceRaster = getRasterData(folderWithData + "../ResistanceCostMapCreation/CostMaps/MEDIUM (100x100m)/resistanceCostMapPitherEtAl2023.tif")
# plt.imshow(resistanceRaster)

# Map with patches

# Retrieve the MPG (Minimal Planar Graph) data made with the grainscape R package
# We need the map of the patches (that contain the patch ID)
# And the links shapefile that we transform into a dictionnary with the list of
# links as key, and the least-cost path value as values.
# We then use these links to explore the possible corridors between patches.
# If we were doing every possible combination/link of patches, this would take
# an incredible amount of time. Plus, Linkage Mapper seems to use a MPG too.

# WARNING : BE SURE TO RUN THE R SCRIPT TO CREATE THE MPG DATA FIRST ! (createMPG.R)
    
patchesRaster = getRasterData(folderWithData + "MinimalPlanarGraph/patchId.tif")
# plt.imshow(patchesRaster)

shapefileMPGLinks = gpd.read_file(folderWithData + "MinimalPlanarGraph/linksCentroid.shp")

# We create a dictionnary of links; it associate the weight value (probability of mouvement)
# of the link to each link "name"
linksDictionnary = dict()

for i in range(0,len(shapefileMPGLinks)):
    linkName = str(shapefileMPGLinks["e1"][i]) + "-" + str(shapefileMPGLinks["e2"][i])
    linksDictionnary[linkName] = shapefileMPGLinks["lcpPerWt"][i]
    
#%% CREATING COST-WEIGHTED RASTERS
# WARNING : This step generates around 10GB of temporary rasters (necessary
# for the next step) !
# It saves them locally to avoid using RAM.
    
# folderToSaveCostWeightedDistanceRasters = folderWithData + "costWeightedDistanceRasters\\"
folderToSaveCostWeightedDistanceRasters = r"./Temp_Rasters_Corridors/"

allPatchID = np.unique(patchesRaster).astype(int)
allPatchID = allPatchID[allPatchID > 0]

# Using the whitebox package. Much much quicker. Needs to create some
# temporary rasters that we remove, but is much better overall.
overwriteRasters = False # Use to re-do the rasters in needed
tempPatchRasterPath = folderToSaveCostWeightedDistanceRasters + "TEMPRASTERPATCH.tif"

dictOfCWD_Rasters_Paths = dict() #Saves the path of each raster created
for patchID in tqdm(allPatchID):
    rasterOutputPath = folderToSaveCostWeightedDistanceRasters + "CWD_Patch_" + str(patchID) + ".tif"
    # No need to create the raster if it already exists
    if not os.path.exists(rasterOutputPath) or not overwriteRasters:
        
        # We create a raster indicating the position of the patches
        maskedArray = np.zeros_like(patchesRaster)
        maskedArray[patchesRaster == patchID] = 1
        
        # We create the temporary raster with the position of the patch, nécéssary for Whitebox
        write_raster_from_array(tempPatchRasterPath,
                                maskedArray, 
                                folderWithData + "../ResistanceCostMapCreation/CostMaps/MEDIUM (100x100m)/resistanceCostMapPitherEtAl2023.tif")
        
        # Here is the function from the Whitebox package.
        # However, the function wants to output a "backlink" raster, which I don't
        # need. I simply delete it right after.
        backlinkRasterOutputPath = folderToSaveCostWeightedDistanceRasters + "CWD_Patch_" + str(patchID) + "_BACK.tif"
        wbt.cost_distance(
        source=tempPatchRasterPath,
        cost=folderWithData + "../ResistanceCostMapCreation/CostMaps/MEDIUM (100x100m)/resistanceCostMapPitherEtAl2023.tif",
        out_accum=rasterOutputPath,
        out_backlink=backlinkRasterOutputPath)
        # We delete the backlink raster
        os.remove(backlinkRasterOutputPath)
        # We save the raster path for latter
        dictOfCWD_Rasters_Paths[patchID] = rasterOutputPath
    
os.remove(tempPatchRasterPath)

#%% CREATE CORRIDOR RASTERS AND COMBINE THEM

# Preparing the final, composite raster of corridors.
finalCorridorsRaster = np.full_like(resistanceRaster, np.inf, dtype=float)

# Getting raster cell size to edit the values of the CWD raster later
# This is because the source code of whitebox computes the CWD by multiplicating every "movement"
# by the distance between cells centroids (in m2). We don't want that here.
with rasterio.open(folderWithData + "../ResistanceCostMapCreation/CostMaps/MEDIUM (100x100m)/resistanceCostMapPitherEtAl2023.tif") as src:
    cell_size_x, cell_size_y = src.res

# For each link of the MPG, we get the two weighted distance rasters made by
# Whitebox, and combine them according to the equation of linkage mapper to get
# a corridor raster between the two patches.
# Then, we immediatly combine this raster into the composite final raster
# by putting the minimal values between the two rasters into the composite final raster.
for link in tqdm(linksDictionnary.keys()):
    # Extracting patch names
    firstPatch = int(link.split("-")[0])
    secondPatch = int(link.split("-")[1])
    
    # Recovering patch CWD raster data for the two patches
    firstPatchCWD_Raster = getRasterData(dictOfCWD_Rasters_Paths[firstPatch])
    secondPatchCWD_Raster = getRasterData(dictOfCWD_Rasters_Paths[secondPatch])

    # Creating corridor raster (removing the addition of the distance cost in meters added by whitebox)
    corridorRaster = (firstPatchCWD_Raster + secondPatchCWD_Raster)/((cell_size_x + cell_size_y)/2)
    
    # WARNING
    # The least cost paths from grainscape don't seem good. I'll recompute them here
    # (see "Functions" section at the top of the script).
    # As such, we keep the MPG links as predicted by MPG; but we in effect re-trace
    # the least cost path for each link. In the absolute, this should not mean any change
    # to the structure of the MPG, but only to the weight of the links.
    linksDictionnary[link] = find_least_cost_path(resistanceRaster,
                                                    patchesRaster,
                                                    firstPatch,
                                                    secondPatch)
    # We compute the corridor raster between the two patches according to the
    # linkage mapper approach; it's the sum of the two cost-weighted distance raster
    # from the patch, substracted in every pixel by the cost of the least-cost path
    # between the two patches.
    corridorRaster = corridorRaster - linksDictionnary[link]
    
    # There should not be any negative values in the resulting corridor raster.
    # If that's the case, there is a problem in the cost-weighted distance raster
    # or in the least-cost path computed. Check that both have the right units of costs.
    if np.any(corridorRaster < 0):
        raise Exception("Negative values detected for link " + str(link))
    
    # Influencing the final corridor merged raster
    finalCorridorsRaster = np.minimum(finalCorridorsRaster, corridorRaster)

# plt.imshow(finalCorridorsRaster, vmax = 1000)
# Outputing the final composite raster of corridors.
write_raster_from_array(folderWithData + "finalCorridorsRaster.tif",
                        finalCorridorsRaster, 
                        folderWithData + "../ResistanceCostMapCreation/CostMaps/MEDIUM (100x100m)/resistanceCostMapPitherEtAl2023.tif")