# -*- coding: utf-8 -*-
"""
@author: Clement Hardy (clem.hardy@pm.me)

The goal of this script is to create a resistance map to movement for 
non-flying terrestrial species (e.g. moose, wolves, caribou, etc.) with the same
methodology as the one used by Pither et al. 2023
(see https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0281980).

Their methodology is simple and clear, and yet seems to be validated on empirical
GPS mouvement data for several species.

This resistance map can then be used for the following analysis (see other folders of the
repository) :

- Circuitscape runs (to generate current through the resistance map)
- Graph/Network analysis (to define the cost of mouvement from node to node through least cost paths
on the resistance map)
- Corridor analysis (the corridor value of the corridor raster map correspond to how passing by a given
pixel is close to the least-cost path between the two closest patches)

The study of Pither et al. 2023 use 23 layers of spatial data 
that represent human or natural features that are obstacles
for the fauna in Canada : cities, roads, agricultural areas; but also slopes, glaciers,
etc. I didn''t look at all of the data layers they used, as I didn't need them in
the study area we considered (in Mauricie). More information on the layers below.

Here, I propose to produce maps of 100x100m, which are in a higher resolution than the one used
by Pither et al. 2023. I tried 20mx20m, but this resulted in very large raster maps that didn't produced
meaningfull differences in our results. The most important is making sure that linear obstacles (such as rivers)
are properly rasterized, meaning with no "holes" in them dues to how lines and polygon are often rasterized (a 
certain quantity of line or polygon has to be in the cell, or it is not rasterized).
"""

#%% LOADING PACKAGES

import os
import importlib
import geopandas as gpd
import rasterio
from rasterio.features import rasterize
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

def RasterizePolygons(polygonsAsGeoPandas, pathToTemplateRaster, valueToRasterizeInPolygons, allTouched = False):
    """
    Uses rasterio  to rasterize polygons loaded with Geopandas
    WARNING : valueToRasterizeInPolygons is the "Burn" value;
                it should be an integer that fits into a int32 raster.
    WARNING : The template raster and the polygons must have the same CRS !
    
    The "allTouched" parameter deals with using the "All Touched" function (i.e.
    even if a small piece of a line or a polygon touches a raster cell, it's burned
    in the raster cell). See https://gis.stackexchange.com/questions/354494/rasterizing-all-intersected-cells-when-rasterizing-line-layer-in-qgis/383867#383867
    for more infos.
    """
    with rasterio.open(pathToTemplateRaster) as src:
        
        out_meta = src.meta.copy()
    
        # Rasterize the polygons with desired values
        shapes = ((geom, valueToRasterizeInPolygons) for geom in polygonsAsGeoPandas.geometry)  # Polygon geometries with value 1000
        if not allTouched:
            rasterizedArray = rasterize(shapes, out_shape=(out_meta['height'], out_meta['width']), transform=out_meta['transform'], fill=1, dtype="int32")
        else:
            rasterizedArray = rasterize(shapes, out_shape=(out_meta['height'], out_meta['width']), transform=out_meta['transform'], all_touched=True, fill=1, dtype="int32")
        
    return(rasterizedArray)

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

#%% EXPLANATION BEFOREHAND

"""
List of layers used in Pither et al. 2023 and the associated resistance costs 
(according to the cost ponderation model C9 they used in the end).

Come from their appendix S2.

I removed data that are for the US only, as I'm only interested to our Quebec area.
CHF = Canadian Human Footprint, a 300x300m raster map covering all of Canada.
The data source for lakes, rails, roads and others are different datasets they
used.
    
                 LAYER                      COUNTRY           COST       
CHF - Built environments                Canada                       1000
CHF - Croplands                         Canada                        100
CHF - Dams                              Canada                       1000
CHF - Forestry (cut between 1985 & 2015)Canada                         10
CHF - Mining                            Canada                       1000
CHF - Nighttime lights                  Canada                       1000
CHF - Oil and gas                       Canada                       1000
CHF - Pasturelands                      Canada                         10
Lakes >= 10ha                           Canada                       1000
Rails                                   Canada                       1000
Roads - minor                           Canada                         10
Roads - two-lane highway                Canada                        100
Roads - multi-lane highways             Canada                       1000
Elevation > 2300m                       Canada & U.S.                1000
Glaciers                                Canada & U.S.                1000
Ocean                                   Canada & U.S.                1000
Rivers > 28m3/sec                       Canada & U.S.                1000
Sea Ice                                 Canada & U.S.                  10
Slopes > 30 degrees                     Canada & U.S.                1000
"""

"""
I propose to use finer data than the ones used by Pither et al., using the
ecoforest maps/polygons of Quebec and the AQ Réseau database.

With the "CO_TER" field of the polygones ecoforestiers, we can identify many of
the environments corresponding to the ones in the layers of Pither et al.

See https://www.donneesquebec.ca/recherche/dataset/carte-ecoforestiere-avec-perturbations
and download "Fiche descriptive des attributs et de leurs domaines de valeurs .xlsx"
for the details about each CO_TER code.

Here are all of the variations of the CO_TER field for the polygons in our study area
area, with their signification (in french) and the category of Pither et al. that
I am assigning them to :

CO_TER Code     Meaning                                                                          Pither et al. 2023 resistance category                    
A                Agricole                                                                        Cropland                                                  
AEP              Aire Empilement et ébranchage                                                   Forestry                                                  
AL               Aulnaie                                                                         None, low cost of traversal                               
ANT              Milieu fortement perturbé par l'activité humaine (non-boisé)                    Built environments                                        
DH               Dénudé et semi-dénudé humide                                                    None, low cost of traversal                               
DS               Dénudé et semi-dénudé sec                                                       None, low cost of traversal                               
EAU              Étendue d'eau, cours d'eau, réservoirs d'origine anthropique et battures        Lakes >= 10ha only if surface is >= 10ha. Rivers will come from another dataset and are ignored (see below).                                         
GR               Gravière                                                                        Mining                                                    
ILE              Île  superficie < 1 ha                                                          None, low cost of traversal                               
INO              Site inondé, site exondé non régénéré                                           None, low cost of traversal                               
LTE              Ligne de transport d'énergie                                                    None, low cost of traversal                               
NF               Milieu faiblement perturbé par l'activité humaine (boisé)                       Forestry                                                  
RO               Route et autoroute (emprise)                                                    Not super reliable, will use AQRéseau+ instead      

In addition, I will use several dataset to match the other categories :
    
# Forest operations
See https://www.donneesquebec.ca/recherche/dataset/recolte-et-reboisement
Will select any forest operation made since 2000. A bit arbitrary, but seems
appropriate. Those will additionaly go into the "forestry" category.

# AQRéseaux+

See https://www.donneesquebec.ca/recherche/dataset/adresses-quebec, "AQréseau+" in shp
Big dataset of roads and rails in Quebec.
Because it is vectorial, I will vectorise everything with the resolution I want
to use : 20m x 20m.
I will use the "all touch" method, where if the line of a road/rail touches a pixel
even slightly, then the pixel is considered "occupied" by the line. This avoids
"holes" in the rails/roads due to the rasterisation process.

- Reseau_ferroviaire.shp contains the rails (correspond directly to the category of Pither et al.)
- Reseau_routier.shp contains the roads.
    - To identify the minor, two-lane highways or multi-lane highways, I will use
      the attributes of the features
    - NbrVoies contains the number of lanes. if >= 2, two-lane highway; if => 3, multi-lane highway.
    - If number of lanes is INC or 1, I'll default to a minor road.

# Others

- No need for elevation because nowere superior to 2300m in our study area (might change for yours)
- No need for glaciers as they are none our region
- No need for the ocean either
- For the rivers, I'll use the same database as Pither et al. : HydroRIVERS
  (see https://www.hydrosheds.org/products/hydrorivers, North and Central America Shapefile)
  - Can use the attribute field "DIS_AV_CMS" that estimate the debit of water
    (with the formula DIS_AV_CMS" >= 28 to match the 28m3/s of Pither et al.'s category).
- See ice is not needed
- Slope will derived from digital elevation model data from https://www.donneesquebec.ca/recherche/dataset/modeles-numeriques-d-altitude-a-l-echelle-de-1-20-000
    - I used the "Slope" raster terrain analysis tool from QGIS to derive the slope from it

SUMMARY OF ALL OF THIS : See DataSourcesAndCostsForResistanceMap.ods
"""

#%% PREPARING OBJECTS

# The working directory for the creation of this resistance map.
# I will then assume that it contains a "SourceFiles" folder where the sources
# files of the spatial data will be located.
resistanceMapFolder = r".\ResistanceCostMapCreation"

# The folder where every output raster will be
folderOfOutputsMaps_MEDIUM = r".\ResistanceCostMapCreation\CostMaps\MEDIUM (100x100m)"

# We create a template raster that will serve as a model for all of the other rasters to export
# The extent of the rasters to be created :
# Corresponds to the extent (in EPSG:32198 - NAD83 / Quebec Lambert) of the
# rasters I've been using for our study area
# (xmin, xmax, ymin, ymax)
# /!\ YOU SHOULD CHANGE THE EXTENT FOR YOUR STUDY AREA HERE
extent = (-551184.1237665037624538,
          -320384.1237665037624538,
          289450.5385142607265152,
          462250.5385142607265152)
resolution = (100, 100)
cols = int((extent[1] - extent[0]) / resolution[0])
rows = int((extent[3] - extent[2]) / resolution[1])
transform = rasterio.Affine(resolution[0], 0, extent[0], 0, -resolution[1], extent[3])

template_raster_100x100 = rasterio.open(
    folderOfOutputsMaps_MEDIUM + '/template_raster_100x100.tif',
    'w',
    driver='GTiff',
    height=rows,
    width=cols,
    count=1,
    dtype=rasterio.uint8,
    crs='EPSG:32198',  # Coordinate Reference System (CRS)
    transform=transform
)

onesDataArray = np.ones((rows, cols), dtype=np.uint8)
template_raster_100x100.write(onesDataArray, indexes=1)

template_raster_100x100.close()

# Check if pyogrio is installed for faster loading of large datasets
try:
    importlib.import_module('pyogrio')
    pyogrio_installed = True
except ImportError:
    pyogrio_installed = False

#%% ECOFORESTED POLYGONS

# /!\ READ THIS : Starting from this section, I'll load the necessary spatial data
# and then use fonctions to compute the cost raster we need and export it.
# I'll indicate everytime how the spatial data to load has been generated so that
# you can do it yourself.

# The data here come from the Carte Ecoforestière Mise à Jour (MAJ) for the area
# Download all data for the region, and then clip the polygons
# to the extent of the region in QGIS.
# From https://www.donneesquebec.ca/recherche/dataset/carte-ecoforestiere-avec-perturbations
if pyogrio_installed:
    ecoforestedPolygons = gpd.read_file(r'.\ResistanceCostMapCreation\SourceFiles\Polygon_ecofor.shp', engine="pyogrio", use_arrow=True)
else:
    ecoforestedPolygons = gpd.read_file(r'.\ResistanceCostMapCreation\SourceFiles\Polygon_ecofor.shp', use_arrow=True)

# RASTER OF BUILT ENVIRONMENTS

# We select the polygons corresponding to the criteria
builtHumanPolygons = ecoforestedPolygons[ecoforestedPolygons['CO_TER'] == "ANT"]

# We rasterize them with the template raster,
# inserting the cost associated to built environments (1000) in Pither et al.
rasterizedBuiltHumanAreasCost = RasterizePolygons(builtHumanPolygons,
                                                  folderOfOutputsMaps_MEDIUM + '/template_raster.tif',
                                                  1000,
                                                  allTouched=True)

# We save the raster
writeNewRasterData(folderOfOutputsMaps_MEDIUM + '/builtAreasCost.tif',
                   rasterizedBuiltHumanAreasCost)

# RASTER OF CROPLANDS

# Same as precedent.
croplandPolygons = ecoforestedPolygons[ecoforestedPolygons['CO_TER'] == "A"]
# Here, cost is 100.
rasterizedCroplandCost = RasterizePolygons(croplandPolygons,
                                            folderOfOutputsMaps_MEDIUM + '/template_raster.tif',
                                            100,
                                            allTouched=True)
writeNewRasterData(folderOfOutputsMaps_MEDIUM + '/croplandsCosts.tif',
                   rasterizedCroplandCost)

# FORESTRY

# First, we take the ecoforest polygons with the right codes
forestryPolygonsEcoforestMap = ecoforestedPolygons[(ecoforestedPolygons['CO_TER'] == "NF") | (ecoforestedPolygons['CO_TER'] == "AEP")]
# We also take all forest operations since 2000
# From https://www.donneesquebec.ca/recherche/dataset/recolte-et-reboisement, "data since 1976"
allForestOperationsPolygons = gpd.read_file(r'.\ResistanceCostMapCreation\SourceFiles\Interventions_Forestières_Zone_Etude.shp', use_arrow=True)
allForestOperationsPolygons['an_origine'] = allForestOperationsPolygons['an_origine'].fillna(0).astype(int)
allForestOperationsPolygons['an_perturb'] = allForestOperationsPolygons['an_perturb'].fillna(0).astype(int)
forestryPolygonsForestOperations = allForestOperationsPolygons[(allForestOperationsPolygons['an_origine'] >= 2000) | (allForestOperationsPolygons['an_perturb'] >= 2000)]
# Here, cost is 10.
rasterizedForestryEcoforestMap = RasterizePolygons(forestryPolygonsEcoforestMap,
                                            folderOfOutputsMaps_MEDIUM + '/template_raster.tif',
                                            10,
                                            allTouched=True)
rasterizedForestryForestOperations = RasterizePolygons(forestryPolygonsForestOperations,
                                            folderOfOutputsMaps_MEDIUM + '/template_raster.tif',
                                            10,
                                            allTouched=True)
# We combine the arrays
forestryCosts = np.maximum(rasterizedForestryEcoforestMap, rasterizedForestryForestOperations)
writeNewRasterData(folderOfOutputsMaps_MEDIUM + '/forestryCosts.tif',
                   forestryCosts)

# MINING

# Same as precedent.
miningPolygons = ecoforestedPolygons[ecoforestedPolygons['CO_TER'] == "GR"]
# Here, cost is 1000.
rasterizedMiningCost = RasterizePolygons(miningPolygons,
                                            folderOfOutputsMaps_MEDIUM + '/template_raster.tif',
                                            1000,
                                            allTouched=True)
writeNewRasterData(folderOfOutputsMaps_MEDIUM + '/miningCosts.tif',
                   rasterizedMiningCost)

# LAKES > 10ha
# Instead of using ecoforest polygons, I use a custom shapefile where rivers
# and narrow parts of lakes have been removed based on a negative 50m buffer.
# The logic is that the ecoforest polygons of the 5th quebec forest inventory
# do not separate cleanly between which is a lake or a rivers; and rivers are
# dealth with below (to identify rivers with sufficient current).
# Plus, lakes often have very narrow passages that are not really rivers, but 
# because there is next to no current, animals can easily traverse them. We want
# to remove these too.

# Method to generate the layer in QGIS (this takes some time, especially the split with lines) :
# - Get all of the surface water elements (RH_S) from the hydro regions of the
# Géobase du réseau hydrographique du Québec (GRHQ) (https://www.donneesquebec.ca/recherche/dataset/grhq)
# - Merge it all together, clip it inside the extent of the study area
# - Create a negative -50m buffer in the resulting polygons and dissolve the polygons of the buffer
# - In the resulting buffer, delete any polygon < 5000 in $area
# - Create a 100x100m polygon grid with the extent of the area in QGIS
# - Use the "Split by lines" tool of QGIS to cut the rivers/lakes with the grid
# - Select all of the resulting polygons that are in contact with the remaining
# polygons of the negative buffer, and merge them
# - Save the resulting file.
# Explanation : By cutting the lakes/rivers into small pieces and selecting them
# with the negative buffer, we end up selecting only polygon of waters that are
# large enough (so, not rivers and narrow passages of lakes). If we don't do the
# splitting by lines, then we'll just end up selecting the whole polygons every time.
lakesWithoutNarrowsPolygons = gpd.read_file(r'.\ResistanceCostMapCreation\SourceFiles\Lakes_WithNarrowsAndRiversRemoved.shp', use_arrow=True)
rasterizedLakesCost = RasterizePolygons(lakesWithoutNarrowsPolygons,
                                        folderOfOutputsMaps_MEDIUM + '/template_raster_100x100.tif',
                                        1000,
                                        allTouched=True)
writeNewRasterData(folderOfOutputsMaps_MEDIUM + '/lakesCosts.tif',
                   rasterizedLakesCost, size = "MEDIUM")

#%% AQRESEAU+ AND OTHERS

# RAILS

# From https://www.donneesquebec.ca/recherche/dataset/adresses-quebec, "AQréseau+" in shp
# We load the rail lines
# Shapefile is simply Reseau_ferroviaire.shp from AQRéseau+
# Reprojected into EPSG:32198 - NAD83 / Quebec Lambert
# and clipped in the study area
railLinesShapefile = gpd.read_file(resistanceMapFolder + ".\ResistanceCostMapCreation\SourceFiles\AQReseau_ReseauFerroviaire_ReprojectedClipped.shp", use_arrow=True)
# We rasterize, but with using the "ALL TOUCHED" method to avoid holes in the path
# of the line.
# Cost here is 1000.
rasterizedRailsCost = RasterizePolygons(railLinesShapefile,
                                            folderOfOutputsMaps_MEDIUM + '/template_raster.tif',
                                            1000,
                                            allTouched=True)
writeNewRasterData(folderOfOutputsMaps_MEDIUM + '/railsCosts.tif',
                   rasterizedRailsCost)


# ROADS

# From https://www.donneesquebec.ca/recherche/dataset/adresses-quebec, "AQréseau+" in shp
# We load the roads
# Shapefile is simply Reseau_routier.shp from AQRéseau+
# Reprojected into EPSG:32198 - NAD83 / Quebec Lambert
# and clipped in the study area
roadsShapefile = gpd.read_file(resistanceMapFolder + ".\ResistanceCostMapCreation\SourceFiles\AQReseau_ReseauRoutier_ReprojectedClipped.shp", use_arrow=True)
# We edit the NbrVoies field to express it as numbers, so that we can select things easily afterward
roadsShapefile['NbrVoies'] = roadsShapefile['NbrVoies'].replace('INC', 0)
roadsShapefile['NbrVoies'] = roadsShapefile['NbrVoies'].fillna(0).astype(int)

# First, we do highways with more than 2 lanes or that are highways - cost 1000
highwaysMore2LanesPolygons = roadsShapefile[(roadsShapefile['NbrVoies'] > 2) | (roadsShapefile['NomRte'].str.contains('Autoroute', case=False)) | (roadsShapefile['ClsRte'] == "Autoroute")]
rasterizedHighwaysMore2Lanes = RasterizePolygons(highwaysMore2LanesPolygons,
                                            folderOfOutputsMaps_MEDIUM + '/template_raster.tif',
                                            1000,
                                            allTouched=True)
writeNewRasterData(folderOfOutputsMaps_MEDIUM + '/highwaysMore2LanesCosts.tif',
                   rasterizedHighwaysMore2Lanes)

# Then, we do highways with 2 lanes - cost 100
highways2LanesPolygons = roadsShapefile[roadsShapefile['NbrVoies'] == 2]
rasterizedHighways2Lanes = RasterizePolygons(highways2LanesPolygons,
                                            folderOfOutputsMaps_MEDIUM + '/template_raster.tif',
                                            100,
                                            allTouched=True)
writeNewRasterData(folderOfOutputsMaps_MEDIUM + '/highways2LanesCosts.tif',
                   rasterizedHighways2Lanes)


# Then, we do all of the other roads - cost 10
allOtherRoadsPolygons = roadsShapefile[roadsShapefile['NbrVoies'] < 2]
rasterizedallOtherRoads = RasterizePolygons(allOtherRoadsPolygons,
                                            folderOfOutputsMaps_MEDIUM + '/template_raster.tif',
                                            10,
                                            allTouched=True)
writeNewRasterData(folderOfOutputsMaps_MEDIUM + '/allOtherRoadsCosts.tif',
                   rasterizedallOtherRoads)

# RIVERS

# From https://www.hydrosheds.org/products/hydrorivers, North and Central America Shapefile
# We load the rivers
# Shapefile is HydroRIVERS_v10_na.shp from the HydroRIVERS dataset
# Reprojected into EPSG:32198 - NAD83 / Quebec Lambert
# and clipped in the study area
riversShapefile = gpd.read_file(resistanceMapFolder + ".\ResistanceCostMapCreation\SourceFiles\HydroRIVERS_ReprojectedClipped.shp", use_arrow=True)
# We get all rivers with a sufficient debit
riversWithSufficientDebit = riversShapefile[riversShapefile['DIS_AV_CMS'] >= 28]
# We rasterize with cost = 1000
riversRasterized = RasterizePolygons(riversWithSufficientDebit,
                                            folderOfOutputsMaps_MEDIUM + '/template_raster.tif',
                                            1000,
                                            allTouched=True)
writeNewRasterData(folderOfOutputsMaps_MEDIUM + '/riversCosts.tif',
                   riversRasterized)


# SLOPE

# First, you need to obtain elevation data
# from https://www.donneesquebec.ca/recherche/dataset/modeles-numeriques-d-altitude-a-l-echelle-de-1-20-000
# Advice : Use a browser plugin like "DownThemAll!" in firefox to easily download all of the files,
# and the "load them all" QGIS plugin to load all raster files. However, the files are in .bin format,
# and "Load Them All" can't deal with them. Simply rename the .bin files to .tif. In windows, 
# use a command prompt, put yourself in a folder with all of the subfolders with the elevation with "cd",
# and then use the command "for /r %x in (*.bin) do ren "%x" *.bin". It will rename all of the .tif in
# the folders and subfolders for you. You can then load them easily with "Load Them All" in QGIS, and then
# merge them.
# Then, reproject the merged raster in a m2 projection (e.g. Quebec Lambert),
# and finally use the "Slope" tool in QGIS to compute the slope. Don't use the slope
# tool without m2 reprojection, or you will have absurd slope vlaues.

# The code here does bilinear resampling to adapt the slope raster to the resolution
# we need for the cost raster (100x100m.)
# Code for resampling from Perplexity; works well.
# Bilinear ressampling is well adapted for elevation data.
with rasterio.open(r".\ResistanceCostMapCreation\SourceFiles\slope.tif") as slopeRasterData:
    # We put the cost of slope 30 to 1000, the rest to 1.
    slopeRasterDataCost = np.where(slopeRasterData >= 30, 1000, 1)
    # We save it
    writeNewRasterData(folderOfOutputsMaps_MEDIUM + '/slopeCosts.tif',
                       slopeRasterDataCost)


#%% COMBINING EVERYTHING INTO THE FINAL RASTER

# We load all of the rasters that we've created in the folder 
extension = ".tif"
rasterFiles = []
for file_name in os.listdir(folderOfOutputsMaps_MEDIUM):
    if file_name.endswith(extension) and "resistanceCostMapPitherEtAl2023" not in file_name:
        file_path = os.path.join(folderOfOutputsMaps_MEDIUM, file_name)
        rasterFiles.append(file_path)
        
dictionnaryOfRasterArrays = dict()
for file in rasterFiles:
    dictionnaryOfRasterArrays[file.split("\\")[-1]] = getRasterData(file)
    # In case a raster doesn't have the right shape, the lines after will return an error.
    # This detects the problematic raster.
    # print("Array shape of " + str(file.split("\\")[-1]) + " : " + str(dictionnaryOfRasterArrays[file.split("\\")[-1]].shape))

arrays_list = list(dictionnaryOfRasterArrays.values())
# Combine the arrays using np.maximum.reduce() : Each pixel becomes the maximum
# value for the pixel between all of our rasters, and so the highest cost of movement.
combined_array = np.maximum.reduce(arrays_list)

# We save the combined array.
writeNewRasterData(folderOfOutputsMaps_MEDIUM + '/resistanceCostMapPitherEtAl2023.tif',
                   combined_array)
# Finally, we save it as a .asc raster (Circuitscape will need this in the next script)
convertTIFToASCRaster(folderOfOutputsMaps_MEDIUM + '/resistanceCostMapPitherEtAl2023.tif')