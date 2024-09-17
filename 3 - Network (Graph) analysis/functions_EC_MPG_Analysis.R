#####################################################
#
# @author: Clement Hardy (clem.hardy@pm.me)
#
# This script contains custom functions to compute the
# probability of connectivity (PC) and Equivalent Connectivity
# (EC) based on a Grainscape object.
#
# See https://www.alexchubaty.com/grainscape/articles/grainscape_vignette.html
# for tutorial about grainscape
#
# WARNING : In previous code, Igraph generated with Grainscape had $patchID and
# $patchArea for the ID and area of patches. Here, it's $name and $Area. I changed
# the code of the functions to adapt to that.
#####################################################



#### LOADING PACKAGES ####

library(grainscape)
library(igraph)
library(raster)
library(knitr)
library(ggplot2)
library(sf)
library(Dict)
library(rlang)
library(pbapply)

#### FUNCTIONS ####

## Adds a attribute to a shapefile
# List of values need to be the same length as the features in the shapefile,
# and must match the order of the features.
# Used to add attributes to the files exported by Grainscape.
# Courtesy of Perplexity (https://www.perplexity.ai/search/create-a-r-function-that-adds-Ay2eIdUaQPO4aVcUdaZymw)
add_attribute_to_shapefile <- function(shapefile, new_attribute_name, new_attribute_values) {
  # Read the shapefile
  shp <- st_read(shapefile)
  
  # Check if the length of new_attribute_values matches the number of features in the shapefile
  if (length(new_attribute_values) != nrow(shp)) {
    stop("The length of new_attribute_values must match the number of features in the shapefile")
  }
  
  # Add the new attribute to the shapefile
  shp[[new_attribute_name]] <- new_attribute_values
  
  # Write the updated shapefile
  output_filename <- paste0(tools::file_path_sans_ext(shapefile), ".shp")
  st_write(shp, output_filename, driver = "ESRI Shapefile", append = FALSE)
  
  cat("Updated shapefile saved as:", output_filename, "\n")
}

## Compute and add BETWEENNESS CENTRALITY INDEX
# The function computes the Betweennees centrality index on a grainscape MPG
# object, and then saves the values in the shapefile exported by grainscape.
addBetweennessToShapefile <- function(IgraphObject,
                                       nodeShapefilePath) {
  mpg_graph <- IgraphObject
  ig <- graph_from_adjacency_matrix(as_adjacency_matrix(mpg_graph), mode = "undirected")
  betweenness_scores <- betweenness(ig)
  add_attribute_to_shapefile(nodeShapefilePath,
                             "Betweenness",
                             betweenness_scores)
}

## Computes a probability of mouvement from a resistance cost associated with
# a path between two nodes
# It's arbitrary.
# See https://sci-hub.st/https://www.sciencedirect.com/science/article/abs/pii/S0169204607000709?via%3Dihub :
# for a similar application.
# Might be updated with a more complex function : https://hal.science/hal-01681596/file/16-671albert_v2%20HAL.pdf
inverse_exponential_probability <- function(distance) {
  # Define the maximum distance
  max_distance <- 10000
  
  # Calculate the probability
  probability <- exp(-distance / max_distance)
  
  return(probability)
}

## EXPORTING LINK RESISTANCE COST VALUE ###
# Computes the probability of mouvement for each link of the network
# and exports it as an attribute to the shapefile
# WARNING : Make sure the file is not opened in QGIS, or it will be corrupted !
add_ProbabilityOfMovement_toShapefile <- function(linksShapefilePath){
  # Read the shapefile
  shp <- st_read(linksShapefilePath)
  
  # Add the new attribute to the shapefile
  shp[["probMouv"]] <- sapply(shp[["lcpPerWt"]], function(x) inverse_exponential_probability(x))
  
  # Write the updated shapefile
  st_write(shp, linksShapefilePath, driver = "ESRI Shapefile", append = FALSE)
  
  cat("Updated shapefile saved as:", linksShapefilePath, "\n")
}


## Computes the PC of the network.
# Takes a Grainscape MPG object and the area of the landscape
# If returnFractions is false, returns a single value of PC
# If true, returns a list with the value of PC, PC_intra, PC_direct and PC_step
computePCOfNetwork <- function(IgraphObject,
                               landscapeAreaM2,
                               returnFractions = TRUE) {
  
  ## Generate a matrix indicating the number of links to go from one node to the other
  matrixOfEdgeNumbersBetweenNodes = distances(
    IgraphObject,
    v = V(IgraphObject),
    to = V(IgraphObject),
    mode = c("all"),
    algorithm = c("unweighted")
  )
  ## We can reclassify this matrix to remove pairs from the final matrix product
  matrixAllPairs <- ifelse(matrixOfEdgeNumbersBetweenNodes == 0, 1, 1) # All are set to 1
  # If we're going to compute the fractions, then we'll prepare the matrices for their
  # computation
  if (returnFractions){
    matrixPairsIntra <- ifelse(matrixOfEdgeNumbersBetweenNodes == 0, 1, 0)
    matrixPairsDirect <- ifelse(matrixOfEdgeNumbersBetweenNodes == 1, 1, 0)
    matrixPairsStep <- ifelse(matrixOfEdgeNumbersBetweenNodes > 1, 1, 0)
  }
  
  ## Making a matrix of patch area
  # Associate the patch ID and patch Area in a temporary dataframe
  patch_df <- data.frame(
    patch_id = V(IgraphObject)$name,
    area = V(IgraphObject)$Area
  )
  # Create the matrix object (allocate memory)
  matrixPatchAreaProduct <- matrix(0, nrow = nrow(patch_df), ncol = nrow(patch_df))
  # Fill the matrix object with the values
  for (i in 1:nrow(patch_df)) {
    for (j in 1:nrow(patch_df)) {
      matrixPatchAreaProduct[i, j] <- patch_df$area[i] * patch_df$area[j]
    }
  }
  # Name the columns of the matrix with the patch ID for easy calling
  rownames(matrixPatchAreaProduct) <- patch_df$patch_id
  colnames(matrixPatchAreaProduct) <- patch_df$patch_id
  
  # We compute the shortest paths
  matrixOfShortestPaths = distances(
    IgraphObject,
    v = V(IgraphObject),
    to = V(IgraphObject),
    mode = c("all"),
    weights = E(IgraphObject)$lcpPerimWeight,
    algorithm = c("dijkstra")
  )

  # We transform the shortest paths (sum of resistance along links) into
  # max probabilities of movement through our negative exponential (see above)
  matrixOfMaxProbMovement <- apply(matrixOfShortestPaths, c(1, 2), inverse_exponential_probability)
  
  # Now, we compute the above component of the PC
  # Product of 3 matrix :
  # - The one selecting the pairs of patches to consider (0 or 1, puts every pairs
  # not of interest to 0, hence a product of 0 for this cell. Here, we take all pairs.)
  # - The one indicating the product of patch areas for the pair (self-product if its
  # a patch with itself)
  # - The one indicating the max prob of mouvement between the two patches (1.0 if patch
  # with itself, 0 if the patches are not connected because different components)
  Sum_PC_Above = sum(matrixAllPairs * matrixPatchAreaProduct * matrixOfMaxProbMovement)
  
  # We divide it by landscape area squared to finally get the PC
  Initial_PC = Sum_PC_Above/landscapeAreaM2**2
  
  if (!returnFractions){
    return(Initial_PC)
  }
  
  if (returnFractions){
    # We do the same but with PC_intra : here, the pair lists are only for the patches with themselves
    # This gives us all of the terms of the sum
    Sum_PC_intra_Above = sum(matrixPairsIntra * matrixPatchAreaProduct * matrixOfMaxProbMovement)
    # We divide it by landscape area squared to finally get the PC_intra
    Initial_PC_intra = Sum_PC_intra_Above/landscapeAreaM2**2
    
    # We do the same but with PC_direct, only for pairs that have a direct connection
    Sum_PC_direct_Above = sum(matrixPairsDirect * matrixPatchAreaProduct * matrixOfMaxProbMovement)
    # We divide it by landscape area squared to finally get the PC_intra
    Initial_PC_direct = Sum_PC_direct_Above/landscapeAreaM2**2
    
    # Finally, we compute dPC_step
    Sum_PC_step_Above = sum(matrixPairsStep * matrixPatchAreaProduct * matrixOfMaxProbMovement)
    Initial_PC_step = Sum_PC_step_Above/landscapeAreaM2**2
    
    # If we return the fractions, we return them as a list
    return(list(
      PC = Initial_PC,
      PC_intra = Initial_PC_intra,
      PC_direct = Initial_PC_direct,
      PC_step = Initial_PC_step
    ))
  }
}

## Computes the dPC for each node of the network, and add them to the shapefile
# of nodes exported by grainscape
add_dPC_toShapefile <- function(IgraphObject,
                                landscapeAreaM2,
                                nodeShapefilePath){
  
  # We compute the initial PC
  initial_PC_andFraction = computePCOfNetwork(IgraphObject, landscapeAreaM2, returnFractions = TRUE)
  
  # We create the dict that will associate a dPC to a node removal
  dictOfDPCvalues = list()
  # We also create the dict that will contain dPC_intra values, meaning
  # only for when we consider the connectivity of a patch with itself (its area)
  dictOfDPC_intra_values = list()
  # We also create the dict that will contain dPC_direct, meaning when we consider only
  # the connectivity to/from a patch in direction connection cases (one step)
  dictOfDPC_direct_values = list()
  # We  make a dictionnary for dPC_step, meaning when we consider only the
  # connectivity to/from a patch with multile steps/stepping stones.
  # We can derive it from the 3 previous values, as dPC = dPC_intra + dPC_direct + dPC_step
  dictOfDPC_step_values = list()
  
  # We initialize a progress bar
  pb <- txtProgressBar(min = 0, max = length(V(IgraphObject)$name), style = 3)
  progress = 0
  
  for (patchID in V(IgraphObject)$name){
    # print("Removing node")
    # We remove the node
    graphWithNodeRemoved = delete_vertices(IgraphObject, which(V(IgraphObject)$name == patchID))
    
    # We recompute the matrices
    # matrices of pairs
    matrixOfEdgeNumbersBetweenNodes = distances(
      graphWithNodeRemoved,
      v = V(graphWithNodeRemoved),
      to = V(graphWithNodeRemoved),
      mode = c("all"),
      algorithm = c("unweighted")
    )
    # we reconpute the matrices of pairs
    matrixAllPairs <- ifelse(matrixOfEdgeNumbersBetweenNodes == 0, 1, 1) # All are set to 1
    matrixPairsIntra <- ifelse(matrixOfEdgeNumbersBetweenNodes == 0, 1, 0)
    matrixPairsDirect <- ifelse(matrixOfEdgeNumbersBetweenNodes == 1, 1, 0)
    matrixPairsStep <- ifelse(matrixOfEdgeNumbersBetweenNodes > 1, 1, 0)
    # We recompute the matrix of patch areas
    patch_df <- data.frame(
      patch_id = V(graphWithNodeRemoved)$name,
      area = V(graphWithNodeRemoved)$Area
    )
    matrixPatchAreaProduct <- matrix(0, nrow = nrow(patch_df), ncol = nrow(patch_df))
    for (i in 1:nrow(patch_df)) {
      for (j in 1:nrow(patch_df)) {
        matrixPatchAreaProduct[i, j] <- patch_df$area[i] * patch_df$area[j]
      }
    }
    rownames(matrixPatchAreaProduct) <- patch_df$patch_id
    colnames(matrixPatchAreaProduct) <- patch_df$patch_id
    # We recompute the Matrix of maxprob
    matrixOfShortestPaths = distances(
      graphWithNodeRemoved,
      v = V(graphWithNodeRemoved),
      to = V(graphWithNodeRemoved),
      mode = c("all"),
      weights = E(graphWithNodeRemoved)$lcpPerimWeight,
      algorithm = c("dijkstra")
    )
    matrixOfMaxProbMovement <- apply(matrixOfShortestPaths, c(1, 2), inverse_exponential_probability)
    
    # We re-compute PC in the new graph and with the new dict
    # print("Computing dPC")
    Sum_PC_Above_nodeRemoved = sum(matrixAllPairs * matrixPatchAreaProduct * matrixOfMaxProbMovement)
    # We divide it by landscape area squared to finally get the PC
    PC_nodeRemoved = Sum_PC_Above_nodeRemoved/landscapeAreaM2**2
    # We compute the difference to the initial PC
    dPC_nodeRemoved = (100 * (initial_PC_andFraction$PC - PC_nodeRemoved))/initial_PC_andFraction$PC 
    # We save the dPC in the list
    dictOfDPCvalues[[patchID]] = dPC_nodeRemoved
    
    # print("Computing dPC_intra")
    # We do the same but with dPC_intra : here, the pair lists are only for the patches with themselves
    # This gives us all of the terms of the sum
    Sum_PC_intra_Above_nodeRemoved = sum(matrixPairsIntra * matrixPatchAreaProduct * matrixOfMaxProbMovement)
    PC_intra_nodeRemoved = Sum_PC_intra_Above_nodeRemoved/landscapeAreaM2**2
    dPC_intra_nodeRemoved = (100 * (initial_PC_andFraction$PC_intra - PC_intra_nodeRemoved))/initial_PC_andFraction$PC_intra
    dictOfDPC_intra_values[[patchID]] = dPC_intra_nodeRemoved
    
    # print("Computing dPC_direct")
    # We do the same but with dPC_direct : here, we only take pairs were the shortest path between the
    # nodes was direct (no stepping stone)
    Sum_PC_direct_Above_nodeRemoved = sum(matrixPairsDirect * matrixPatchAreaProduct * matrixOfMaxProbMovement)
    # We divide it by landscape area squared to finally get the PC
    PC_direct_nodeRemoved = Sum_PC_direct_Above_nodeRemoved/landscapeAreaM2**2
    dPC_direct_nodeRemoved = (100 * (initial_PC_andFraction$PC_direct - PC_direct_nodeRemoved))/initial_PC_andFraction$PC_direct
    dictOfDPC_direct_values[[patchID]] = dPC_direct_nodeRemoved
    
    # We finish with dPC_step
    Sum_PC_direct_Above_nodeRemoved = sum(matrixPairsStep * matrixPatchAreaProduct * matrixOfMaxProbMovement)
    # We divide it by landscape area squared to finally get the PC
    PC_step_nodeRemoved = Sum_PC_direct_Above_nodeRemoved/landscapeAreaM2**2
    dPC_step_nodeRemoved = (100 * (initial_PC_andFraction$PC_step - PC_step_nodeRemoved))/initial_PC_andFraction$PC_step
    dictOfDPC_step_values[[patchID]] = dPC_step_nodeRemoved
    
    # Updates the progress bar
    progress = progress + 1
    setTxtProgressBar(pb, progress)
  }
  
  ### Exporting everything to the .csv file ###
  dPCRawValues <- unlist(dictOfDPCvalues, use.names = FALSE)
  
  add_attribute_to_shapefile(nodeShapefilePath,
                             "dPC",
                             dPCRawValues)
  
  dPC_intra_RawValues <- unlist(dictOfDPC_intra_values, use.names = FALSE)
  
  add_attribute_to_shapefile(nodeShapefilePath,
                             "dPC_intra",
                             dPC_intra_RawValues)
  
  dPC_direct_RawValues <- unlist(dictOfDPC_direct_values, use.names = FALSE)
  
  add_attribute_to_shapefile(nodeShapefilePath,
                             "dPC_direct",
                             dPC_direct_RawValues)
  
  dPC_step_RawValues <- unlist(dictOfDPC_step_values, use.names = FALSE)
  
  add_attribute_to_shapefile(nodeShapefilePath,
                             "dPC_step",
                             dPC_step_RawValues)
  
}

## Compute the Equivalent Connectivity for the network
computeECOfNetwork <- function(IgraphObject,
                               landscapeAreaM2) {
  
  ## Generate a matrix indicating the number of links to go from one node to the other
  matrixOfEdgeNumbersBetweenNodes = distances(
    IgraphObject,
    v = V(IgraphObject),
    to = V(IgraphObject),
    mode = c("all"),
    algorithm = c("unweighted")
  )
  ## We can reclassify this matrix to remove pairs from the final matrix product
  matrixAllPairs <- ifelse(matrixOfEdgeNumbersBetweenNodes == 0, 1, 1) # All are set to 1
  # If we're going to compute the fractions, then we'll prepare the matrices for their
  # computation
  
  ## Making a matrix of patch area
  # Associate the patch ID and patch Area in a temporary dataframe
  patch_df <- data.frame(
    patch_id = V(IgraphObject)$name,
    area = V(IgraphObject)$Area
  )
  # Create the matrix object (allocate memory)
  matrixPatchAreaProduct <- matrix(0, nrow = nrow(patch_df), ncol = nrow(patch_df))
  # Fill the matrix object with the values
  for (i in 1:nrow(patch_df)) {
    for (j in 1:nrow(patch_df)) {
      matrixPatchAreaProduct[i, j] <- patch_df$area[i] * patch_df$area[j]
    }
  }
  # Name the columns of the matrix with the patch ID for easy calling
  rownames(matrixPatchAreaProduct) <- patch_df$patch_id
  colnames(matrixPatchAreaProduct) <- patch_df$patch_id
  
  # We compute the shortest paths
  matrixOfShortestPaths = distances(
    IgraphObject,
    v = V(IgraphObject),
    to = V(IgraphObject),
    mode = c("all"),
    weights = E(IgraphObject)$lcpPerimWeight,
    algorithm = c("dijkstra")
  )
  
  # We transform the shortest paths (sum of resistance along links) into
  # max probabilities of movement through our negative exponential (see above)
  matrixOfMaxProbMovement <- apply(matrixOfShortestPaths, c(1, 2), inverse_exponential_probability)
  
  # Now, we compute the above component of the PC
  # Product of 3 matrix :
  # - The one selecting the pairs of patches to consider (0 or 1, puts every pairs
  # not of interest to 0, hence a product of 0 for this cell. Here, we take all pairs.)
  # - The one indicating the product of patch areas for the pair (self-product if its
  # a patch with itself)
  # - The one indicating the max prob of mouvement between the two patches (1.0 if patch
  # with itself, 0 if the patches are not connected because different components)
  Sum_EC_BeforeSqrt = sum(matrixAllPairs * matrixPatchAreaProduct * matrixOfMaxProbMovement)
  
  # We do the square root of it to get the EC (or ECA)
  EC_network = sqrt(Sum_EC_BeforeSqrt)
  
  return(EC_network)
}

#### Compute PC and EC based on a resistance map and a shapefile layer representing
# patches of habitats.
# Uses Graphab.
# Returns a list with PC, its fractions, and EC
computePC_and_EC_graphab <- function(listOfShapefileLayers,
                                     resistanceMapPath,
                                     temporaryRasterPath,
                                     graphabProjectName){
  start_time <- Sys.time()
  
  print("Loading resistance raster")
  # We load the resistance raster
  resistanceMap <- raster(resistanceMapPath)
  
  # Creating the raster for the patches
  # Mixing sacred mountains and all protected areas
  # Merging the shapefiles
  print("Merging shapefiles")
  merged_shapefile <- bind_rows(listOfShapefileLayers)
  # Rasterizing them, taking into account any cell touched by a polygon
  # Using the terra package for touches = TRUE (seem to approximate GDAL's ALL_TOUCHED;
  # see https://rdrr.io/cran/terra/man/rasterize.html)
  print("Rasterizing shapefiles")
  rasterized <- rasterize(merged_shapefile, resistanceMap,
                          touches = TRUE, field = 1,
                          filename = temporaryRasterPath,
                          overwrite = TRUE,
                          datatype="INT4S")
  # plot(rasterized)
  # Writing temporary raster
  # writeRaster(rasterized, filename=temporaryRasterPath,
  # format="GTiff", datatype="INT4S", overwrite=TRUE)
  
  # plot(rasterized)
  print("Creating Graphab project")
  # We create the graphab project
  graphab_project(proj_name = graphabProjectName,
                  temporaryRasterPath,
                  habitat = c(1),
                  minarea = 0)
  
  print("Creating Graphab links")
  # WE create the shortest paths for the complete graph
  # Around 2 minutes to create all the links for the network
  graphab_link(proj_name = graphabProjectName,
               distance = "cost",
               cost = resistanceMapPath,
               name = "complete_links_LeastResistance",
               topo = "complete",
               remcrosspath = FALSE)
  
  print("Creating iGraph object")
  # We create an igraph object from the graphab project so that I can use custom
  # functions to compute PC and EC for the network.
  igraphGraph =  graphab_to_igraph(
    proj_name = graphabProjectName,
    linkset = "complete_links_LeastResistance",
    nodes = "patches",
    weight = "cost",
    fig = FALSE,
    crds = FALSE
  )
  
  print("Computing PC and EC")
  listOfResults = list()
  # We compute PC
  # As we will see, PC_step = 0. It's because since we're in a complete graph,
  # All shortest paths are of 1 edge (since every patch is directly connected to
  # one another through a shortest path).
  resultsPC = computePCOfNetwork(igraphGraph, landscapeAreaM2)
  listOfResults["PC"] = resultsPC["PC"]
  listOfResults["PC_intra"] = resultsPC["PC_intra"]
  listOfResults["PC_direct"] = resultsPC["PC_direct"]
  listOfResults["PC_step"] = resultsPC["PC_step"]
  # We compute EC
  listOfResults["EC"] = computeECOfNetwork(igraphGraph, landscapeAreaM2)
  
  end_time <- Sys.time()
  print("Time taken :")
  print(end_time - start_time)
  
  return(listOfResults)
}
  

#### Compute dPC and betweenness based on a resistance map and a raster layer representing
# patches of habitats.
#
# WARNING : This function uses functions coming from the graph4lg package that use
# the "graphab" java program (made by the same authors). It's very useful to compute
# the complete graph for the different scenarios without takig too much time.
#
# However, I had a couple
# of errors when first trying to use them. 
# For these graphab functions to work, needed 2 things :
# 1. Be in a folder with no spaces or complex chracters (' or others) in the path
# 2. Replace the file C:\Users\Clement\AppData\Local\graph4lg_jar\graphab-2.8.jar
# downloaded by the graph4lg R package when it was installed
# with a .jar downloaded from the Graphab website and renamed graphab-2.8.jar
# as what I had seemed to be corrupted.
#
#
# The function will create a folder with shapefiles (graphabProjectName) in the working directory
# that will have dPC and Betweeness associated to the nodes
# Can choose werever an MPG or a complete graph is created
computedPC_and_Betweenness_graphab <- function(patchRasterPath,
                                     resistanceMapPath,
                                     graphabProjectName,
                                     MPG = FALSE){
  start_time <- Sys.time()
  
  # Landscape area : extent of the cost raster
  resistanceMap <- raster(resistanceMapPath)
  landscapeAreaM2 = res(resistanceMap)[1]^2 * (ncell(resistanceMap))
  
  # plot(rasterized)
  print("Creating Graphab project")
  # We create the graphab project
  graphab_project(proj_name = graphabProjectName,
                  patchRasterPath,
                  habitat = c(1),
                  minarea = 0)
  
  print("Creating Graphab links")
  # WE create the shortest paths for the complete graph
  # Around 2 minutes to create all the links for the network
  if (!MPG)
  {
    graphab_link(proj_name = graphabProjectName,
                 distance = "cost",
                 cost = resistanceMapPath,
                 name = "linkset",
                 topo = "complete",
                 remcrosspath = FALSE)
  }
  else
  {
    graphab_link(proj_name = graphabProjectName,
                 distance = "cost",
                 cost = resistanceMapPath,
                 name = "linkset",
                 topo = "planar",
                 remcrosspath = TRUE)
  }

  
  print("Creating iGraph object")
  # We create an igraph object from the graphab project so that I can use custom
  # functions to compute PC and EC for the network.
  igraphGraph =  graphab_to_igraph(
    proj_name = graphabProjectName,
    linkset = "linkset",
    nodes = "patches",
    weight = "cost",
    fig = FALSE,
    crds = FALSE
  )
  
  
  
  print("Computing dPC and betweenness")
  add_dPC_toShapefile(igraphGraph,
                      landscapeArea,
                      paste("./", graphabProjectName, "/patches.shp", sep = ""))
  
  addBetweennessToShapefile(igraphGraph,
                            paste("./", graphabProjectName, "/patches.shp", sep = ""))
  
  print("Everything is done !")
  print("Attribute of shapefile patches.shp :")
  shapefile <- st_read(paste("./", graphabProjectName, "/patches.shp", sep = ""))
  print(shapefile)
}


