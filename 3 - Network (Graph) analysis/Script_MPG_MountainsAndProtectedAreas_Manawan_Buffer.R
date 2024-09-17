#####################################################
#
# @author: Clement Hardy (clem.hardy@pm.me)
#
# This script generates a Minimal Planar Graph (MPG) for the
# a set of protected areas, and a resistance map that has been
# created by a previous script (see CreatingResistanceMap_PitherEtAl2023.py). 
#
# The MPG is created using the grainscape package.
# It is then exported as a shapefile and some rasters.
# Then, the Probability of Connectivity (PC) is calculted for the whole network.
# It is then re-calculated by removing one node of the network each time,
# allowing us to calculate the dPC of the node. The fractions of the dPC (intra, direct, step)
# are also calculated.
# For more information about the PC, what it means and how it is measured and what these fractions of dPC are, see
# https://besjournals.onlinelibrary.wiley.com/doi/10.1111/2041-210X.12799 for examples and equations.
#####################################################

#### LOADING PACKAGES ####

# install.packages("grainscape")
library(grainscape)
library(igraph)
library(raster)
library(knitr)
library(ggplot2)
library(sf)
library(Dict)
library(rlang)
library(pbapply)


#### LOADING DATA ####

# Resistance map (see CreatingResistanceMap_PitherEtAl2023.py for details)

directoryWithFiles = "./Networks_Analysis"

resistanceMap <- raster("./ResistanceCostMapCreation/CostMaps/MEDIUM (100x100m)/resistanceCostMapPitherEtAl2023.tif")

# Patches raster map - This map should contain the location of your protected area of other patches of interest.
# Their position should be indicating with 1, and the rest of the map should be 0.

patchs_SacredMountainsAndProtectedAreas <- raster(paste(directoryWithFiles, "raster_ProtectedAreas.tif", sep = ""))

# Combining the resistance and patch maps for the grainscape functions :
# Patches will have the value 99999. 

resistanceMapAndPatches = resistanceMap
resistanceMapAndPatches[patchs_SacredMountainsAndProtectedAreas == 1] <- 99999

# To check if all is OK

plot(resistanceMapAndPatches)



#### CREATING AND DISPLAYING THE MPG ####

patchyMPG <- MPG(resistanceMapAndPatches, patch = resistanceMapAndPatches == 99999)


plot(patchyMPG, quick = "mpgPlot", theme = FALSE)

# Exports the data as shapefiles for visualisation
pathOfExport = "./Networks_Analysis/"
directoryOfExport = "MPG_Outputs"
fullPathOfExport = paste(pathOfExport, directoryOfExport, sep = "")
export(patchyMPG, dirname = directoryOfExport,
       path = pathOfExport,
       vorBound = FALSE, overwrite = TRUE)

# Since the export function does not currently like it when we add attributes
# (see https://github.com/achubaty/grainscape/issues/71),
# to the nodes, here is a function to add the attributs to nodes or links
# once the shapefiles have been created
nodeShapefilePath = paste(pathOfExport,
                          directoryOfExport,
                          "nodes.shp",
                          sep = "/")
linksShapefilePath = paste(pathOfExport,
                          directoryOfExport,
                          "linksCentroid.shp",
                          sep = "/")

# Function to easily add new attribus to the polygon files
# Will be used to add an attribute with the value of dPC later.
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

#### BETWEENNESS CENTRALITY INDEX
# Computes the betweeness centrality index (BCI) for the polygons, and add it as an attribute.
# For more info on the BCI, see https://link.springer.com/article/10.1007/s10980-006-9015-0 .
# WARNING : Make sure the files are not opened in QGIS, or they will be corrupted !
mpg_graph <- patchyMPG$mpg
ig <- graph_from_adjacency_matrix(as_adjacency_matrix(mpg_graph), mode = "undirected")
betweenness_scores <- betweenness(ig)
add_attribute_to_shapefile(nodeShapefilePath,
                           "Betweenness",
                           betweenness_scores)




#### dPC (CHANGE IN PC WHEN ANY NODE IS TAKEN AWAY) ####

### PREPARATIONS BEFORE THE COMPUTATION ###

## Generate a matrix indicating the number of links to go from one node to the other
matrixOfEdgeNumbersBetweenNodes = distances(
                                  patchyMPG$mpg,
                                  v = V(patchyMPG$mpg),
                                  to = V(patchyMPG$mpg),
                                  mode = c("all"),
                                  algorithm = c("unweighted")
                                )
## We can reclassify this matrix to remove pairs from the final matrix product
# These different matrices are used to select the pairs of nodes we consider in the different
# computations later. All pairs takes...all pairs; Pairs intra only take one node with itself;
# Direct takes only the nodes that have one connection between one and the other (they are directly
# connected); Step are the connection between two nodes that use multiple intermediary nodes.
matrixAllPairs <- ifelse(matrixOfEdgeNumbersBetweenNodes == 0, 1, 1) # All are set to 1
matrixPairsIntra <- ifelse(matrixOfEdgeNumbersBetweenNodes == 0, 1, 0)
matrixPairsDirect <- ifelse(matrixOfEdgeNumbersBetweenNodes == 1, 1, 0)
matrixPairsStep <- ifelse(matrixOfEdgeNumbersBetweenNodes > 1, 1, 0)


### Define probability of movement associated with each link ###
# In the MPG we created, every link between nodes is a least cost-path between the two nodes
# based on our resistance map.
# The "Weight" of the link is therefore the sum of the resistance along the least-cost path.
# But to compute the PC, we need to transform this into a probability of connectivity.
# To my knowledge, this is a pretty arbitrary and tricky part. Many methods exists.
# Here, I'll use a simple one : a negative exponential function that decrease severely after
# a certain threshold.
# See https://www.sciencedirect.com/science/article/abs/pii/S0169204607000709
# for an exemple in the litterature.
# distance : Distance close to 0 = 1, 50 000= close to 0.
# I recommend you look at this to make a more nuanced version, based on a better parameter : see https://hal.science/hal-01681596/file/16-671albert_v2%20HAL.pdf .
# Remember that the values of resistance used in the resistance map we generated with CreatingResistanceMap_PitherEtAl2023.py are pretty arbitrary and vary
# logarithmically (1, 10, 100 or 1000 depending on the obstacle).
inverse_exponential_probability <- function(distance) {
  # Define the maximum cumulative resistance distance before which the probability drops severely
  # /!\ YOU SHOULD CHANGE THIS VALUE ACCORDING TO YOUR SPECIES OR RESISTANCE MAP
  max_distance <- 10000
  
  # Calculate the probability
  probability <- exp(-distance / max_distance)
  
  return(probability)
}

### Define total landscape area (in same unit as patches area) ###

# Needed to compute the PC, as PC is the probability that two points at random
# in the landscape are in the same patch, or in patches connected together.
# Landscape area : extent of the cost raster
landscapeAreaM2 = res(resistanceMap)[1]^2 * (ncell(resistanceMap))

## Making a matrix of patch area
# Associate the patch ID and patch Area in a temporary dataframe
patch_df <- data.frame(
  patch_id = V(patchyMPG$mpg)$patchId,
  area = V(patchyMPG$mpg)$patchArea
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

# We make dictionnaries that give us the ID or name of the nodes.
# This is because there is a distinction between the two in in the igraph
# object that contains the nodes, and we'll need one or the other at different point In
# the coming functions.
dictNameToID = list()
dictIDToName = list()
for (patchID in V(patchyMPG$mpg)$patchId){
  dictNameToID[[patchID]] = which(V(patchyMPG$mpg)$name == patchID)
  dictIDToName[[which(V(patchyMPG$mpg)$name == patchID)]] = patchID
}


### Calculating the initial PC (with all patches) ###

# We compute the shortest paths
matrixOfShortestPaths = distances(
                        patchyMPG$mpg,
                        v = V(patchyMPG$mpg),
                        to = V(patchyMPG$mpg),
                        mode = c("all"),
                        weights = E(patchyMPG$mpg)$lcpPerimWeight,
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



### CALCULATING dPC BY REMOVE NODES ###

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

# We initialize a progress bar.
pb <- txtProgressBar(min = 0, max = length(V(patchyMPG$mpg)$patchId), style = 3)
progress = 0

for (patchID in V(patchyMPG$mpg)$patchId){
  # print("Removing node")
  # We remove the node
  graphWithNodeRemoved = delete_vertices(patchyMPG$mpg, dictNameToID[[patchID]])
  
  # We recompute the matrices in the new graph (without the removed node)
  # matrices of pairs
  matrixOfEdgeNumbersBetweenNodes = distances(
    graphWithNodeRemoved,
    v = V(graphWithNodeRemoved),
    to = V(graphWithNodeRemoved),
    mode = c("all"),
    algorithm = c("unweighted")
  )
  matrixAllPairs <- ifelse(matrixOfEdgeNumbersBetweenNodes == 0, 1, 1) # All are set to 1
  matrixPairsIntra <- ifelse(matrixOfEdgeNumbersBetweenNodes == 0, 1, 0)
  matrixPairsDirect <- ifelse(matrixOfEdgeNumbersBetweenNodes == 1, 1, 0)
  matrixPairsStep <- ifelse(matrixOfEdgeNumbersBetweenNodes > 1, 1, 0)
  # We recompute the matrix of patch areas
  patch_df <- data.frame(
    patch_id = V(graphWithNodeRemoved)$patchId,
    area = V(graphWithNodeRemoved)$patchArea
  )
  matrixPatchAreaProduct <- matrix(0, nrow = nrow(patch_df), ncol = nrow(patch_df))
  for (i in 1:nrow(patch_df)) {
    for (j in 1:nrow(patch_df)) {
      matrixPatchAreaProduct[i, j] <- patch_df$area[i] * patch_df$area[j]
    }
  }
  rownames(matrixPatchAreaProduct) <- patch_df$patch_id
  colnames(matrixPatchAreaProduct) <- patch_df$patch_id
  # We also redo the matrix of probability of mouvement on the paths of the nodes
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
  dPC_nodeRemoved = (100 * (Initial_PC - PC_nodeRemoved))/Initial_PC
  # We save the dPC in the list
  dictOfDPCvalues[[patchID]] = dPC_nodeRemoved
  
  # print("Computing dPC_intra")
  # We do the same but with dPC_intra : here, the pair lists are only for the patches with themselves
  # This gives us all of the terms of the sum
  Sum_PC_intra_Above_nodeRemoved = sum(matrixPairsIntra * matrixPatchAreaProduct * matrixOfMaxProbMovement)
  PC_intra_nodeRemoved = Sum_PC_intra_Above_nodeRemoved/landscapeAreaM2**2
  dPC_intra_nodeRemoved = (100 * (Initial_PC_intra - PC_intra_nodeRemoved))/Initial_PC_intra
  dictOfDPC_intra_values[[patchID]] = dPC_intra_nodeRemoved
  
  # print("Computing dPC_direct")
  # We do the same but with dPC_direct : here, we only take pairs were the shortest path between the
  # nodes was direct (no stepping stone)
  Sum_PC_direct_Above_nodeRemoved = sum(matrixPairsDirect * matrixPatchAreaProduct * matrixOfMaxProbMovement)
  # We divide it by landscape area squared to finally get the PC
  PC_direct_nodeRemoved = Sum_PC_direct_Above_nodeRemoved/landscapeAreaM2**2
  dPC_direct_nodeRemoved = (100 * (Initial_PC_direct - PC_direct_nodeRemoved))/Initial_PC_direct
  dictOfDPC_direct_values[[patchID]] = dPC_direct_nodeRemoved
  
  # We finish with dPC_step
  Sum_PC_direct_Above_nodeRemoved = sum(matrixPairsStep * matrixPatchAreaProduct * matrixOfMaxProbMovement)
  # We divide it by landscape area squared to finally get the PC
  PC_step_nodeRemoved = Sum_PC_direct_Above_nodeRemoved/landscapeAreaM2**2
  dPC_step_nodeRemoved = (100 * (Initial_PC_step - PC_step_nodeRemoved))/Initial_PC_step
  dictOfDPC_step_values[[patchID]] = dPC_step_nodeRemoved
  
  # Updates the progress bar
  progress = progress + 1
  setTxtProgressBar(pb, progress)
}


### Exporting everything to the .csv file and as an attribute in the shapefile ###
# WARNING : Make sure the files are not opened in QGIS, or they will be corrupted !
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


### EXPORTING LINK RESISTANCE COST VALUE ###

# We compute the probability for each link and export it to the shapefile
# WARNING : Make sure the file is not opened in QGIS, or it will be corrupted !

# Read the shapefile
shp <- st_read(linksShapefilePath)

# Add the new attribute to the shapefile
shp[["probMouv"]] <- sapply(shp[["lcpPerWt"]], function(x) inverse_exponential_probability(x))

# Write the updated shapefile
st_write(shp, linksShapefilePath, driver = "ESRI Shapefile", append = FALSE)

cat("Updated shapefile saved as:", linksShapefilePath, "\n")


#### END REMARK ####

# You can now use QGIS to visualize your MPG
# The file linksCentroid.shp contains the links going from node centroid to node centroid,
# with the link probability of mouvement as attributes
# The file nodes.shp contains the position of the nodes and their Area, BCI and dPC + fractions associated to them.







































#### DRAFTS ####

generate_unique_pairs <- function(listID) {
  # Get all combinations of pairs
  pairs <- combn(listID, 2)
  
  # Convert to a list of vectors (each vector represents a pair)
  pair_list <- lapply(seq_len(ncol(pairs)), function(i) pairs[,i])
  
  return(pair_list)
}

pairsOfPatches = generate_unique_pairs(V(patchyMPG$mpg)$patchId)

### Make a dictionnary containing these probabilities by pair of patch ###
# for easier access latter

# Make a dataframe with each row for an edge containing their ID, origin/end points and weight
# We also use the function defined earlier to calculate the probability of movement
# associated to each link based on the distance
edge_weights <- edge_attr(patchyMPG$mpg, "lcpPerimWeight")
all_edge_ids <- E(patchyMPG$mpg)

edge_info <- data.frame(
  edge_id = seq_len(gsize(patchyMPG$mpg)),
  get.edgelist(patchyMPG$mpg),
  edge_attr(patchyMPG$mpg, "lcpPerimWeight"),
  sapply(edge_attr(patchyMPG$mpg, "lcpPerimWeight"), inverse_exponential_probability)
)
colnames(edge_info) <- c("edge_id", "from", "to", "Weight (Sum of resistance)", "Probability of movement")
print(edge_info)

# We read each row to create the dictionnary
dictEdgeCaracteristics = list()
for (i in 1:nrow(edge_info)) {
  # We get the data for the edge from the dataframe
  origin = edge_info[[2]][i]
  end = edge_info[[3]][i]
  edgeId = edge_info[[1]][i]
  edgeProb = edge_info[[5]][i]
  # We create the nested dictionnary
  if (!(origin %in% names(dictEdgeCaracteristics))) {
    dictEdgeCaracteristics[[origin]] = list()
  }
  if (!(end %in% names(dictEdgeCaracteristics))) {
    dictEdgeCaracteristics[[end]] = list()
  }
  # We create the two links (as our links are bi-directional)
  dictEdgeCaracteristics[[origin]][[end]] = list()
  dictEdgeCaracteristics[[origin]][[end]][["ID"]] = edgeId
  dictEdgeCaracteristics[[origin]][[end]][["Prob"]] = edgeProb
  dictEdgeCaracteristics[[end]][[origin]] = list()
  dictEdgeCaracteristics[[end]][[origin]][["ID"]] = edgeId
  dictEdgeCaracteristics[[end]][[origin]][["Prob"]] = edgeProb
}

### Define area of each patch (accessible through dictionnary) ###
# Same thing as with edges/links before.
coreAreas <- vertex_attr(patchyMPG$mpg, "coreArea")
patchArea <- vertex_attr(patchyMPG$mpg, "patchArea")
all_vertex_ids <- V(patchyMPG$mpg)

vertex_infos <- data.frame(
  vertex_id = all_vertex_ids,
  patchArea,
  coreAreas
)
print(vertex_infos)

dictVertexCharacteristics = list()
for (i in 1:nrow(vertex_infos)) {
  # We get the data for the edge from the dataframe
  vertexID = vertex_infos[[1]][i]
  patchArea = vertex_infos[[2]][i]
  coreArea = vertex_infos[[3]][i]
  
  dictVertexCharacteristics[[vertexID]] = list()
  dictVertexCharacteristics[[vertexID]][["patchArea"]] = patchArea
  dictVertexCharacteristics[[vertexID]][["coreArea"]] = coreArea
}


sumPC = 0
pb <- txtProgressBar(min = 0, max = length(V(patchyMPG$mpg)$patchId), style = 3)
i = 0

# We create the tuples for the pair of patches
# Perplexity is doing this, as I'm at a loss at doing these things in a reasonable
# time in R.
# all_pairs <- expand.grid(V(patchyMPG$mpg)$patchId, V(patchyMPG$mpg)$patchId)
# result <- apply(all_pairs, 2, as.integer)
all_pairs <- combn(V(patchyMPG$mpg)$patchId, 2)
all_pairs_list <- lapply(seq_len(ncol(all_pairs)), function(i) as.integer(c(all_pairs[1, i], all_pairs[2, i])))

computeTermOfPCForPairOfPatch <- function(graphObject, pairOfPatches){
  areaOfPatch1 = V(graphObject)$patchArea[which(V(graphObject)$name == pairOfPatches[1])]
  areaOfPatch2 = V(graphObject)$patchArea[which(V(graphObject)$name == pairOfPatches[2])]
  probPatch1To2 = maxProbabilityConnectivityBetweenNodes(graphObject, pairOfPatches[1], pairOfPatches[2])
  
  return(areaOfPatch1*areaOfPatch2*probPatch1To2)
}

results <- pbsapply(all_pairs_list, computeTermOfPCForPairOfPatch, graphObject = patchyMPG$mpg)

sum(results)*2

for (patchID in V(patchyMPG$mpg)$patchId){
  for (otherPatchID in V(patchyMPG$mpg)$patchId){
    
  }
  i = i + 1
  setTxtProgressBar(pb, i)
}
# Then we divide by the landscape area squared
PCWithAllPatches = sumPC/(landscapeAreaM2**2)


# Create a list of tuples
tuples <- list(c(2, 3), c(4, 5), c(6, 7))

# Use sapply() to calculate the product of each tuple and sum them up
sum <- 0
result <- sapply(tuples, function(x) {
  product <- x[1] * x[2]
  sum <<- sum + product
  return(product)
})

# Print the result
print(result)
print(sum)




### Calculating the PC ###

# Now, we can finally compute the PC !

# For each pair of patch, we sum to get the top of the equation
sumPC = 0
pb <- txtProgressBar(min = 0, max = length(V(patchyMPG$mpg)$patchId), style = 3)
i = 0
for (patchID in V(patchyMPG$mpg)$patchId){
  for (otherPatchID in V(patchyMPG$mpg)$patchId){
    areaOfPatch1 = V(patchyMPG$mpg)$patchArea[which(V(patchyMPG@mpg)$name == patchID)]
    areaOfPatch2 = V(patchyMPG$mpg)$patchArea[which(V(patchyMPG@mpg)$name == otherPatchID)]
    probPatch1To2 = dictShortestPathAllPatches[[patchID]][[otherPatchID]][["maxProb"]]
    
    sumPC = sumPC + (areaOfPatch1*areaOfPatch2*probPatch1To2)
  }
  i = i + 1
  setTxtProgressBar(pb, i)
}
# Then we divide by the landscape area squared
PCWithAllPatches = sumPC/(landscapeAreaM2**2)



# Assuming you have already created your MPG object called 'mpg'

# Function to calculate PC of whole graph
calculate_pc <- function(graph) {
  n <- vcount(graph)
  
  # Get the adjacency matrix with edge weights
  adj_matrix <- as_adjacency_matrix(graph, attr="weight", sparse=FALSE)
  
  # Convert weights to probabilities if they aren't already
  adj_matrix <- 1 - adj_matrix  # Assuming weights represent distances
  
  # Calculate shortest paths
  paths <- shortest.paths(graph, weights=E(graph)$weight)
  
  # Calculate PC
  pc_sum <- 0
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      if (is.finite(paths[i,j])) {
        path <- get.shortest.paths(graph, i, j, weights=E(graph)$weight)$vpath[[1]]
        prob <- 1
        for (k in 1:(length(path)-1)) {
          prob <- prob * adj_matrix[path[k], path[k+1]]
        }
        pc_sum <- pc_sum + prob
      }
    }
  }
  
  pc <- (2 * pc_sum) / (n * (n - 1))
  return(pc)
}


# Initial PC
initial_PC <- calculate_PC(mpg@mpg)

# Calculate dPC for each node
dPC <- sapply(V(mpg@mpg), function(node) {
  # Remove node
  temp_graph <- delete_vertices(mpg@mpg, node)
  
  # Recalculate PC
  new_PC <- calculate_PC(temp_graph)
  
  # Calculate percent change
  dPC <- (initial_PC - new_PC) / initial_PC * 100
  
  return(dPC)
})

# Add dPC as a node attribute
V(mpg@mpg)$dPC <- dPC