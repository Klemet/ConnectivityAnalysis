#####################################################
#
# @author: Clement Hardy (clem.hardy@pm.me)
#
# This script is used to compute the Ecological Connectivity (EC) of
# a complete graph/network of protected areas in a landscape (complete means
# that all possible paths between patches are taken into account).
#
# I use Complete graphs and not Minimal planar graphs as they are better adapted
# to really compute all of the connectivity in a landscape.
# 
# WARNING : This script uses functions coming from the graph4lg package that use
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
# IN CASE OF other ERRORS WITH THE GRAPHAB FUNCTIONS : open a command prompt in windows,
# paste the command run by the graph4lg package into the command prompt, and check the error message.
#####################################################



#### LOADING PACKAGES ####
library(fs)
library(grainscape)
library(igraph)
library(raster)
library(knitr)
library(ggplot2)
library(sf)
library(Dict)
library(rlang)
library(pbapply)
library(graph4lg)
library(dplyr)
library(terra)

#### LOADING FUNCTIONS ####

source("./Networks_Analysis/functions_EC_MPG_Analysis.r")

#### SET THE DIRECTORIES AND FILE LOCATIONS ####

# WARNING : BE SURE TO USE A PATH WITHOUT ACCENTS; graphab doesn't like that.
setwd("./Networks_Analysis/GraphabProjects")
temporaryRasterPath <- "./Networks_Analysis/GraphabProjects/TemporaryRaster.tif"

#### LOADING DATA ####

# We set the folder with the inputs

# Resistance map (see CreatingResistanceMap_PitherEtAl2023.py for details)
resistanceMapPath = "./ResistanceCostMapCreation/CostMaps/MEDIUM (100x100m)/resistanceCostMapPitherEtAl2023.tif"
resistanceMap <- raster(resistanceMapPath)

# We read all vector files; these must contain the position of the different protected area of your different scenarios.
# The polygons in the different files can superpose.
# These must be situated inside the same extent as your resistance map, or they won't be taken into account later.
shp_currentProtectedArea <- st_read("./Networks_Analysis/InputFiles/currentProtectedAreas.shp")
shp_AdditionalProtectedArea1 <- st_read("./Networks_Analysis/InputFiles/AdditionalProtectedArea1.shp")
shp_AdditionalProtectedArea2 <- st_read("./Networks_Analysis/InputFiles/AdditionalProtectedArea2.shp")

# Landscape area : extent of the cost raster
landscapeAreaM2 = res(resistanceMap)[1]^2 * (ncell(resistanceMap))


#### COMPUTING EC ####
# These coming section will compute the Ecological Connectivity (EC)
# for the different configuration of polygons in the Landscape
# The polygons are transformed into patches for the network/graph, and
# the complete graph is then computed by graphab in every case.
# The Probability of connectivity and Ecological Connectivity are then computed
# based on the outputs of graphab.
#
# You can edit or manipulate your vector layers with R functions before passing them
# as a list to the computePC_and_EC_graphab.
#
# The details of the function computePC_and_EC_graphab are in another file (functions_EC_MPG_Analysis.R).

#### FIRST CASE : CURRENT LANDSCAPE ####
graphabProjectName = "CurrentLandscape"
listOfShapefileLayers = list(shp_currentProtectedArea)

results_CurrentLandscape = computePC_and_EC_graphab(listOfShapefileLayers,
                                                     resistanceMapPath,
                                                     temporaryRasterPath,
                                                     graphabProjectName)

#### ALTERNATIVES SCENARIOS ####

##########################################################
## Adding protected area 1
graphabProjectName = "Network_AdditionalProtectedArea1"
listOfShapefileLayers = list(shp_currentProtectedArea,
                             shp_AdditionalProtectedArea1)
results_AdditionalProtectedArea1 = computePC_and_EC_graphab(listOfShapefileLayers,
                                                    resistanceMapPath,
                                                    temporaryRasterPath,
                                                    graphabProjectName)

##########################################################
## Adding protected area 1
graphabProjectName = "Network_AdditionalProtectedArea2"
listOfShapefileLayers = list(shp_currentProtectedArea,
                             shp_AdditionalProtectedArea2)
results_AdditionalProtectedArea2 = computePC_and_EC_graphab(listOfShapefileLayers,
                                                    resistanceMapPath,
                                                    temporaryRasterPath,
                                                    graphabProjectName)

#### SUMMARIZING RESULTS ####

listOfResults = list(
  "Current landscape" = results_CurrentLandscape,
  "Adding the protected area 1" = results_AdditionalProtectedArea1,
  "Adding the protected area 2" = results_AdditionalProtectedArea2
  )
PC_values = c()
EC_values = c()
for (scenario in names(listOfResults))
{
  PC_values = c(PC_values, listOfResults[scenario][[1]]$PC)
  EC_values = c(EC_values, listOfResults[scenario][[1]]$EC)
}


#### CREATING FIGURE ####
# Code comes from https://r-graph-gallery.com/web-horizontal-barplot-with-labels-the-economist.html
# 
# Code will output a .PNG file called plotEC_Variation_Scenarios.png in the working directory.
# 
## Load packages
library(grid)
library(tidyverse)
library(shadowtext)
library(extrafont)

# Initialize fonts
font_import()
loadfonts(device = "win")

# Create Data
# names <- c(
#   "Hantavirus", "Tularemia", "Dengue", "Ebola", "E. coli", 
#   "Tuberculosis", "Salmonella", "Vaccinia", "Brucella"
# )
names <- names(listOfResults)

# Name is an ordered factor. We do this to ensure the bars are sorted.
data <- data.frame(
  PC = PC_values,
  EC = EC_values,
  name = factor(names, levels = names),
  y = seq(length(names)) * 0.9
)

# We make EC smaller (in hectares instead of m2)
data$EC = data$EC/10000

# We further make EC smaller (in 100 of hectares instead of ha)
data$EC = data$EC/100

# We sort the data frame and reset the index of rows
data <- data[order(data$EC, decreasing = TRUE), ]
rownames(data) <- NULL

# Define colors
# The colors
BLUE <- "#5e81ac"
RED <- "#bf616a"
BLACK <- "#2e3440"
GREY <- "#d8dee9"
GREEN <- "#a3be8c"

# We indicate the color we want for each bar in the data frame
data$color = rep(BLUE, length(data$EC))
data[data$name == "Current Landscape",]$color = GREEN


# Basic Barchart
plt <- ggplot(data) +
  geom_col(aes(EC, reorder(name, -EC)), fill = data$color, width = 0.6) 

plt

# Customize layout
plt <- plt + 
  scale_x_continuous(
    limits = c(0, round(max(data$EC*1.01), -3)),
    breaks = seq(0, round(max(data$EC*1.01), -3), by = round(max(data$EC)/5, -3)), 
    expand = c(0, 0), # The horizontal axis does not extend to either side
    position = "top"  # Labels are located on the top
  ) +
  # The vertical axis only extends upwards 
  scale_y_discrete(expand = expansion(add = c(0, 0.5))) +
  theme(
    # Set background color to white
    panel.background = element_rect(fill = "white"),
    # Set the color and the width of the grid lines for the horizontal axis
    panel.grid.major.x = element_line(color = BLUE, size = 0.3),
    # Remove tick marks by setting their length to 0
    axis.ticks.length = unit(0, "mm"),
    # Remove the title for both axes
    axis.title = element_blank(),
    # Only left line of the vertical axis is painted in black
    axis.line.y.left = element_line(color = "black"),
    # Remove labels from the vertical axis
    axis.text.y = element_blank(),
    # But customize labels for the horizontal axis
    axis.text.x = element_text(family = "Lato", size = 16)
  )

plt

# Add labels
plt <- plt + 
  geom_shadowtext(
    data = subset(data, EC < max(data$EC)*0.3749666),
    aes(EC, y = name, label = name),
    hjust = 0,
    nudge_x = max(data$EC)*0.01071333,
    colour = subset(data, EC < max(data$EC)*0.3749666)$color,
    bg.colour = "white",
    bg.r = 0.2,
    family = "Lato",
    size = 3
  ) + 
  geom_text(
    data = subset(data, EC >= max(data$EC)*0.3749666),
    aes(0, y = name, label = name),
    hjust = 0,
    nudge_x = max(data$EC)*0.01071333,
    colour = "white",
    family = "Lato",
    size = 3
  )

plt

# Add anotations and final tweaks
plt <- plt +
  labs(
    title = "Equivalent Connectivity (EC) in the landscape for different scenarios",
    subtitle = "In hundred of hectares of equivalent connectivity"
  ) + 
  theme(
    plot.title = element_text(
      family = "Lato", 
      face = "bold",
      size = 22
    ),
    plot.subtitle = element_text(
      family = "Lato",
      size = 20
    )
  )

# Make room for annotations
plt <- plt + 
  theme(
    plot.margin = margin(0.05, 0, 0.05, 0.01, "npc")
  )

# Print the ggplot2 plot
plt

# Add horizontal line on top
# It goes from x = 0 (left) to x = 1 (right) on the very top of the chart (y = 1)
# You can think of 'gp' and 'gpar' as 'graphical parameters'.
# There we indicate the line color and width
grid.lines(
  x = c(0, 1),
  y = 1,
  gp = gpar(col = RED, lwd = 4)
)

# Add rectangle on top-left
# lwd = 0 means the rectangle does not have an outer line
# 'just' gives the horizontal and vertical justification
grid.rect(
  x = 0,
  y = 1,
  width = 0.05,
  height = 0.025,
  just = c("left", "top"),
  gp = gpar(fill = RED, lwd = 0)
)

# We have two captions, so we use grid.text instead of 
# the caption provided by  ggplot2.
# grid.text(
#   "Sources: Laboratory-Acquired Infection Database; American Biological Safety Association", 
#   x = 0.005, 
#   y = 0.06, 
#   just = c("left", "bottom"),
#   gp = gpar(
#     col = GREY,
#     fontsize = 16,
#     fontfamily = "Lato"
#   )
# )
# grid.text(
#   "The Economist", 
#   x = 0.005, 
#   y = 0.005, 
#   just = c("left", "bottom"),
#   gp = gpar(
#     col = GREY,
#     fontsize = 16,
#     fontfamily = "Lato"
#   )
# )

# The extra mile

plt <- plt + 
  labs(title = NULL, subtitle = NULL) +
  theme(
    plot.margin = margin(0.15, 0, 0.05, 0.01, "npc")
  )

plt 


# We open the connection to save the plot
# Necessary to add these few elements

# Open file to store the plot
png(paste("plotEC_Variation_Scenarios.png", sep = ""),
    width = 10, height = 5, units = "in", res = 300)

# Print the plot
plt 

grid.text(
  "Equivalent Connectivity (EC) in the landscape for different scenarios", 
  0.006, 
  0.925,
  just = c("left", "bottom"),
  gp = gpar(
    fontsize = 15,
    fontface = "bold",
    fontfamily = "Lato"
  )
)

grid.text(
  "In hundred of hectares of Equivalent Connectivity", 
  0.006, 
  0.875,
  just = c("left", "bottom"),
  gp = gpar(
    fontsize = 13,
    fontfamily = "Lato"
  )
)

grid.lines(
  x = c(0, 1),
  y = 1,
  gp = gpar(col = RED, lwd = 4)
)

grid.rect(
  x = 0,
  y = 1,
  width = 0.05,
  height = 0.025,
  just = c("left", "top"),
  gp = gpar(fill = RED, col = NA, lwd = 0)
)

# Close connection to file
dev.off()