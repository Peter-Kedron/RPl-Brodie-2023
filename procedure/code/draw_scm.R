## -------------------------------------------------------------------
## Script name: draw_scm
## Purpose of script: Do the structural Causal Modeling.

## Usage
## Better to use in an interactive environment, but you could run it by
## Rscript path/to/draw_scm.R
## -------------------------------------------------------------------


# Load libraries
library(groundhog)
pkgs <- c("dagitty", "ggdag")
groundhog.library(pkgs, "2023-12-05")

# ------------------------------------------------------------------------------
# Structural Causal Modelling  
# ------------------------------------------------------------------------------

# Create the directed acyclic graph (DAG) of potential causal pathways with the 
# addition of connectivity as a moderator of the PA effect. Assign PA as 
# the exposure and diversity as the outcome.

dagBrodie <- dagitty("dag {
  PA -> Diversity
  Connectivity -> Diversity
  ForestStructure -> PA
  SiteAccessibility -> PA
  Bioclimate -> ForestStructure -> PA -> Diversity
  Bioclimate -> UnderstoryDensity -> Diversity
  Bioclimate -> Diversity
  Elevation -> Bioclimate
  Elevation -> ForestStructure
  Elevation -> UnderstoryDensity
  Elevation -> SiteAccessibility
  Topography -> UnderstoryDensity
  Topography -> ForestStructure
  Topography -> SiteAccessibility
  Topography -> Diversity
  HDI -> Diversity
  PA [exposure]
  Diversity [outcome]
               }"
)

# Organize the data into a visual hierarchy
coordinates( dagBrodie ) <-  list(x = c(Diversity = 3, 
                                        UnderstoryDensity = 1, 
                                        ForestStructure = 2, 
                                        PA = 3, 
                                        SiteAccessibility = 4, 
                                        HDI = 5, 
                                        Bioclimate = 2, 
                                        Elevation = 3, 
                                        Topography = 4,
                                        Connectivity = 6),
                                  y = c(Diversity = 3, 
                                        UnderstoryDensity = 2, 
                                        ForestStructure = 2, 
                                        PA = 2, 
                                        SiteAccessibility = 2, 
                                        HDI = 2, 
                                        Bioclimate = 1, 
                                        Elevation = 1, 
                                        Topography = 1, 
                                        Connectivity = 2))

# Plot the DAG to confirm the visual structure and exposure
ggdag_status(dagBrodie) + theme_dag()

# Identify the set of adjustment variables needed to identify 
# the effect of PA on biodiversity and Test whether connectivity fulfills the 
# adjustment criterion
adjustmentSets(dagBrodie)
isAdjustmentSet(dagBrodie, c("Connectivity"))

# Plot the alternative adjustment sets
ggdag_adjustment_set(dagBrodie, 
                     node_size = 20, 
                     text_col = "black"
) + theme(legend.position = "bottom")
