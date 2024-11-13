## -----------------------------------------------------------------------------
## prep_dm
## -----------------------------------------------------------------------------
#
# Author:       Lei Song, Peter Kedron
# Date Created: 2024-08-13
# Last Update:  2024-08-25
# Email:        lsong@ucsb.edu
#
# Import from package: sf, dplyr, terra
#
## Purpose of script: ---------------------------------------------------------- 
# Prepare data model for the calculation of the connectivity metrics flux and 
# area weighted flux for a set of points based on four assumptions.
#
## Connectivity assumptions:
##  Vector based method 1 (*_ptg in outputs): Treat each species point 
##  (pts in inputs) as point and surrounding PAs as polygons. The awf are 
##  calculated based on the distance between this point and polygons.
##
##  If there is no surrounding PA polygons, all values are set to 0. For tiny 
##  polygons: The function uses `touches == TRUE` for rasterization to make sure 
##  using all polygons.

# Inputs: ---------------------------------------------------------------------
# pas (sf): The sf of protected areas. 
#
# pts (sf): The sf of points to process, at least have column station template 
#           (SpatRaster): The grid template for the calculation.
#
# pa_groups (sf): The PAs groups or NULL. Yes for mammal and NULL for birds.
#                 Check clean_pa.R for more details.
#
# med_dist (integer): The median dispersal distance in km.
#
# buffer_size (integer): The buffer size in the same unit as pas (e.g. degree or
#                         meter). Default is 5 * med_dist * 1000.
#
# dst_path (character): The path to save out the result or NULL to directly 
#                       return the result.
#
# Outputs: --------------------------------------------------------------------
# A data.frame of the metrics or save out metrics as a csv file.
#
# metrics:
##  station: station id
##  awf_ptg: area weighted flux (awf) for point to polygons
##  med_dist: The median dispersal distance in km. It is the same as input. 

# Usage example: ---------------------------------------------------------------
# See function calc_conn for a complete example.
#
## -----------------------------------------------------------------------------

# Define the function to calculate the fluxes and area weighted fluxes
prep_dm <- function(pas, 
                    pts, 
                    pa_groups = NULL,
                    med_dist = 10,
                    # 5 times, the flux is around 0.0#
                    buffer_size = 5 * med_dist * 1000,
                    dst_path = NULL){
  # If use PA groups (mammal), spatial join points with PA groups.
  if (!is.null(pa_groups)){
    pts <- pts %>% st_join(pa_groups)
  }
  
  # Calculate k in the flux/awf equation
  k <- -log(0.5) / (med_dist)
  
  metrics <- do.call(rbind, pbmclapply(unique(pts$station), function(sta_id){
    # Get involved polygons
    pts_this <- pts %>% filter(station == sta_id)
    pas_included <- suppressMessages(st_intersection(st_buffer(pts_this, buffer_size), pas))
    
    # Remove pas not at the same group
    if (!is.null(pa_groups)){
      pas_included <- pas_included %>% filter(group %in% pts_this$group)
    }
    
    # Calculate the metrics
    if (nrow(pas_included) > 0){
        plys <- rbind(pts_this %>% select(),
                      pas_included %>% select())
        dists <- st_distance(plys) %>% units::set_units("km") %>% 
            units::drop_units()
        
        areas <- as.numeric(c(1, pas_included$REP_AREA))
        
        for (i in 1:nrow(dists)){
            for (j in 1:ncol(dists)){
                dists[i, j] <- exp(-k * dists[i, j]) * areas[i] * areas[j]
            }
        }
        diag(dists) <- 0
        
      # gather results
      data.frame(station = sta_id, awf_ptg = sum(dists))
    } else data.frame(station = sta_id, awf_ptg = 0)
  }, mc.cores = 6))
  
  # add common parameters
  metrics <- metrics %>% mutate(med_dist = med_dist)
  
  # save out
  if (!is.null(dst_path)){
    write.csv(metrics, dst_path, row.names = FALSE)
  }
  
  # return
  metrics
}
