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
##  1. Vector based method 1 (*_ptg in outputs): Treat each species point 
##     (pts in inputs) as point and surrounding PAs as polygons. The flux and  
##     awf are calculated based on the distance between this point and polygons.
##  2. Vector based method 2 (*_gtg in outputs): Extract the PA polygon 
##     with each species point (pts in inputs) as the representation 
##     of this species point. The flux and awf are calculated based on the  
##     distance between this polygon and surrounding polygons.
##  3. Raster based method 1 (*_rst_ptp in outputs): Treat each species point
##     (pts in inputs) as a pixel in template grid and the surrounding PA 
##     polygons as a set of pixels. The flux and awf are calculated based on the  
##     distance between the center of this species pixel and surrounding pixels.
##  4. Raster based method 2 (*_rst_ptp2 in outputs): The logic of converting 
##     vectors to rasters is the same as raster based method 1. The only 
##     difference is to exclude PA polygon that the species point is within from
##     the surrounding PA polygons.
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
##  flux_ptg: flux for point to polygons
##  awf_ptg: area weighted flux (awf) for point to polygons
##  flux_gtg: flux for polygon (if has point) to polygons
##  awf_gtg: awf for for polygon (if has point) to polygons
##  flux_rst_ptp: flux for pixel to pixels
##  awf_rst_ptp: awf for pixel to pixels
##  flux_rst_ptp2: flux for pixel to other PA pixels
##  awf_rst_ptp2: awf for pixel to other PA pixels
##  med_dist: The median dispersal distance in km. It is the same as input. 



# Usage example: ---------------------------------------------------------------
# See function calc_conn for a complete example.
#
## -----------------------------------------------------------------------------

# Define the function to calculate the fluxes and area weighted fluxes
# Points inside of PA could (both are ecologically meaningful):
# - the distance of this points to other surrounding PAs.
# - or the boundary distance between the PA having this point and other PAs.
prep_dm <- function(pas, 
                    pts, 
                    template, 
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
  
  metrics <- do.call(rbind, pblapply(unique(pts$station), function(sta_id){
    # Get involved polygons
    pts_this <- pts %>% filter(station == sta_id)
    pas_included <- suppressMessages(st_intersection(st_buffer(pts_this, buffer_size), pas))
    
    # Remove pas not at the same group
    if (!is.null(pa_groups)){
      pas_included <- pas_included %>% filter(group %in% pts_this$group)
    }
    
    # Calculate the metrics
    if (nrow(pas_included) > 0){
      # Vector method 1: points to polygons, flux and whole area weighted flux
      dists <- as.vector(st_distance(pts_this, pas_included) %>% 
                           units::set_units("km"))
      flux_ptg <- sum(exp(- k * dists))
      awf_ptg <- sum(exp(- k * dists) * as.numeric(pas_included$REP_AREA))
      
      ## Check if this point is inside of a PA
      id_pfp <- unlist(suppressMessages(st_intersects(pts_this, pas_included)))
      
      if(length(id_pfp) > 0){
        # PA point
        # Raster method 1: Exclude this point from the grid
        ## Remove the point from the raster
        pas_included_grid <- rasterize(
          pas_included, template, touches = TRUE) %>% trim()
        cell_pts <- terra::extract(
          pas_included_grid, pts_this, cells = TRUE) %>% pull(cell)
        pas_included_grid[cell_pts] <- NA
        
        ## Calculate the values
        pas_included_grid <- pas_included_grid %>% 
          as.points() %>% st_as_sf()
        dists <- as.vector(st_distance(pts_this, pas_included_grid) %>% 
                             units::set_units("km"))
        flux_rst_ptp <- sum(exp(- k * dists))
        # Use the number of pixels to indicate area thing
        awf_rst_ptp <- sum(exp(- k * dists)) * nrow(pas_included_grid)
        
        if (nrow(pas_included) != length(id_pfp)){
          # Vector method 2:
          # if the point is inside of PA, then use the boundary distance from
          # this polygon to others.
          dists <- as.vector(st_distance(pas_included[id_pfp, ], 
                                         pas_included[-id_pfp, ]) %>% 
                               units::set_units("km"))
          flux_gtg <- sum(exp(- k * dists))
          awf_gtg <- sum(exp(- k * dists) * 
                           as.numeric(pas_included$REP_AREA[-id_pfp]))
          
          # Raster method 2: Exclude the whole polygon from the grid
          ## Remove the polygon
          pas_included_grid <- rasterize(
            pas_included[-id_pfp, ], template, touches = TRUE) %>% 
            trim() %>% as.points() %>% st_as_sf()
          
          ## Calculate the values
          dists <- as.vector(st_distance(pts_this, pas_included_grid) %>% 
                               units::set_units("km"))
          flux_rst_ptp2 <- sum(exp(- k * dists))
          awf_rst_ptp2 <- sum(exp(- k * dists)) * nrow(pas_included_grid) 
        } else{
          # The same as there is no PA involved
          flux_gtg <- awf_gtg <- flux_rst_ptp2 <- awf_rst_ptp2 <- 0
        }
        
      } else {
        # Non-PA point, the same as point to polygons
        flux_gtg <- flux_ptg
        awf_gtg <- awf_ptg
        
        # Raster method 1 and 2, the same
        pas_included_grid <- rasterize(
          pas_included, template, touches = TRUE) %>% trim() %>% 
          as.points() %>% st_as_sf()
        dists <- as.vector(st_distance(pts_this, pas_included_grid) %>% 
                             units::set_units("km"))
        flux_rst_ptp <- flux_rst_ptp2 <- sum(exp(- k * dists))
        awf_rst_ptp <- awf_rst_ptp2 <- sum(exp(- k * dists)) * 
          nrow(pas_included_grid)
      }
      
      # gather results
      data.frame(station = sta_id, flux_ptg = flux_ptg, awf_ptg = awf_ptg,
                 flux_gtg = flux_gtg, awf_gtg = awf_gtg,
                 flux_rst_ptp = flux_rst_ptp, awf_rst_ptp = awf_rst_ptp,
                 flux_rst_ptp2 = flux_rst_ptp2, awf_rst_ptp2 = awf_rst_ptp2)
    } else {
      data.frame(station = sta_id, flux_ptg = 0, awf_ptg = 0,
                 flux_gtg = 0, awf_gtg = 0,
                 flux_rst_ptp = 0, awf_rst_ptp = 0,
                 flux_rst_ptp2 = 0, awf_rst_ptp2 = 0)
    }
  }))
  
  # add common parameters
  metrics <- metrics %>% mutate(med_dist = med_dist)
  
  # save out
  if (!is.null(dst_path)){
    write.csv(metrics, dst_path, row.names = FALSE)
  }
  
  # return
  metrics
}
