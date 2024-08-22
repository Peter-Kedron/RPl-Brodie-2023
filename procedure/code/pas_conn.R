## -------------------------------------------------------------------
## Script name: pas_conn
## Purpose of script: The function to and calculate connectivity 
## metrics of flux and area weighted flux.It calculates the 
## fluxes, areas, and area weighted fluxes for points inside of
## PA could (both are ecologically meaningful):
## - the distance of this points to other surrounding PAs.
## - or the boundary distance between the PA having this point and 
## other PAs.
## Author: Lei Song
## Date Created: 2024-08-13
## Email: lsong@ucsb.edu

## Import from package:
## sf, dplyr, terra

## Inputs:
## pas (vector): The protected area file. pas has two layers: first is 
## the index of PAs, the second is the area
## pts (vector): An sf object of points to process, should at least have 
## a column for station ids.
## template (raster): A raster template.
# pa_groups (boolean): PAs groups or NULL. Yes for mammal and NULL for birds.
# med_dist (numeric): The median dispersal distance in KM.
# buffer_size (numeric): The buffer size in a CORRECT unit, e.g. degree or meter.
# dst_path (character): the path to save result

## Outputs:
# flux_ptg: flux for point to polygons
# awf_ptg: awf for point to polygons
# flux_gtg: flux for polygon (if has point) to polygons
# awf_gtg: awf for for polygon (if has point) to polygons
# flux_rst_ptp: flux for pixel to pixels
# awf_rst_ptp: awf for pixel to pixels
# flux_rst_ptp2: flux for pixel to other PA pixels
# awf_rst_ptp2: awf for pixel to other PA pixels

## Usage example:
# library(optparse)

# 
# # Directories and paths
# src_dir <- opt$src_dir
# dst_dir <- opt$dst_dir
# taxon <- opt$taxon
## -------------------------------------------------------------------
pas_conn <- function(pas, pts, template, 
                     pa_groups = NULL,
                     med_dist = 10,
                     # 5 times, the flux is around 0.0#
                     buffer_size = 5 * med_dist * 1000,
                     dst_path = NULL){
    
    # Q3: how to deal with very small polygons for raster based method?
    # now the code uses `touches == TRUE` for rasterization to use all polygons
    
    # Calculate
    if (!is.null(pa_groups)){
        pts <- pts %>% st_join(pa_groups)
    }
    
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