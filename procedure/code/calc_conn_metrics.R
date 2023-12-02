# Load libraries
library(sf)
library(dplyr)
library(terra)
library(stringr)
library(pbmcapply)
sf_use_s2(FALSE) # deal with buffering odd

# Directories
src_dir <- "data/raw/public"

# Read samples and raster template
fnames <- list.files(file.path(src_dir, "training"), full.names = TRUE)
pts <- do.call(rbind, lapply(fnames, function(fname){
    read.csv(fname) %>% 
        select(all_of(names(.)[
            str_detect(names(.), "station|country|lat|long")])) %>% 
        st_as_sf(coords = c(3, 2), crs = 4326)
}))

template <- rast(
    file.path(src_dir, "GEDIv002_20190417to20220413_cover_krig.tif")) %>% 
    extend(c(100, 100)) # add a buffer
values(template) <- 1:ncell(template)

# Get WDPA, too lazy so just download the gdb and put it in the right place
st_layers(
    file.path(src_dir, "protected_area", "WDPA_WDOECM_Nov2023_Public_AS",
              "WDPA_WDOECM_Nov2023_Public_AS.gdb"))
pas <- read_sf(
    file.path(src_dir, "protected_area", "WDPA_WDOECM_Nov2023_Public_AS",
              "WDPA_WDOECM_Nov2023_Public_AS.gdb"),
    layer = "WDPA_WDOECM_poly_Nov2023_AS")

# Q1: should we only consider terrestrial PAs?
pas <- pas %>% 
    filter(ISO3 %in% c("KHM", "CHN", "IDN", "LAO", "MYS", 
                       "SGP", "THA", "VDR", "SVR", "BRN")) %>% 
    select(WDPAID, NAME, REP_AREA, GIS_AREA, REP_M_AREA, GIS_M_AREA) %>% 
    rename(Geometry = SHAPE) %>% 
    slice(unique(unlist(st_intersects(st_as_sfc(st_bbox(template)), .)))) %>% 
    mutate(id = 1:nrow(.))

# Split MULTIPOLYGON to POLYGONS because each segment should be treated independently
# Hmm... then REP_AREA is no longer usable, so calculate the GIS area
# Q2: I think I did right thing here, didn't I?
pas <- st_cast(pas, "POLYGON") %>% st_make_valid() %>% 
    mutate(REP_AREA = st_area(.) %>% units::set_units("km2"))
# Warning: some PAs may overlap (sometimes highly).

# Okay, ready to calculate the flux and area weighted flux
## Get cell ids, NOTE that one cell could have multiple points
## Also remove points outside of study area
pts_clean <- extract(template, pts, cells = TRUE, bind = TRUE) %>% 
    st_as_sf() %>% 
    filter(!is.na(cell)) %>% 
    select(station, country)

# Define the function to calculate the fluxes, areas, and area weighted fluxes
# Points inside of PA could (both are ecologically meaningful):
# - the distance of this points to other surrounding PAs.
# - or the boundary distance between the PA having this point and other PAs.
pas_conn <- function(pas, pts, template, buffer_size = 0.083, max_dist = 10){
    # pas has two layers: first is the index of PAs, the second is the area
    # pts is a sf of points to process, at least have column station
    # template: grid template
    # buffer_size is the buffer size in right unit, e.g. degree or meter
    # max_dist is the maximum dispersal distance in km
    
    # result
    # flux_ptg: flux for point to polygons
    # awf_ptg: awf for point to polygons
    # flux_gtg: flux for polygon (if has point) to polygons
    # awf_gtg: awf for for polygon (if has point) to polygons
    # flux_rst_ptp: flux for pixel to pixels
    # awf_rst_ptp: awf for pixel to pixels
    # flux_rst_ptp2: flux for pixel to other PA pixels
    # awf_rst_ptp2: awf for pixel to other PA pixels
    
    # Q3: how to deal with very small polygons for raster based method?
    # now the code uses `touches == TRUE` for rasterization to use all polygons
    
    # Calculate
    # wenxin: I think this is max_dist here because 10 km is already the median dispersal distance
    k <- -log(0.5) / (max_dist / 2)
    
    metrics <- pbmclapply(unique(pts$station), function(sta_id){
        message(sta_id)
        # Get involved polygons
        pts_this <- pts %>% filter(station == sta_id)
        pas_included <- st_intersection(st_buffer(pts_this, buffer_size), pas)
        
        # Calculate the metrics
        if (nrow(pas_included) > 0){
            # Vector method 1: points to polygons, flux and whole area weighted flux
            dists <- as.vector(st_distance(pts_this, pas_included) %>% 
                                   units::set_units("km"))
            flux_ptg <- sum(exp(- k * dists))
            awf_ptg <- sum(exp(- k * dists) * as.numeric(pas_included$REP_AREA))
            
            ## Check if this point is inside of a PA
            id_pfp <- unlist(st_intersects(pts_this, pas_included))
            
            if(length(id_pfp) > 0){
                # PA point
                # Raster method 1: Exclude this point from the grid
                ## Remove the point from the raster
                pas_included_grid <- rasterize(
                    pas_included, template, touches = TRUE) %>% trim()
                cell_pts <- extract(
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
    }, mc.cores = detectCores()) %>% bind_rows()
    
    # return
    metrics %>% mutate(buffer_size = buffer_size, max_dist = max_dist)
}

conn_metrics <- pas_conn(pas, pts, template, buffer_size = 0.083, max_dist = 10)
