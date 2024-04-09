## -------------------------------------------------------------------
## Script name: calc_conn_metrics
## Purpose of script: Calculate connectivity metrics in different ways.
## Author: Lei Song
## Date Created: 2023-12-14
## Email: lsong@ucsb.edu

## Usage
## Rscript path/to/calc_conn_metrics.R --src_dir data/raw/public 
## --dst_path path/for/result.csv
## -------------------------------------------------------------------

# Load libraries, easy to switch to use groundhog
pkgs <- c("sf", "dplyr", "terra", "stringr", "pbapply", "optparse")
sapply(pkgs, require, character.only = TRUE)
sf_use_s2(FALSE) # deal with buffering odd

# Command line inputs
option_list <- list(
    make_option(c("-t", "--taxon"), 
                action = "store", default = "mammal", type = 'character',
                help = "The taxon group in [bird, mammal] to process. [default %default]."),
    make_option(c("-s", "--src_dir"), 
                action = "store", default = "data/raw/public", type = 'character',
                help = "The source directory for reading data [default %default]."),
    make_option(c("-d", "--dst_dir"), 
                action = "store", default = 'data/derived/public', type = 'character',
                help = paste0("The path to save the csv [default %default].")))
opt <- parse_args(OptionParser(option_list = option_list))

# Directories and paths
src_dir <- opt$src_dir
dst_dir <- opt$dst_dir
taxon <- opt$taxon

# Load PAs and groups
fnames <- file.path(dst_dir, c("clean_pas.geojson", "pa_groups.shp"))
if (!all(file.exists(fnames))){
    stop("No cleaned pas and groups found, run clean_pas.R to get them.")
} else {
    clean_pas <- st_read(file.path(dst_dir, "clean_pas.geojson"))
    if (taxon == "mammal"){
        pa_groups <- st_read(file.path(dst_dir, "pa_groups.shp"))
    } else pa_groups <- NULL
}

# Read samples and raster template
## According to the pairs of lat/lon and east/north, they used
## UTM Zone 46 (EPSG:32646) for projection, so here we will use the same one.
## (Not mean it is the most right one to use)
fname <- list.files(file.path(src_dir, "training"), full.names = TRUE, 
                     pattern = taxon)
pts <- read.csv(fname) %>% 
    select(all_of(names(.)[
        str_detect(names(.), "station|country|lat|long")])) %>% 
    st_as_sf(coords = c(3, 2), crs = 4326) %>% 
    st_transform(st_crs(clean_pas))

template <- rast(
    file.path(src_dir, "GEDIv002_20190417to20220413_cover_krig.tif")) %>% 
    extend(c(100, 100)) %>%  # add a buffer
    project(crs(clean_pas))
values(template) <- 1:ncell(template)

# Now it is near perfect to calculate the flux and area weighted flux
## Get cell ids, NOTE that one cell could have multiple points
## Also remove points outside of study area
pts_clean <- terra::extract(template, pts, cells = TRUE, bind = TRUE) %>% 
    st_as_sf() %>% 
    filter(!is.na(cell)) %>% 
    select(station, country)

# Define the function to calculate the fluxes, areas, and area weighted fluxes
# Points inside of PA could (both are ecologically meaningful):
# - the distance of this points to other surrounding PAs.
# - or the boundary distance between the PA having this point and other PAs.
pas_conn <- function(pas, pts, template, 
                     pa_groups = NULL,
                     med_dist = 10,
                     # 5 times, the flux is around 0.0#
                     buffer_size = 5 * med_dist * 1000,
                     dst_path = NULL){
    # pas has two layers: first is the index of PAs, the second is the area
    # pts is a sf of points to process, at least have column station
    # template: grid template
    # pa_groups: PAs groups or NULL. Yes for mammal and NULL for birds.
    # med_dist is the median dispersal distance in KM.
    # buffer_size is the buffer size in RIGHT unit, e.g. degree or meter
    # dst_path: the path to save result
    
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
    if (!is.null(pa_groups)){
        pts <- pts %>% st_join(pa_groups)
    }
    
    # wenxin: I think this is max_dist here because 10 km is already the median dispersal distance
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

# Calculate and save out the result
conn_metrics <- do.call(
    rbind, lapply(seq(10, 150, 10), function(med_dist){
        message(sprintf("Calculate flux for median dispersal distance: %s", med_dist))
        dst_path <- file.path(dst_dir, sprintf("conn_flux_%s_%s.csv", taxon, med_dist))
        pas_conn(pas = clean_pas, pts = pts_clean, template = template, 
                 pa_groups = pa_groups, med_dist = med_dist, dst_path = dst_path)
    }))
dst_path <- file.path(dst_dir, sprintf("conn_flux_%s.csv", taxon))
write.csv(conn_metrics, dst_path, row.names = FALSE)
