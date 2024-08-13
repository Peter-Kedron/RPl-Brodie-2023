## -------------------------------------------------------------------
## Script name: calc_conn
## Purpose of script: The main function to load inputs and calculate 
## connectivity metrics of flux and area weighted flux.
## Author: Lei Song
## Date Created: 2023-12-14
## Email: lsong@ucsb.edu

## Import from package:
## sf, dplyr, terra

## Inputs:
## taxon (character): The taxon name. Either mammal or bird.
## med_dists (vector): A vector of median dispersal distance in km.
## See details in function prep_dm.
## src_dir (character): The directory of data source.
## dst_dir (character): The directory to save files to.

## Outputs:
## A combined data.frame of output from function prep_dm calculated with
## each med_dists. See details in function prep_dm.

## Usage example:
# library(optparse)
# option_list <- list(
#   make_option(c("-t", "--taxon"), 
#               action = "store", default = "mammal", type = 'character',
#               help = "The taxon group in [bird, mammal] to process. [default %default]."),
#   make_option(c("-s", "--src_dir"), 
#               action = "store", default = "data/raw/public", type = 'character',
#               help = "The source directory for reading data [default %default]."),
#   make_option(c("-d", "--dst_dir"), 
#               action = "store", default = 'data/derived/public', type = 'character',
#               help = paste0("The path to save the csv [default %default].")))
# opt <- parse_args(OptionParser(option_list = option_list))
# 
# # Directories and paths
# src_dir <- opt$src_dir
# dst_dir <- opt$dst_dir
# taxon <- opt$taxon
# calc_conn(taxon = taxon, src_dir = src_dir, dst_dir = dst_dir)
## -------------------------------------------------------------------

# Start the function
calc_conn <- function(taxon,
                      med_dists = seq(10, 150, 10),
                      src_dir,
                      dst_dir){
  # Check inputs
  if (!taxon %in% c("bird", "mammal")){
    stop("Taxon must be either bird or mammal.")
  }
  
  # Load prep_dm function
  source(here("procedure/code/prep_dm.R"))
  
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
  fname <- file.path(
    src_dir, "training", sprintf("%s_data_corrected_240122.csv", taxon))
  pts <- read.csv(fname) %>% 
    select(all_of(names(.)[
      str_detect(names(.), "station|utm")])) %>% 
    st_as_sf(coords = c(2, 3), crs = 32646)
  
  # Assume the users already downloaded these layers. 
  fname <- file.path(src_dir, "GEDIv002_20190417to20220413_cover_krig.tif")
  if (file.exists()){
    template <- rast(fname) %>% 
      extend(c(100, 100)) # add a buffer
    values(template) <- 1:ncell(template)
  } else{
    stop("Please download GEDI layers first.")
  }
  
  # Now it is near perfect to calculate the flux and area weighted flux
  ## Get cell ids, NOTE that one cell could have multiple points
  ## Also remove points outside of study area
  pts_clean <- terra::extract(template, pts, cells = TRUE, bind = TRUE) %>% 
    st_as_sf() %>% 
    filter(!is.na(cell)) %>% 
    select(station)
  
  # Calculate and save out the result
  conn_metrics <- do.call(
    rbind, lapply(med_dists, function(med_dist){
      message(sprintf("Calculate flux for median dispersal distance: %s", med_dist))
      # Save out files in case the workflow is broken.
      dst_path <- file.path(dst_dir, sprintf("conn_flux_%s_%s.csv", taxon, med_dist))
      pas_conn(pas = clean_pas, pts = pts_clean, template = template, 
               pa_groups = pa_groups, med_dist = med_dist, dst_path = dst_path)
    }))
  
  # Save out
  dst_path <- file.path(dst_dir, sprintf("conn_flux_%s.csv", taxon))
  write.csv(conn_metrics, dst_path, row.names = FALSE)
}
