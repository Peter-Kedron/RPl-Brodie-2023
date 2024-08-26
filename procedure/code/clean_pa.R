## -----------------------------------------------------------------------------
## clean_pa
## -----------------------------------------------------------------------------
#
# Purpose of script: -----------------------------------------------------------
# Clean the Protected areas.
#
# Author:       Lei Song, Peter Kedron
# Date Created: 2024-08-13
# Last Update:  2024-08-25
# Email:        lsong@ucsb.edu

# Import from package: sf, dplyr, terra, wdpar, countrycode

# Inputs: ---------------------------------------------------------------------
# geometry_precision (integer): The precision to deal with geometry. See details 
#                               in function st_set_precision in sf package. 
#                               Default is 100,000.
#
# dst_crs (character): The CRS to use. The default is "EPSG:32646".
#
# bbox_buffer (integer): The size set to buffer the bounding box to clip PAs.
#                        Unit should be the same as dst_crs. Default is 800,000 
#                        in meters.
#
# src_dir (character): The directory of data source.
#
# dst_dir (character): The directory to save files to.
#
# Outputs: --------------------------------------------------------------------
## No return values. Two files are created and saved:
## 1. clean_pas.geojson: cleaned PAs
## 2. pa_groups.shp: PA groups for mammal connectivity

# Usage example: ---------------------------------------------------------------
# 
# library(optparse)
# option_list <- list(
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
# clean_pa(src_dir = src_dir, dst_dir = dst_dir)
#
## -----------------------------------------------------------------------------

clean_pa <- function(geometry_precision = 100000, 
                     dst_crs = "EPSG:32646",
                     bbox_buffer = 800000,
                     src_dir, 
                     dst_dir){
  ################# Setting and inputs #################
  # 1. Switch off s2 to avoid buffering issues
  org_setting <- sf_use_s2()
  sf_use_s2(FALSE)
  
  # 2. Read raster template
  fname <- file.path(src_dir, "GEDIv002_20190417to20220413_cover_krig.tif")
  if (file.exists()){
    template <- rast(fname) %>% 
      extend(c(100, 100)) # add a buffer
    values(template) <- 1:ncell(template)
  } else{
    stop("Please download GEDI layers first.")
  }
  
  # 3. Get the associated country names
  fnames <- list.files(file.path(src_dir, "training"), full.names = TRUE)
  cnts <- sapply(fnames, function(fname){
    read.csv(fname) %>% pull(country)
  }) %>% unlist() %>% unique()
  cnts <- countrycode(cnts, origin = 'country.name', destination = 'iso3c')
  rm(fnames)
  
  ################# Query and clean Protected areas (PAs) #################
  ## Within this part, most of the issues related to PAs are solved 
  ## According to the pairs of lat/lon and east/north, they used
  ## UTM Zone 46 (EPSG:32646) for projection.
  ## 1. Clean and project the PAs or relevant layers to use precise distance.
  ## 2. Trim the PAs to terrestrial only.
  ## 3. Separate or union polygons and re-index them.
  
  # Query PAs
  raw_pas <- lapply(
    cnts, wdpa_fetch, wait = TRUE,
    download_dir = rappdirs::user_data_dir("wdpar")) %>%
    bind_rows()
  
  # Clean and project the PAs, see details in function wdpa_clean
  ## Note: PAs include both polygons and point with buffer of designed size.
  clean_pas <- wdpa_clean(raw_pas, crs = dst_crs,
                        geometry_precision = geometry_precision)
  
  # Trim the PAs to terrestrial only
  ## Remove marine PAs
  clean_pas <- clean_pas %>% filter(MARINE != "marine")
  ## Note: now the No is 1638
  
  ## Crop the PAs that cover both land and marine to only land
  ### Clip all not marine polygons to administrative border.
  ### Use the same data source in examples of R package wdpar.
  raw_adm_bry <- lapply(
    cnts, function(iso){
        file_path <- tempfile(fileext = ".gpkg")
        lk <- "https://geodata.ucdavis.edu/gadm/gadm4.1/gpkg"
        download.file(sprintf("%s/gadm41_%s.gpkg", lk, iso), file_path)
        read_sf(file_path, "ADM_ADM_0")}) %>% bind_rows()
  
  ### Process the boundary a bit
  adm_bry <- raw_adm_bry %>%
    st_set_precision(geometry_precision) %>%
    sf::st_make_valid() %>%
    st_set_precision(geometry_precision) %>%
    st_combine() %>%
    st_union() %>%
    st_set_precision(geometry_precision) %>%
    sf::st_make_valid() %>%
    st_transform(st_crs(clean_pas)) %>%
    sf::st_make_valid()
  
  ### Clip it to boundary
  clean_pas <- clean_pas %>% st_intersection(adm_bry)
  rm(raw_pas, adm_bry)
  ### No of PAs drop from 1638 to 1618
  
  # Separate or union polygons and re-index them.
  ## Separate or union
  ## Have to drop the attributes. Shouldn't be a problem though.
  ## QUESTION: should we consider the year of status?
  clean_pas <- st_cast(st_union(clean_pas), "POLYGON") %>% 
    st_as_sf() %>% rename(geometry = x)
  
  ## Crop it to extent
  bbox <- st_as_sfc(st_bbox(template)) %>% st_transform(st_crs(clean_pas)) %>% 
    st_buffer(bbox_buffer)
  clean_pas <- clean_pas %>% 
    slice(unique(unlist(suppressMessages(st_intersects(bbox, .))))) %>% 
    mutate(index = 1:nrow(.))
  
  ## Re-calculate the area
  clean_pas <- clean_pas %>% 
    mutate(REP_AREA = st_area(.) %>% units::set_units("km2"))
  
  ## Clean the tiny ones with "partial" marine type due to a high chance to be the
  ## noisy remaining of the clip
  false_pas <- st_join(clean_pas, raw_pas %>% select(MARINE)) %>% 
    filter(MARINE != "terrestrial" & REP_AREA < 0.01 %>% units::set_units("km2"))
  
  clean_pas <- clean_pas %>% filter(!index %in% unique(false_pas$index))
  
  ################# Cluster Protected areas (PAs) #################
  ## Cluster the PAs that are connected. 
  ## The assumption here is that it is difficult (if not possible) for mammals 
  ## to pass the ocean or huge waterbodies.
  ## Warning: this is only for mammal species
  adm_reorg <- raw_adm_bry %>% 
    st_transform(st_crs(clean_pas)) %>% 
    st_buffer(10) %>% st_union() %>% st_buffer(-10) %>% st_cast("POLYGON") %>% 
    st_intersection(bbox) %>% st_as_sf() %>% mutate(group = 1:nrow(.)) %>% 
    mutate(area = st_area(.) %>% units::set_units("km2"))
  
  clean_pas <- st_join(clean_pas, adm_reorg)
  
  ## Clean further
  false_pas <- clean_pas %>% 
    filter(area > 5e+05 %>% units::set_units("km2") & 
             REP_AREA < 1 %>% units::set_units("km2"))
  clean_pas <- clean_pas %>% filter(!index %in% unique(false_pas$index))
  
  ## Re-index
  clean_pas <- clean_pas %>% mutate(index = 1:nrow(.)) %>% 
    select(-area)
  
  adm_reorg <- adm_reorg %>% filter(group %in% clean_pas$group) %>% 
    rename(geometry = x) %>% select(group)
  
  # Save out
  st_write(clean_pas, file.path(dst_dir, "clean_pas.geojson"))
  ## Note: save out as shp because of a weird bug in geojson.
  st_write(adm_reorg, file.path(dst_dir, "pa_groups.shp"))
  
  # Switch back
  sf_use_s2(org_setting)
}
