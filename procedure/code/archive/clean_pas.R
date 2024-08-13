## -------------------------------------------------------------------
## Script name: clean_pas
## Purpose of script: Clean the Protected areas.
## Author: Lei Song
## Date Created: 2023-12-14
## Email: lsong@ucsb.edu

## Usage
## Rscript path/to/clean_pas.R --src_dir data/raw/public 
## --dst_dir path/for/result

## Results
## Two files to save out:
## 1. clean_pas.geojson: cleaned PAs
## 2. pa_groups.shp: PA groups for mammal connectivity
## -------------------------------------------------------------------

# Load libraries, easy to switch to use groundhog
pkgs <- c("sf", "dplyr", "terra", "optparse", "wdpar")
sapply(pkgs, require, character.only = TRUE)
sf_use_s2(FALSE) # deal with buffering odd

# Command line inputs
option_list <- list(
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

# Read samples and raster template
## According to the pairs of lat/lon and east/north, they used
## UTM Zone 46 (EPSG:32646) for projection, so here we will use the same one.
## (Not mean it is the most right one to use)
fnames <- list.files(file.path(src_dir, "training"), full.names = TRUE)
pts <- do.call(rbind, lapply(fnames, function(fname){
    read.csv(fname) %>% 
        select(all_of(names(.)[
            stringr::str_detect(names(.), "station|country|lat|long")])) %>% 
        st_as_sf(coords = c(3, 2), crs = 4326)
}))

template <- rast(
    file.path(src_dir, "GEDIv002_20190417to20220413_cover_krig.tif")) %>% 
    extend(c(100, 100)) # add a buffer
values(template) <- 1:ncell(template)

########################## Query Protected areas (PAs) #########################
## Within this part, most of the issues related to PAs are solved 
## 1. Project the PAs or relevant layers to use precise distance.
## 2. Trim the PAs to terrestrial only.
## 3. Separate or union polygons and re-index them.

# Query and clean PAs
raw_pas <- c("KHM", "CHN", "IDN", "LAO", "MYS", 
             "SGP", "THA", "VNM", "BRN") %>%
    lapply(wdpa_fetch, wait = TRUE,
           download_dir = rappdirs::user_data_dir("wdpar")) %>%
    bind_rows()
raw_pas <- wdpa_clean(raw_pas, crs = "EPSG:32646",
                      geometry_precision = 100000)
raw_pas <- raw_pas %>% filter(MARINE != "marine")
## Note: now the No is 1638

# Trim the polygons to terrestrial only
## Clip all not marine polygons to administrative border.
## Use the same data source in examples of R package wdpar.
raw_adm_bry <- lapply(
    c("KHM", "CHN", "IDN", "LAO", "MYS", 
      "SGP", "THA", "VNM", "BRN"), function(iso){
          file_path <- tempfile(fileext = "rds")
          lk <- "https://biogeo.ucdavis.edu/data/gadm3.6/Rsf"
          download.file(sprintf("%s/gadm36_%s_0_sf.rds", lk, iso), file_path)
          readRDS(file_path)}) %>% bind_rows()

## Process the boundary a bit
adm_bry <- raw_adm_bry %>%
    st_set_precision(100000) %>%
    sf::st_make_valid() %>%
    st_set_precision(100000) %>%
    st_combine() %>%
    st_union() %>%
    st_set_precision(100000) %>%
    sf::st_make_valid() %>%
    st_transform(st_crs(raw_pas)) %>%
    sf::st_make_valid()

## Clip it to boundary
clean_pas <- raw_pas %>% st_intersection(adm_bry)
## No of PAs drop from 1638 to 1618

# Crop it to extent
## WARNING: remember to re-run this step again after union/separate the polygons
bbox <- st_as_sfc(st_bbox(template)) %>% st_transform(st_crs(clean_pas))
clean_pas <- clean_pas %>% 
    slice(unique(unlist(suppressMessages(st_intersects(bbox, .)))))
## Note: No of PAs drop further to 1260

# Union/separate polygons.
# Have to drop the attributes. Shouldn't be a problem though.
clean_pas <- st_cast(st_union(clean_pas), "POLYGON") %>% 
    st_as_sf() %>% rename(geometry = x)
## No of PAs change from 1260 to 4270 (lose administrative meaning)

# Crop to extent again
clean_pas <- clean_pas %>% 
    slice(unique(unlist(suppressMessages(st_intersects(bbox, .))))) %>% 
    mutate(index = 1:nrow(.))
## No of PAs change from 4270 to 4259

# Re-calculate the area
clean_pas <- clean_pas %>% 
    mutate(REP_AREA = st_area(.) %>% units::set_units("km2"))

# Clean the tiny ones with "partial" marine type due to a high chance to be the
# noisy remaining of the clip
false_pas <- st_join(clean_pas, raw_pas %>% select(MARINE)) %>% 
    filter(MARINE != "terrestrial" & REP_AREA < 0.01 %>% units::set_units("km2"))

clean_pas <- clean_pas %>% filter(!index %in% unique(false_pas$index))

######################### Cluster Protected areas (PAs) ########################
## 4. Cluster the PAs that are connected. 
## The assumption here is that it is difficult (if not possible) for mammals 
## to pass the ocean or huge waterbodies.
## Warning: this is only for mammal species
adm_reorg <- raw_adm_bry %>% 
    st_transform(st_crs(clean_pas)) %>% 
    st_buffer(10) %>% st_union() %>% st_buffer(-10) %>% st_cast("POLYGON") %>% 
    st_intersection(bbox) %>% st_as_sf() %>% mutate(group = 1:nrow(.)) %>% 
    mutate(area = st_area(.) %>% units::set_units("km2"))

clean_pas <- st_join(clean_pas, adm_reorg) # No is 2873

## Clean further
false_pas <- clean_pas %>% 
    filter(area > 5e+05 %>% units::set_units("km2") & 
               REP_AREA < 1 %>% units::set_units("km2"))
clean_pas <- clean_pas %>% filter(!index %in% unique(false_pas$index)) # 2198

## Re-index
clean_pas <- clean_pas %>% mutate(index = 1:nrow(.)) %>% 
    select(-area)

adm_reorg <- adm_reorg %>% filter(group %in% clean_pas$group) %>% 
    rename(geometry = x) %>% select(group)

# Save out
st_write(clean_pas, file.path(dst_dir, "clean_pas.geojson"))
st_write(adm_reorg, file.path(dst_dir, "pa_groups.shp"))
