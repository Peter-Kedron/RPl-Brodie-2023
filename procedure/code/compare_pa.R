## -------------------------------------------------------------------
## Specification checks
## -------------------------------------------------------------------

# Load libraries
pkgs <- c("sf", "dplyr", "terra", "stringr", "ggplot2")
sapply(pkgs, require, character.only = TRUE)

# Set directories
src_dir <- "data/raw/public"
dst_dir <- 'data/derived/public'

# Load datasets
pts_bird <- st_read(file.path(dst_dir, "pts_bird.geojson"))
pts_mammal <- st_read(file.path(dst_dir, "pts_mammal.geojson"))
pts_all <- rbind(pts_bird %>% mutate(taxon = "bird"), 
                 pts_mammal %>% mutate(taxon = "mammal"))
clean_pas <- st_read(file.path(dst_dir, "clean_pas.geojson"))

########################### Check PA attributes ################################
# Check if point is in/out PAs
pts_all <- st_join(pts_all, clean_pas %>% select(index))
pts_all <- pts_all %>% mutate(PA_ours = ifelse(is.na(index), 0, 1)) %>% 
    select(-index)

# Check if point is not/close to big PA
ids_pa <- st_nearest_feature(pts_all, clean_pas)
dists <- sapply(1:nrow(pts_all), function(id){
    st_distance(pts_all[id, ], clean_pas[ids_pa[id], ]) %>% units::set_units("km")
})

pts_all <- pts_all %>% 
    mutate(dist_to_PA_ours = dists,
           PA_size_km2_ours = clean_pas[ids_pa, ]$REP_AREA)

# Load Brodie data
dat_brodie <- read.csv(
    file.path(src_dir, "training/bird_data_230326.csv"), 
    header = T) %>% select(station, PA, dist_to_PA, PA_size_km2) %>% 
    rbind(read.csv(
        file.path(src_dir, "training/mammal_data_230326.csv"), 
        header = T) %>% select(station, PA, dist_to_PA, PA_size_km2))

pts_all <- left_join(pts_all, dat_brodie, by = "station")

# Add BigPA and CloseToPA
pts_all <- pts_all %>% 
    mutate(BigPA_ours = ifelse(PA_size_km2_ours < 500, 0, 1),
           BigPA = ifelse(PA_size_km2 < 500, 0, 1)) %>% 
    mutate(CloseToPA_ours = ifelse(dist_to_PA_ours > 2, 0, 1),
           CloseToPA = ifelse(dist_to_PA > 2, 0, 1))

st_write(pts_all, file.path(dst_dir, 'PA_attributes_compare.geojson'))

################## Check Connectivity index distribution #######################
## Use our calculation to keep consistent.

# Read flux
awf <- lapply(c("bird", "mammal"), function(taxon){
    read.csv(file.path(dst_dir, sprintf("conn_flux_%s_10_160.csv", taxon))) %>% 
        select(station, awf_ptg, med_dist)
}) %>% bind_rows()

awf <- left_join(awf, pts_all %>% st_drop_geometry() %>% 
                     select(station, taxon, PA_ours),
                 by = "station")
ggplot(awf, aes(x = as.factor(med_dist), y = awf_ptg)) + 
    geom_boxplot(aes(fill = as.factor(PA_ours))) +
    # geom_jitter(shape = 16, position = position_jitter(0.2)) +
    facet_wrap(~taxon, nrow = 2)

