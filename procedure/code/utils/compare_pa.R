## -------------------------------------------------------------------
## Specification checks
## -------------------------------------------------------------------

# Load libraries
pkgs <- c("sf", "dplyr", "terra", "stringr", "ggplot2", 
          "tidyverse", "rnaturalearth")
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
pa_groups <- st_read(file.path(dst_dir, "pa_groups.shp"))

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

# Query country boundaries
cnts <- ne_countries(continent = "Asia") %>% 
    st_transform(st_crs(pts_all)) %>% 
    filter(subunit %in% unique(pts_all$country)) %>% 
    st_intersection(st_convex_hull(st_union(clean_pas)))

# Add BigPA and CloseToPA,
# and keep values to just outside of PAs.
pts_all <- pts_all %>% 
    mutate(PA_size_km2_ours = ifelse(PA_ours == 1, NA, PA_size_km2_ours),
           PA_size_km2 = ifelse(PA == 1, NA, PA_size_km2),
           dist_to_PA_ours = ifelse(PA_ours == 1, NA, dist_to_PA_ours),
           dist_to_PA = ifelse(PA == 1, NA, dist_to_PA)) %>% 
    mutate(BigPA_ours = ifelse(PA_ours == 1, NA, ifelse(PA_size_km2_ours < 500, 0, 1)),
           BigPA = ifelse(PA == 1, NA, ifelse(PA_size_km2 < 500, 0, 1))) %>% 
    mutate(CloseToPA_ours = ifelse(PA_ours == 1, NA, ifelse(dist_to_PA_ours > 2, 0, 1)),
           CloseToPA = ifelse(PA == 1, NA, ifelse(dist_to_PA > 2, 0, 1)))

# Compare with original values
pa_sizes <- pts_all %>% st_drop_geometry() %>% 
    select(taxon, PA_size_km2_ours, PA_size_km2) %>% 
    pivot_longer(2:3, names_to = "Type") %>% 
    mutate(Type = factor(Type, levels = c("PA_size_km2", "PA_size_km2_ours"),
                         labels = c("Original", "Updated")))

ggdensity(pa_sizes, x = "value",
          add = "mean", rug = TRUE,
          color = "Type", fill = "Type",
          palette = c("#00AFBB", "#E7B800"),
          na.rm = TRUE) +
    xlab("PA size (square km)") +
    ylab("Density") +
    theme(legend.text = element_text(size = 12, color = "black"))
     
ggsave("results/figures/pa_size_compare.png", width = 6, height = 6)

pa_dists <- pts_all %>% st_drop_geometry() %>% 
    select(taxon, dist_to_PA_ours, dist_to_PA) %>% 
    pivot_longer(2:3, names_to = "Type") %>% 
    mutate(Type = factor(Type, levels = c("dist_to_PA", "dist_to_PA_ours"),
                         labels = c("Original", "Updated")))

ggdensity(pa_dists, x = "value",
          add = "mean", rug = TRUE,
          color = "Type", fill = "Type",
          palette = c("#00AFBB", "#E7B800"),
          na.rm = TRUE) +
    xlab("Distance to PA (km)") +
    ylab("Density") +
    theme(legend.text = element_text(size = 12, color = "black"))
ggsave("results/figures/dist_pa_compare.png", width = 6, height = 6)

# PA or not PA
sum(pts_all$PA_ours != pts_all$PA) / nrow(pts_all) * 100

vals <- pts_all %>% st_drop_geometry() %>% 
    select(PA, PA_ours) %>% 
    mutate(PA = as.factor(PA), PA_ours = as.factor(PA_ours))

cm <- confusionMatrix(vals$PA, vals$PA_ours)
cm$table[1, 2] / nrow(pts_all) * 100 # from 1 to 0
cm$table[2, 1] / nrow(pts_all) * 100 # from 0 to 1

# BigPA or not BigPA
sum(paste(pts_all$BigPA_ours) != paste(pts_all$BigPA)) / nrow(pts_all) * 100

vals <- pts_all %>% st_drop_geometry() %>% 
    select(BigPA, BigPA_ours) %>% 
    mutate(BigPA = as.factor(paste(BigPA)), 
           BigPA_ours = as.factor(paste(BigPA_ours)))

cm <- confusionMatrix(vals$BigPA, vals$BigPA_ours)
cm$table[1, 2] / nrow(pts_all) * 100 # from 1 to 0
cm$table[3, 2] / nrow(pts_all) * 100 # from 1 to NA
cm$table[2, 1] / nrow(pts_all) * 100 # from 0 to 1
cm$table[3, 1] / nrow(pts_all) * 100 # from 0 to NA
cm$table[1, 3] / nrow(pts_all) * 100 # from NA to 0
cm$table[2, 3] / nrow(pts_all) * 100 # from NA to 1

# CloseToPA or not CloseToPA
sum(paste(pts_all$CloseToPA_ours) != paste(pts_all$CloseToPA)) / nrow(pts_all) * 100

vals <- pts_all %>% st_drop_geometry() %>% 
    select(CloseToPA, CloseToPA_ours) %>% 
    mutate(CloseToPA = as.factor(paste(CloseToPA)), 
           CloseToPA_ours = as.factor(paste(CloseToPA_ours)))

cm <- confusionMatrix(vals$CloseToPA, vals$CloseToPA_ours)
cm$table[1, 2] / nrow(pts_all) * 100 # from 1 to 0
cm$table[3, 2] / nrow(pts_all) * 100 # from 1 to NA
cm$table[2, 1] / nrow(pts_all) * 100 # from 0 to 1
cm$table[3, 1] / nrow(pts_all) * 100 # from 0 to NA
cm$table[1, 3] / nrow(pts_all) * 100 # from NA to 0
cm$table[2, 3] / nrow(pts_all) * 100 # from NA to 1
