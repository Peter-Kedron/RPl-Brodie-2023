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
    read.csv(file.path(dst_dir, sprintf("conn_flux_%s_10_150.csv", taxon))) %>% 
        select(station, awf_ptg, med_dist)
}) %>% bind_rows()

awf <- left_join(awf, pts_all %>% st_drop_geometry() %>% 
                     select(station, taxon, PA_ours),
                 by = "station")

write.csv(awf, file.path(dst_dir, "awf_ptg_collection.csv"), row.names = FALSE)

# Check the distribution
ggplot(awf, aes(x = as.factor(med_dist), y = awf_ptg), outlier.size = 0.8) + 
    geom_boxplot(aes(fill = as.factor(PA_ours))) +
    facet_wrap(~taxon, nrow = 2) +
    xlab("Median dispersal distance (km)") +
    ylab("Area weighted flux (point to polygon)") +
    labs(fill = "Inside PA\n(0|1, out | in)") +
    theme_light() +
    theme(text = element_text(size = 10),
          strip.text = element_text(
              size = 12, color = "blue"))
ggsave("sd_awf_pa.png", width = 8, height = 10)

# PAs distribution
ggplot() +
    geom_sf(data = pa_groups, aes(fill = as.factor(group)), 
            lwd = 0, show.legend = FALSE) +
    geom_sf(data = clean_pas, fill = "transparent", color = 'black') +
    geom_sf(data = pts_bird %>% mutate(Taxon = "Bird") %>% 
                rbind(pts_mammal %>% mutate(Taxon = "Mammal")), 
            aes(color = Taxon), size = 0.6) +
    scale_color_manual(values = c("yellow", "blue")) +
    theme_bw()+
    theme(text = element_text(size = 10))
ggsave("pas_stations.png", width = 8, height = 8)
