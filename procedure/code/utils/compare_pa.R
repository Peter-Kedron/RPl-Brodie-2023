## -------------------------------------------------------------------
## Specification checks
## -------------------------------------------------------------------

# Load libraries
pkgs <- c("sf", "dplyr", "terra", "stringr", "ggplot2", "tidyverse")
sapply(pkgs, require, character.only = TRUE)

# Set directories
src_dir <- "data/raw/public"
dst_dir <- 'data/derived/public'

# Load datasets
pts_bird <- st_read(file.path(dst_dir, "pts_bird.geojson"))
pts_mammal <- st_read(file.path(dst_dir, "pts_mammal.geojson"))
clean_pas <- st_read(file.path(dst_dir, "clean_pas.geojson"))
pa_groups <- st_read(file.path(dst_dir, "pa_groups.shp"))

########################### Check PA attributes ################################
# Check if point is not/close to big PA
## bird
ids_pa <- st_nearest_feature(pts_bird, clean_pas)
dists <- sapply(1:nrow(pts_bird), function(id){
    st_distance(pts_bird[id, ], clean_pas[ids_pa[id], ]) %>% units::set_units("km")
})

pts_bird <- pts_bird %>% 
    mutate(taxon = "bird",
           dist_to_PA_ours = dists,
           PA_size_km2_ours = clean_pas[ids_pa, ]$REP_AREA)
pts_bird <- st_join(pts_bird, clean_pas %>% select(index))
pts_bird <- pts_bird %>% mutate(PA_ours = ifelse(is.na(index), 0, 1)) %>% 
    select(-index)

## Mammal
pts_mammal <- pts_mammal %>% st_join(pa_groups)
pts_mammal <- lapply(1:nrow(pts_mammal), function(i){
    pts_this <- pts_mammal %>% slice(i)
    pa_this <- clean_pas %>% filter(group %in% pts_this$group)
    id_this <- st_nearest_feature(pts_this, pa_this)
    
    dist <- st_distance(pts_this, pa_this[id_this, ]) %>% units::set_units("km")
    
    pts_this <- pts_this %>% 
        mutate(taxon = "mammal",
               dist_to_PA_ours = dist,
               PA_size_km2_ours = pa_this[id_this, ]$REP_AREA)
    
    pts_this <- st_join(pts_this, pa_this %>% select(index))
    pts_this %>% mutate(PA_ours = ifelse(is.na(index), 0, 1)) %>% 
        select(-index)
}) %>% bind_rows()

# Check if point is in/out PAs
pts_all <- rbind(pts_bird, pts_mammal %>% select(-group))

# Load Brodie data
dat_brodie <- read.csv(
    file.path(src_dir, "training/bird_data_corrected_240122.csv"), 
    header = T) %>% select(station, PA, dist_to_PA, PA_size_km2) %>% 
    rbind(read.csv(
        file.path(src_dir, "training/mammal_data_corrected_240122.csv"), 
        header = T) %>% select(station, PA, dist_to_PA, PA_size_km2))

pts_all <- left_join(pts_all, dat_brodie, by = "station")

# Add BigPA and CloseToPA
pts_all <- pts_all %>% 
    mutate(BigPA_ours = ifelse(PA_size_km2_ours < 500, 0, 1),
           BigPA = ifelse(PA_size_km2 < 500, 0, 1)) %>% 
    mutate(CloseToPA_ours = ifelse(dist_to_PA_ours > 2, 0, 1),
           CloseToPA = ifelse(dist_to_PA > 2, 0, 1))

################## Check Connectivity index distribution #######################
# Compare with original values
pa_sizes <- pts_all %>% st_drop_geometry() %>% 
    select(taxon, PA_ours, PA, PA_size_km2_ours, PA_size_km2) 
pa_sizes <- data.frame(value = pa_sizes$PA_size_km2_ours[pa_sizes$PA_ours == 0],
                  type = "PA_size_km2_ours",
                  taxon = pa_sizes$taxon[pa_sizes$PA_ours == 0]) %>% 
    rbind(data.frame(value = pa_sizes$PA_size_km2[pa_sizes$PA == 0],
                     type = "PA_size_km2",
                     taxon = pa_sizes$taxon[pa_sizes$PA == 0])) %>% 
    mutate(type = factor(type, levels = c("PA_size_km2", "PA_size_km2_ours"),
                         labels = c("Original", "Updated")))

ggplot(pa_sizes, aes(x = taxon, y = value, fill = type)) +
    geom_boxplot(outliers = FALSE) +
    scale_fill_brewer(name = "", palette = "Dark2") +
    xlab("Taxon group") +
    ylab("PA size (square km)") +
    theme_classic()+
    theme(panel.spacing.y = unit(2, "lines"),
          axis.text = element_text(size = 12, color = "black"),
          axis.title = element_text(size = 12, color = "black"),
          legend.text = element_text(size = 12, color = "black"),
          legend.position = "top")
ggsave("results/figures/pa_size_compare.png", width = 6, height = 6)

pa_dists <- pts_all %>% st_drop_geometry() %>% 
    select(taxon, PA_ours, PA, dist_to_PA_ours, dist_to_PA)

pa_dists <- data.frame(value = pa_dists$dist_to_PA_ours[pa_dists$PA_ours == 0],
                       type = "dist_to_PA_ours",
                       taxon = pa_dists$taxon[pa_dists$PA_ours == 0]) %>% 
    rbind(data.frame(value = pa_dists$dist_to_PA[pa_dists$PA == 0],
                     type = "dist_to_PA",
                     taxon = pa_dists$taxon[pa_dists$PA == 0])) %>% 
    mutate(type = factor(type, levels = c("dist_to_PA", "dist_to_PA_ours"),
                         labels = c("Original", "Updated")))

ggplot(pa_dists, aes(x = taxon, y = value, fill = type)) +
    geom_boxplot(outliers = FALSE) +
    scale_fill_brewer(name = "", palette = "Dark2") +
    xlab("Taxon group") +
    ylab("Distance to PA (km)") +
    theme_classic()+
    theme(panel.spacing.y = unit(2, "lines"),
          axis.text = element_text(size = 12, color = "black"),
          axis.title = element_text(size = 12, color = "black"),
          legend.text = element_text(size = 12, color = "black"),
          legend.position = "top")
ggsave("results/figures/dist_pa_compare.png", width = 6, height = 6)

pts_all %>% st_drop_geometry() %>% 
    select(taxon, PA_ours, PA) %>% 
    mutate(change = PA_ours - PA) %>% 
    group_by(change) %>% 
    summarise(n = n())

pts_all %>% st_drop_geometry() %>% 
    select(taxon, PA_ours, PA) %>% 
    pivot_longer(2:3, names_to = "type") %>%
    group_by(taxon, type, value) %>% 
    summarise(n = n()) %>% 
    mutate(type = factor(type, levels = c("PA", "PA_ours"),
                         labels = c("Original", "Updated")))
# Not too different. Good! Likely the increase in bird because of the added 
# circular PAs. Decrease in mammals results from the reduction because of 
# not same group.

big_pa <- pts_all %>% st_drop_geometry() %>% 
    select(taxon, PA_ours, PA, BigPA_ours, BigPA)

big_pa <- data.frame(value = big_pa$BigPA_ours[big_pa$PA_ours == 0],
                       type = "BigPA_ours",
                       taxon = big_pa$taxon[big_pa$PA_ours == 0]) %>% 
    rbind(data.frame(value = big_pa$BigPA[big_pa$PA == 0],
                     type = "BigPA",
                     taxon = big_pa$taxon[big_pa$PA == 0])) %>% 
    group_by(taxon, type, value) %>% 
    summarise(n = n()) %>%
    mutate(type = factor(type, levels = c("BigPA", "BigPA_ours"),
                         labels = c("Original", "Updated")))
# 1 increase, make sense

closeto_pa <- pts_all %>% st_drop_geometry() %>% 
    select(taxon, PA_ours, PA, CloseToPA_ours, CloseToPA)

closeto_pa <- data.frame(value = closeto_pa$CloseToPA_ours[closeto_pa$PA_ours == 0],
                     type = "CloseToPA_ours",
                     taxon = closeto_pa$taxon[closeto_pa$PA_ours == 0]) %>% 
    rbind(data.frame(value = closeto_pa$CloseToPA[closeto_pa$PA == 0],
                     type = "CloseToPA",
                     taxon = closeto_pa$taxon[closeto_pa$PA == 0])) %>% 
    group_by(taxon, type, value) %>% 
    summarise(n = n()) %>%
    mutate(type = factor(type, levels = c("CloseToPA", "CloseToPA_ours"),
                         labels = c("Original", "Updated")))
# Consistent with PA, good
