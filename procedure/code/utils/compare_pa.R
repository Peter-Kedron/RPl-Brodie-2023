## -------------------------------------------------------------------
## Specification checks
## -------------------------------------------------------------------

# Load libraries
pkgs <- c("sf", "dplyr", "terra", "stringr", "ggplot2", 
          "tidyverse", "rnaturalearth", "ggpubr", "caret", "colorspace")
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

v1 <- ggviolin(pa_sizes, x = "Type", y = "value", fill = "Type",
               palette = c("#00AFBB", "#E7B800"),
               add.params = list(fill = "white"), add = 'mean_ci') + 
    theme(legend.text = element_text(size = 12, color = "black")) +
    ylim(0, 3000) + ylab("PA size (square km)") + xlab("")

pa_dists <- pts_all %>% st_drop_geometry() %>% 
    select(taxon, dist_to_PA_ours, dist_to_PA) %>% 
    pivot_longer(2:3, names_to = "Type") %>% 
    mutate(Type = factor(Type, levels = c("dist_to_PA", "dist_to_PA_ours"),
                         labels = c("Original", "Updated")))

v2 <- ggviolin(pa_dists, x = "Type", y = "value", fill = "Type",
               palette = c("#00AFBB", "#E7B800"),
               add.params = list(fill = "white"), add = 'mean_ci') +
    theme(legend.text = element_text(size = 12, color = "black")) +
    ylab("Distance to PA (square km)") + xlab("")

ggarrange(v1, v2, ncol = 2, common.legend = TRUE)

ggsave("results/figures/violin_pa_size_dist.png", 
       width = 6, height = 4, bg = "white")

# PA or not PA
sum(pts_all$PA_ours != pts_all$PA) / nrow(pts_all) * 100

vals <- pts_all %>% st_drop_geometry() %>% 
    select(PA, PA_ours) %>% 
    mutate(PA = as.factor(PA), PA_ours = as.factor(PA_ours))
cm <- confusionMatrix(vals$PA, vals$PA_ours, dnn = c("Updated", "Original"))

# stats
cm$table[1, 2] / nrow(pts_all) * 100 # from 1 to 0
cm$table[2, 1] / nrow(pts_all) * 100 # from 0 to 1

# Figure
plt <- as.data.frame(cm$table)
levels(plt$Updated) <- c("Outside PA", "Inside PA")
levels(plt$Original) <- c("Outside PA", "Inside PA")
plt$Original <- factor(plt$Original, levels = rev(levels(plt$Original)))

c1 <- ggplot(plt, aes(x = Original, y = Updated, fill= Freq)) +
    geom_tile(color = "white") + coord_equal() +
    geom_text(aes(Original, Updated, label = round(Freq, 0)), 
              color = "black", size = 4, fontface = "bold") +
    scale_fill_continuous_sequential(name = "Number", palette = "Heat") +
    xlab("Brodie et al. (2023)") + ylab("Reanalysis") +
    theme_minimal() +
    theme(panel.grid.major = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.ticks = element_blank(),
          legend.spacing.y = unit(0.1, 'cm'),
          legend.margin = margin(rep(0, 4)), 
          text = element_text(size = 10, color = "black"),
          legend.text = element_text(size = 10),
          legend.title = element_text(size = 10),
          axis.text.x = element_text(color= "black", size = 10),
          axis.text.y = element_text(size = 10, color = "black", 
                                     angle = 90, hjust = 0.5),
          axis.title = element_text(size = 10, face = "italic", color = "blue"),
          axis.title.x = element_text(vjust = -2),
          axis.title.y = element_text(vjust = 2))

# BigPA or not BigPA
sum(paste(pts_all$BigPA_ours) != paste(pts_all$BigPA)) / nrow(pts_all) * 100

vals <- pts_all %>% st_drop_geometry() %>% 
    select(BigPA, BigPA_ours) %>% 
    mutate(BigPA = as.factor(paste(BigPA)), 
           BigPA_ours = as.factor(paste(BigPA_ours)))

cm <- confusionMatrix(vals$BigPA, vals$BigPA_ours, 
                      dnn = c("Updated", "Original"))

# stats
cm$table[1, 2] / nrow(pts_all) * 100 # from 1 to 0
cm$table[3, 2] / nrow(pts_all) * 100 # from 1 to NA
cm$table[2, 1] / nrow(pts_all) * 100 # from 0 to 1
cm$table[3, 1] / nrow(pts_all) * 100 # from 0 to NA
cm$table[1, 3] / nrow(pts_all) * 100 # from NA to 0
cm$table[2, 3] / nrow(pts_all) * 100 # from NA to 1

# Figure
plt <- as.data.frame(cm$table)
levels(plt$Updated) <- c(expression("<"~500~km^2), expression(">="~500~km^2), "Inside~PA")
levels(plt$Original) <- c(expression("<"~500~km^2), expression(">="~500~km^2), "Inside~PA")
plt$Updated <- factor(plt$Updated, levels = rev(levels(plt$Updated)))

c2 <- ggplot(plt, aes(Original, Updated, fill= Freq)) +
    geom_tile(color = "white") + coord_equal() + 
    geom_text(aes(Original, Updated, label = round(Freq, 0)), 
              color = "black", size = 4, fontface = "bold") +
    scale_fill_continuous_sequential(name = "Number", palette = "Heat") +
    xlab("Brodie et al. (2023)") + ylab("Reanalysis") +
    scale_x_discrete(labels= str2expression(levels(plt$Original))) + 
    scale_y_discrete(labels= str2expression(levels(plt$Updated))) +
    theme_minimal() +
    theme(panel.grid.major = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.ticks = element_blank(),
          legend.position = 'bottom',
          legend.direction = "horizontal",
          legend.spacing.y = unit(0.1, 'cm'),
          legend.margin = margin(rep(0, 4)), 
          text = element_text(size = 10, color = "black"),
          legend.text = element_text(size = 10),
          legend.title = element_text(size = 10),
          axis.text.x = element_text(color= "black", size = 10),
          axis.text.y = element_text(size = 10, color = "black", 
                                     angle = 90, hjust = 0.5),
          axis.title = element_text(size = 10, face = "italic", color = "blue"),
          axis.title.x = element_text(vjust = -2),
          axis.title.y = element_text(vjust = 2))

# CloseToPA or not CloseToPA
sum(paste(pts_all$CloseToPA_ours) != paste(pts_all$CloseToPA)) / nrow(pts_all) * 100

vals <- pts_all %>% st_drop_geometry() %>% 
    select(CloseToPA, CloseToPA_ours) %>% 
    mutate(CloseToPA = as.factor(paste(CloseToPA)), 
           CloseToPA_ours = as.factor(paste(CloseToPA_ours)))

cm <- confusionMatrix(vals$CloseToPA, vals$CloseToPA_ours,
                      dnn = c("Updated", "Original"))

# stats
cm$table[1, 2] / nrow(pts_all) * 100 # from 1 to 0
cm$table[3, 2] / nrow(pts_all) * 100 # from 1 to NA
cm$table[2, 1] / nrow(pts_all) * 100 # from 0 to 1
cm$table[3, 1] / nrow(pts_all) * 100 # from 0 to NA
cm$table[1, 3] / nrow(pts_all) * 100 # from NA to 0
cm$table[2, 3] / nrow(pts_all) * 100 # from NA to 1

# Figure
plt <- as.data.frame(cm$table)
levels(plt$Updated) <- c("< 2 km", ">= 2 km", "Inside PA")
levels(plt$Original) <- c("< 2 km", ">= 2 km", "Inside PA")
plt$Updated <- factor(plt$Updated, levels = rev(levels(plt$Updated)))

c3 <- ggplot(plt, aes(Original, Updated, fill= Freq)) +
    geom_tile(color = "white") + coord_equal() + labs(x = "", y = "") +
    geom_text(aes(Updated, Original, label = round(Freq, 0)), 
              color = "black", size = 4, fontface = "bold") +
    scale_fill_continuous_sequential(name = "Number", palette = "Heat") +
    theme_minimal() +
    xlab("Brodie et al. (2023)") + ylab("Reanalysis") +
    theme_minimal() +
    theme(panel.grid.major = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.ticks = element_blank(),
          legend.spacing.y = unit(0.1, 'cm'),
          legend.margin = margin(rep(0, 4)), 
          text = element_text(size = 10, color = "black"),
          legend.text = element_text(size = 10),
          legend.title = element_text(size = 10),
          axis.text.x = element_text(color= "black", size = 10),
          axis.text.y = element_text(size = 10, color = "black", 
                                     angle = 90, hjust = 0.5),
          axis.title = element_text(size = 10, face = "italic", color = "blue"),
          axis.title.x = element_text(vjust = -2),
          axis.title.y = element_text(vjust = 2))

ggarrange(c1, c2, c3, ncol = 3, common.legend = TRUE, legend = "right")

ggsave("results/figures/pa_change_confusion_matrics.png", 
       width = 9, height = 3, bg = "white")
