# ====== General Information ======
# Task of the script: produce supplemental information
# Author: Wenxin Yang
# Date created: 01/19/2024

# ======= load packages & set up dir =======
library(sf)
library(dplyr)
library(terra)
library(ggplot2)
library(units)
library(cowplot)
library(nngeo)
library(ggplot2)
sf_use_s2(FALSE) # deal with buffering odd

dir_pa <- '/Users/wenxinyang/Desktop/Dissertation/DATA'

# ======= read in files ======
# read in raster data for bounding box
dir_ras <- file.path(dir_pa, 'BrodieData')
template <- rast(
    file.path(dir_ras, "GEDIv002_20190417to20220413_cover_krig.tiff")) %>% 
    extend(c(100, 100)) # add a buffer
values(template) <- 1:ncell(template)

# read in WDPA and filter to study area - copying code from calc_conn_metrics.R
pas <- read_sf(
    file.path(dir_pa, "WDPA_WDOECM_Nov2023_Public_AS",
              "WDPA_WDOECM_Nov2023_Public_AS.gdb"),
    layer = "WDPA_WDOECM_poly_Nov2023_AS")

# Q1: should we only consider terrestrial PAs?
pas <- pas %>% 
    filter(ISO3 %in% c("KHM", "CHN", "IDN", "LAO", "MYS", 
                       "SGP", "THA", "VDR", "SVR", "BRN")) %>% 
    select(WDPAID, ISO3, NAME, REP_AREA, GIS_AREA, REP_M_AREA, GIS_M_AREA) %>% 
    rename(Geometry = SHAPE) %>% 
    slice(unique(unlist(suppressMessages(st_intersects(st_as_sfc(st_bbox(template)), .))))) %>% 
    mutate(id = 1:nrow(.))

# Split MULTIPOLYGON to POLYGONS because each segment should be treated independently
# Hmm... then REP_AREA is no longer usable, so calculate the GIS area
# Q2: I think I did right thing here, didn't I?
pas <- st_cast(pas, "POLYGON") %>% st_make_valid() %>% 
    mutate(REP_AREA = st_area(.) %>% units::set_units("km2"))
# Warning: some PAs may overlap (sometimes highly).
# Wenxin: Should we dissolve them?

st_write(pas, file.path(dir_ras, 'studyarea_pa.shp'))

# ======= produce viz ======
# histogram for PA size
li_iso <- unique(pas$ISO3)
# BRN 47, IDN 920, SGP 6, CHN 13, KHM 128, MYS 978, THA 411, LAO 40

list <-lapply(1:length(li_iso),
              function(iso) 
                  ggplot(pas[pas$ISO3==li_iso[[iso]],], aes(x=REP_AREA))+
                  geom_histogram()+
                  labs(title=li_iso[[iso]], x='PA Area', y='Frequency'))

cowplot::plot_grid(plotlist = list)

ggplot(pas, aes(x=REP_AREA))+
    geom_histogram()+
    labs(title='All PA area', x='PA Area', y='Frequency')

# compute nearest neighbor for all PAs in the study area
# I used arcpro because it's much faster
pas_1 <- read_sf(file.path(dir_ras, 'studyarea_pa.shp'))
pas_1$nnd_km <- pas_1$NEAR_DIST/1000


ggplot(pas_1, aes(x=nnd_km))+
    geom_histogram(binwidth = 10)+
    labs(title='Distance to nearest neighbor for all PAs',
         x='Nearest distance (km)',
         y='Frequency')

# ======== supp 01/26 ======
dat_brodie <- data.frame(read.csv(
                          here("data/raw/public/training/bird_data_230326.csv"), 
                          header = T))
# Simplify the variable names of site identifier and geographic coordinates
names(dat_brodie)[names(dat_brodie) == "site"] <- "station"
names(dat_brodie)[names(dat_brodie) == "lat_wgs84"] <- "lat"
names(dat_brodie)[names(dat_brodie) == "long_wgs84"] <- "long"

# Search for HDI in the column names 
grep("hdi", names(dat_brodie), value = TRUE)
grep("HDI", names(dat_brodie), value = TRUE)

# Assign stations the HDI value of its country
# Reference: Human Development Report 2020: The Next Frontier—Human Development 
# and the Anthropocene (United Nations Development Programme, 2020).
# https://hdr.undp.org/data-center/human-development-index#/indicies/HDI
# https://hdr.undp.org/sites/default/files/2021-22_HDR/HDR21-22_Statistical_Annex_HDI_Table.xlsx

dat_HDI <- data.frame(
    country = unique(dat_brodie$country), 
    HDI = c(0.593, 0.768, 0.705, 0.607, 0.803, 0.939, 0.800, 
            0.703, 0.829))
dat_brodie <- left_join(dat_brodie, dat_HDI, by = "country")

# Create dataframe containing the subset of variable used in the analysis
dat <- dat_brodie %>% select(station, country, PA, long, lat, 
                             Hansen_recentloss,access_log10, HDI, dist_to_PA, 
                             PA_size_km2, rh_95_a0.pred, pavd_0_5.pred, 
                             pai_a0.pred, fhd_pai_1m_a0.pred, cover_a0.pred, 
                             agbd_a0.pred,asymptPD, maxFRic, SR.mean)

# Add the connectivity variables for each station calculated with 
# calc_conn_metrics.R
dat_conn_metrics <- data.frame(read.csv(
    here("data/derived/public/conn_flux.csv"), 
    header = T))
dat <- left_join(dat, dat_conn_metrics, by = "station")

dat100 <- dat[dat$med_dist == 100, ]
ggplot()+
    geom_histogram(aes(dat100$awf_ptg), binwidth = 100, fill="#69b3a2")

write.csv(dat100, here('data/raw/public/dat100awf.csv'))

ggplot()+
    geom_point(aes(dat100$awf_ptg, dat100$PA))
# awf_ptg


# ======= residual analysis =====
ggplot()+
    geom_point(aes(dat_matched_PD$awf_rst_ptp2.z, residuals$resid.mod_PD_efficacy.))

residuals <- data.frame(resid(mod_PD_efficacy))

a <- mod_PD_efficacy$fitted
b <- data.frame(mod_PD_efficacy$residuals)
c <- mod_PD_efficacy$data$asymptPD