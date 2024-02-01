# ====== General Information ======
# Task of the script: produce supplemental information
# Author: Wenxin Yang
# Date created: 01/19/2024

# ======= load packages & set up dir =======
pkgs <- c('sf', 'dplyr', 'terra', 'ggplot2', 'units', 'cowplot', 'nngeo',
          'tidyverse', 'hrbrthemes', 'viridis', 'here', 'spdep', 'tmap')
lapply(pkgs, library, character.only=TRUE)
sf_use_s2(FALSE) # deal with buffering odd

dir_pa <- '/Users/wenxinyang/Desktop/Dissertation/DATA'

# ======= prep ======
modelfolder <- '/Users/wenxinyang/Desktop/Dissertation/DATA/BrodieData/model_data'
clean_pas <- st_read(file.path(modelfolder, "clean_pas.geojson"))
pas <- clean_pas

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
# note: studyarea_pa.shp is not the latest version, clean_pa is
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
dat_brodie_mammal <- data.frame(read.csv(
    here("data/raw/public/training/mammal_data_230326.csv"),
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

# create a dataframe to get geometry
dat_geom_bird <- dat_brodie %>% select(long, lat, utm_east, utm_north)
dat_geom_mammal <- dat_brodie_mammal %>% select(long_wgs84, lat_wgs84, utm_east, utm_north)
colnames(dat_geom_mammal) <- c('long', 'lat', 'utm_east', 'utm_north')
dat_for_geom <- rbind(dat_geom_bird, dat_geom_mammal)
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

# write.csv(dat100, here('data/raw/public/dat100awf.csv'))

ggplot()+
    geom_point(aes(dat100$awf_ptg, dat100$PA))
# awf_ptg


# ======= residual analysis 01/31 =====
# ======= 1) spatial patterns of Brodie residuals =========
models <- list.files(path=modelfolder, pattern = "\\.rda$")
# rdmodels <- lapply(models, function(x) load(file=file.path(modelfolder, x)))

load_object <- function(file) {
    tmp <- new.env()
    load(file = file, envir = tmp)
    tmp[[ls(tmp)[1]]]
}

getBrodieDf <- function(species, ind){
    
    modname <- paste(species, ind, "models.rda", sep="_")
    modpath <- file.path(modelfolder, modname)
    mod <- load_object(modpath)
    
    bromod <- mod$brodie
    df <- cbind(bromod$data, data.frame(resid(bromod)))
    df$PA <- as.factor(df$PA)
    
    df <- merge(df, dat_for_geom, by=c('utm_east', 'utm_north'), all.x = TRUE)
    
    df$species <- species
    df$ind <- ind
    
    return(df)
}

getBrodieAllDfSpecies <- function(species){
    dfspecies <- data.frame(matrix(ncol=19, nrow=0))
    for(ind in c('sr', 'fr', 'pd')){
        dfind <- getBrodieDf(species, ind)
        dfspecies <- rbind(dfspecies, dfind)
    }
    dfspecies$category <- as.factor(paste(dfspecies$ind, dfspecies$PA, sep='-'))
    return(dfspecies)
}

getMoransI <- function(dataf){
    
    df_moran <- data.frame(matrix(ncol=5, nrow=0))
    
    for(ind in c('sr', 'fr', 'pd')){
        dftmp <- dataf[dataf$ind == ind, ]
        df.dist <- as.matrix(dist(cbind(dftmp$long, dftmp$lat)))
        df.dist.inv <- 1/df.dist
        diag(df.dist.inv) <- 0
        df_moran_ind <- data.frame(Moran.I(dftmp$resid.bromod., df.dist.inv))
        df_moran_ind$ind <- ind
        df_moran <- rbind(df_moran, df_moran_ind)
    }

    return(df_moran)
}

getMap <- function(species, ind){

    df <- getBrodieDf(species, ind)
    geodf <- st_as_sf(df, coords=c("long", "lat"), crs="epsg:4326")
   
     # it's taking too long to produce a map with PA boundary
    ggplot() +
        geom_sf(data=geodf[geodf$resid.bromod.<400, ], # can remove outlier to improve visual effect
                aes(color=resid.bromod., shape=PA), size=0.5)+
        scale_shape_manual(values=c(0, 1))+
        scale_color_viridis()+
        # geom_sf(data = clean_pas, fill = "transparent", color = 'black')+
        theme_bw()
}


# analyze residuals
BrodieMammal <- getBrodieAllDfSpecies('mammal')
# map
BrodieMammal %>%
    ggplot( aes(x=category, y=resid.bromod., fill=PA)) +
    geom_boxplot() +
    scale_fill_viridis(discrete = TRUE, alpha=0.6) +
    # geom_jitter(color="black", size=0.4, alpha=0.9) +
    theme_ipsum() +
    theme(
        legend.position="none",
        plot.title = element_text(size=11)
    ) +
    ggtitle(paste("Brodie mammal")) +
    xlab("Outside (0) or within PA (1)") +
    ylab("Residuals")
# summary stats
tapply(BrodieMammal$resid.bromod., BrodieMammal$category, summary)
# Moran's I
BrodieMammalMoran <- getMoransI(BrodieMammal)
BrodieMammalMoran$species <- 'mammal'


BrodieBird <- getBrodieAllDfSpecies('bird')
BrodieBird %>%
    ggplot( aes(x=category, y=resid.bromod., fill=PA)) +
    geom_boxplot() +
    scale_fill_viridis(discrete = TRUE, alpha=0.6) +
    # geom_jitter(color="black", size=0.4, alpha=0.9) +
    theme_ipsum() +
    theme(
        legend.position="none",
        plot.title = element_text(size=11)
    ) +
    ggtitle(paste("Brodie bird")) +
    xlab("Outside (0) or within PA (1)") +
    ylab("Residuals")
tapply(BrodieBird$resid.bromod., BrodieBird$category, summary)
BrodieBirdMoran <- getMoransI(BrodieBird)
BrodieBirdMoran$species <- 'bird'

BrodieMoran <- rbind(BrodieBirdMoran, BrodieMammalMoran)

# ======= 2) Residuals of all models ======
# create a dataframe that stores all the residuals for all 
# for birds
species <- 'bird'
ind <- 'sr'
mod <- load_object(file.path(modelfolder, 'bird_sr_models.rda'))

getResidSpeciesInd <- function(species, ind){
    modname <- paste(species, ind, 'models.rda', sep='_')
    modpath <- file.path(modelfolder, modname)
    mods <- load_object(modpath)
    
    li_modnames <- names(mods)
    
    if(species == 'bird'){
        dfresid <- dat_geom_bird[c('long', 'lat')]
    }else{
        dfresid <- dat_geom_mammal[c('long', 'lat')]
    }
    
    for(i in 1:length(mods)){
        modi <- mods[[i]]
        df <- cbind(modi$data, data.frame(resid(modi)))
        df$PA <- as.factor(df$PA)
        cols <- c('long', 'lat', 'resid.modi.')
        df <- merge(df, dat_for_geom, by=c('utm_east', 'utm_north'), all.x = TRUE)[cols]
        names(df)[names(df) == "resid.modi."] <- paste(species, ind, li_modnames[i], sep='_')
        
        dfresid <- merge(dfresid, df, by=c('long', 'lat'), all.x=TRUE)
    }
    
    return(dfresid)
}

getResidSpecies <- function(species){
    dfresidsr <- getResidSpeciesInd(species, 'sr')
    dfresidfr <- getResidSpeciesInd(species, 'fr')
    dfresidpd <- getResidSpeciesInd(species, 'pd')
    
    dfresidall <- merge(dfresidsr, dfresidfr, by=c('long', 'lat'))
    dfresidall <- merge(dfresidall, dfresidpd, by=c('long', 'lat'))
    
    return(dfresidall)
}

dfBirdResid <- getResidSpecies('bird')
write.csv(dfBirdResid, file.path(modelfolder, 'BirdsResidualsAllModels.csv'))

dfMammalResid <- getResidSpecies('mammal')
write.csv(dfMammalResid, file.path(modelfolder, 'MammalsResidualsAllModels.csv'))
