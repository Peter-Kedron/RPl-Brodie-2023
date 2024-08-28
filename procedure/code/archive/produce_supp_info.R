# ====== General Information ======
# Task of the script: produce supplemental information
# Author: Wenxin Yang
# Date created: 01/19/2024
# Updated: 04/05/2024

# ======= load packages & set up dir =======
pkgs <- c('sf', 'dplyr', 'terra', 'ggplot2', 'units', 'cowplot', 'nngeo',
          'tidyverse', 'hrbrthemes', 'viridis', 'here', 'spdep', 'tmap',
          'data.table', 'stringr', 'broom', 'haven', 'extrafont', 'nlme',
          'broom', 'broom.mixed', 'showtext', 'ggthemes')


lapply(pkgs, library, character.only=TRUE)
sf_use_s2(FALSE) # deal with buffering odd
theme_set(theme_bw())
dir_pa <- '/Users/wenxinyang/Desktop/Dissertation/DATA'


load_object <- function(file) {
    tmp <- new.env()
    load(file = file, envir = tmp)
    tmp[[ls(tmp)[1]]]
}

modelfolder <- '/Users/wenxinyang/Desktop/Dissertation/DATA/BrodieData/model_data'
modelfolder <- 'C:/Users/peterkedron/Desktop/Brodie_models'
# ======= prep ======

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

# cowplot::plot_grid(plotlist = list)

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
                          here("data/raw/public/training/bird_data_corrected_240122.csv"), 
                          header = T))
dat_brodie_mammal <- data.frame(read.csv(
    here("data/raw/public/training/mammal_data_corrected_240122.csv"),
    header = T))
# Simplify the variable names of site identifier and geographic coordinates
names(dat_brodie)[names(dat_brodie) == "site"] <- "station"
names(dat_brodie)[names(dat_brodie) == "lat_wgs84"] <- "lat"
names(dat_brodie)[names(dat_brodie) == "long_wgs84"] <- "long"

names(dat_brodie_mammal)[names(dat_brodie_mammal) == "lat_wgs84"] <- "lat"
names(dat_brodie_mammal)[names(dat_brodie_mammal) == "long_wgs84"] <- "long"


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
dat_geom_mammal <- dat_brodie_mammal %>% select(long, lat, utm_east, utm_north)
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
#write.csv(dfBirdResid, file.path(modelfolder, 'BirdsResidualsAllModels.csv'))
#dfBirdResid <- read.csv(file.path(modelfolder, 'BirdsResidualsAllModels.csv'))

dfMammalResid <- getResidSpecies('mammal')
#write.csv(dfMammalResid, file.path(modelfolder, 'MammalsResidualsAllModels.csv'))
#dfMammalResid <- read.csv(file.path(modelfolder, 'MammalsresidualsAllModels.csv'))

# ======= 3) Statistics summary of residuals in relation to dist to PAs ======
colsTask3 <- c('long', 'lat', 'dist_to_PA', 'PA')

prepDfTask3 <- function(init_df, resid_df){
    dftask3 <- init_df %>% select(long, lat, dist_to_PA, PA)
    dftask3 <- merge(resid_df, dftask3, by=c('long', 'lat'), all.y=TRUE)
    
    dfSumDisp <- data.frame(matrix(ncol=7, nrow=0))
    
    for(c in 1:ncol(dftask3)){
        col <- colnames(dftask3)[c]
        if(! col %in% colsTask3){
            dftmp <- dftask3[append(colsTask3, col)]
            colnames(dftmp) <- append(colsTask3, 'resid')
            dftmp$model <- col
            
            strs <- strsplit(col, '_')[[1]]
            ind <- strs[2]
            if(length(strs)>3){
                dftmp$dispersaldist <- as.integer(strs[4])
            } else{
                dftmp$dispersaldist <- 0
            }
            dfSumDisp <- rbind(dfSumDisp, dftmp)
        }
    }
    dfSumDisp$ifWithin <- dfSumDisp$dist_to_PA < dfSumDisp$dispersaldist
    return(dfSumDisp)
}

dfBirdTask3 <- prepDfTask3(dat_brodie, dfBirdResid)
dfBirdTask3.outPA <- dfBirdTask3[dfBirdTask3$PA == 0, ]
dfBirdTask3.outPA$category <- paste(dfBirdTask3.outPA$model, dfBirdTask3.outPA$ifWithin, 
                                    sep='_')
a <- tapply(dfBirdTask3.outPA$resid,
            dfBirdTask3.outPA$category,
            summary)

dfMammalTask3 <- prepDfTask3(dat_brodie_mammal, dfMammalResid)
dfBirdTask3.outPA <- dfBirdTask3[dfBirdTask3$PA == 0, ]


# ====== local Moran's I (in progress) ======
# a little suspicious --> I used GeoDa instead.
# inverse distance weighting
geom_birdresid <- vect(dfBirdResid, geom=c('long', 'lat'), crs='epsg:4326')
df.dist <- as.matrix(dist(cbind(dfBirdResid$long, dfBirdResid$lat)))
df.dist.inv <- 1/df.dist
diag(df.dist.inv) <- 0

geom_birdresid$brodsrlcMoranI <- autocor(geom_birdresid$bird_sr_brodie, df.dist.inv, "locmor")
plot(geom_birdresid, "brodsrlcMoranI")

geom_birdresid$sr10lcMoranI <- autocor(geom_birdresid$bird_sr_connec_10, df.dist.inv, "locmor")
plot(geom_birdresid, "sr10lcMoranI")

# ====== 4) tie fighter pie plot =====
# PD for birds and mammals: PA estimates and awf_ptg.z estimates
loadfonts(device=c("all"))
import_econ_sans() # you will have to install the font before calling it
title_theme <- theme(
    plot.title = element_text(
        family = "Econ Sans Cnd",
        face = "bold",
        size = 12
    ),
    plot.subtitle = element_text(
        family = "Econ Sans Cnd",
        size = 10,
    ),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank()
)

birds_pd_models = load_object(file.path(modelfolder, 'bird_pd_models.rda'))
mammals_pd_models = load_object(file.path(modelfolder, 'mammal_pd_models.rda'))

getEachMod <- function(mod, vars, modname){
    modinfo <- broom.mixed::tidy(mod, effects='fixed',conf.int = TRUE)
    modinfo <- modinfo[modinfo$term %in% vars,]
    modinfo <- modinfo[c('conf.low', 'estimate', 'conf.high', 'term')]
    colnames(modinfo) <- c('lower', 'est.', 'upper', 'var')
    modinfo$model <- modname
    modinfo$dist <- as.integer(strsplit(modname, '_')[[1]][2])
    
    return(modinfo)
}

getAllMod <- function(species, ind, vars){
    mods <- load_object(file.path(modelfolder, paste(species, ind, 'models.rda', sep='_')))
    li_modnames <- names(mods)
    dfplot <- data.frame(matrix(ncol=6, nrow=0))
    
    for(i in 2:length(li_modnames)){ # removing brodie model
        mod <- mods[[i]]
        modname <- li_modnames[i]
        dfplot <- rbind(dfplot, getEachMod(mod, vars, modname))
        
    }
    rownames(dfplot) <- seq(1:nrow(dfplot))
    if(ind=='pd'){
        ind == 'Phylogenetic Diversity'
    }
    dfplot$ind <- toupper(ind)
    return(dfplot)
    
}


getFightPiePlot <- function(species, inds, vars, varnames, 
                            title_string, subtitle_string){
    dftmp <- data.frame(matrix(ncol=6, nrow=0))
    for (i in inds){
        dftmp <- rbind(dftmp, getAllMod(species, i, vars))
    }
    
    dfvarnames <- data.frame(vars, varnames)
    colnames(dfvarnames) <- c('var', 'varname')
    dftmp <- merge(dftmp, dfvarnames, by='var')
    
    dftmp %>% ggplot(aes(x=dist, y=est.))+
        # geom_point(position=position_dodge(width=2))+
        geom_pointrange(
            size=0.2,
            aes(ymin=lower, ymax=upper),
                      position=position_dodge(width=2), width=1.5)+
        scale_x_continuous(breaks=seq(10, 150, by=10))+
        # ggtitle(paste('Plot for', species, ind, sep=' '))+
        xlab('Dispersal distance (km)') +
        ylab('Eestimated effect')+
        labs(
            title=title_string,
            subtitle=subtitle_string
        )+
        # geom_hline(yintercept=0.38, linetype='dashed', color = "red")+
        # ylab(paste(var, 'effect', sep=' '))+
        # facet_wrap(~ind+varname, scales = 'free_y')+
        theme_tufte()+
        theme(
            #strip.text.x = element_text(angle = 0, vjust = 0),
              panel.border = element_rect(colour='white', fill=NA))+
        title_theme
}



ttstr <- 'Replication Outcome'
subttstr <- 'Connectivity moderates the effect of PA on Phylogenetic Diversity'
getFightPiePlot('bird', c('pd'), c('PA:awf_ptg.z'), c('Interaction'),
                ttstr, subttstr)

li_vars <- c('PA', 'awf_ptg.z', 'PA:awf_ptg.z')
li_inds <- c('sr', 'fr', 'pd')
li_varnames <- c('PA', 'Connectivity', 'Interaction')


# ====== 5) dot whisker plots for beta estimates for the three bird models =====
# 04/05/2024
# use the broom.mixed package to tidy model results up
# Note: the R package is computing confidence interval using a different way from estimate+2*std.err but they are close
# for consistency I'm using estimate+2*std.err for both
getConfInt <- function(df){
    df$conf.low = df$estimate-2*df$std.error
    df$conf.high = df$estimate+2*df$std.error
    
    return(df)
}

cleanModOutput <- function(modfolder, modname, biodivind, rawtab){
    # read in the raw model
    rawmod <- load_object(file.path(modfolder, modname))
    # clean up model results
    df <- broom.mixed::tidy(rawmod, conf.int=FALSE)# if all results were saved as models, we can set this to TRUE and get the error bars automatically
    
    # further clean up to match with Brodie's reported table
    group_rpr <- paste0("rpr-", biodivind, sep='')
    df$term <- c('Intercept', 'Forest Canopy Height', 'Site Accessibility',
                 'HDI', 'PA', 'sd_int', 'sd_obs')
    df <- df[df$effect=='fixed', ]
    df <- getConfInt(df)
    df$df <- NULL
    df$statistic <- NULL
    df$p.value <- NULL
    df$group <- "Reproduction"
    
    # merging brodie's reported table
    group_raw <- paste0("raw-", biodivind, sep='')
    rawtabind <- rawtab[rawtab$group==group_raw, ]
    rawtabind$group <- "Brodie"
    tabind <- getConfInt(rawtabind)
    
    df <- rbind(df, tabind)
    df$ind <- toupper(biodivind)
    
    return(df)
}

brodietable <- read.csv('data/derived/public/Brodietable.csv')
bird_pd_fn <- cleanModOutput(modelfolder, 'bird_pd_model_orig.rda', 'pd', brodietable)
bird_sr_fn <- cleanModOutput(modelfolder, 'bird_sr_model_orig.rda', 'sr', brodietable)
bird_fr_fn <- cleanModOutput(modelfolder, 'bird_fr_model_orig.rda', 'fr', brodietable)

allbird <- rbind(bird_sr_fn, bird_fr_fn, bird_pd_fn)

unique(allbird$term)
termname <- c("PA", "Forest Canopy Height", "Site Accessibility", "HDI", "Intercept")
termorder <- c(1, 2, 3, 4, 5)
termdf <- data.frame(termname, termorder)
colnames(termdf) <- c("term", "order")
allbird <- merge(allbird, termdf, by=c("term"))

inddf <- data.frame(c("SR", "PD", "FR"), c("Species Richness", "Phylogenetic Diversity", "Functional Richness"))
colnames(inddf) <- c("ind", "indname")
allbird <- merge(allbird, inddf, by=c("ind"))

allbird %>%
    # reorder the coefficients so that the largest is at the top of the plot
    # PJK - Wenxin please reorder the variables in the following order top to bottom 
    # - PA, Forest, Access, HDI, Intercept. We want the PA on the top because it is the main effect
    mutate(term = fct_reorder(term, order, .desc = TRUE)) %>%
    ggplot(aes(estimate, term, col=group)) +
    # add in a dotted line at zero
    geom_vline(xintercept = 0, lty = 2) +
    # geom_point() +
    scale_colour_manual(values=c("black", "red"))+
    geom_pointrange(size=0.2, aes(xmin = conf.low, xmax = conf.high), 
                    position = position_dodge(width = .1)) +
    # errorbar has ends
    labs(
        x = "Estimate of effect",
        y = NULL, 
        title = "Computational Reproduction",
        subtitle = "Matching effect estimates observed across studies",
        colour = "Group",
    )+
    facet_wrap(~indname, scales = "free_x")+# scales='free_x' allows xlim to differ by subplot
    theme_tufte()+
    theme(strip.text.x = element_text(angle = 0, hjust = 0),
          panel.border = element_rect(colour='white', fill=NA),
          legend.position = "none")+
    title_theme +
    guides(fill = FALSE) # I believe this will remove the legend on the right

ggsave(here("results/figures/comp_repo_compare.png"), width = 8, height = 6.5, bg = "white")

# adding a geom_rangeframe() is weird

# Figure of Connectivity Moderates PA Efficacy
## Note: use 100 km dispersal distance, a bit hardcoded
## Load function, or you can function prepare_data load manually
source(here("procedure/code/modeling_functions.R"))

## Load models
load(here("data/derived/public/bird_pd_models.rda"))
intercept <- bird_pd_models$connec_100$coefficients$fixed["(Intercept)"]
B_conn <- bird_pd_models$connec_100$coefficients$fixed["awf_ptg.z"]
B_connPA <- B_conn + bird_pd_models$connec_100$coefficients$fixed["PA:awf_ptg.z"]

## Reproduce the same data for the target model
dat_clean <- prepare_data("bird", here("data/raw/public"),
                          here("data/derived/public"))

# Recreate the matched dataset (we actually could simply grab it from model result)
med_dispersal_dist <- 100
outliers <- c("L921125", "L2422371", "L4331944", "L13465594")
dat_clean <- subset(dat_clean, med_dist == med_dispersal_dist)
dat_clean <- dat_clean[!dat_clean$station %in% outliers, ]
dat_clean <- dat_clean %>% 
    select(asymptPD, PA, country, utm_east, 
           utm_north, utm_east.z, utm_north.z, forest_structure, 
           access_log10.z, HDI.z, awf_ptg.z)
dat_clean <- dat_clean[complete.cases(dat_clean), ]
match_mod <- matchit(PA ~ utm_east.z + utm_north.z + forest_structure + 
                         access_log10.z + HDI.z,	
                     data = dat_clean, method = "full", 
                     distance = "glm", link = "probit", replace = F)
dat_matched <- match.data(match_mod)

dat_matched <- dat_matched %>%  
    arrange(PA) %>% # make sure the PA points show front
    mutate(PA = factor(PA, levels = c(0, 1), labels = c(0, 1)))

# Making a Economist inspired theme for the plot title
title_theme <- theme(
    plot.title = element_text(
        family = "Econ Sans Cnd", 
        face = "bold",
        size = 12
    ),
    plot.subtitle = element_text(
        family = "Econ Sans Cnd",
        size = 10,
    ),
    legend.position="none"
)

ggplot(data = dat_matched) +
    geom_point(aes(x = awf_ptg.z, y = asymptPD, 
                   color = PA, size = weights)) +
    scale_color_manual(values = c("grey80", "#009E73")) +
    geom_abline(slope = B_conn, intercept = intercept, 
                color = "grey80", lwd = 0.8) +
    geom_abline(slope = B_connPA, intercept = intercept, 
                color = "#009E73", lwd = 0.8) +
    geom_rangeframe() +
    ylim(1.5, 6) +
    labs( x = "Connectivity",
          y = "Phylogenetic Diversity",
          title = "Connectivity Moderates PA Efficacy",
          subtitle = "On phylogenetic diversity, at a dispersal distance of 100km") +
    theme_tufte() + title_theme

ggsave(here("results/figures/PD_conn2.png"), width = 8, height = 5, bg = "white")

# Figure of comparison between new model and brodie model
conn_mod <- coef(summary(bird_pd_models$connec_100)) %>% 
    as.data.frame() %>% 
    mutate(item = rownames(.)) %>% 
    select(item, Value, Std.Error) %>% 
    mutate(group = "Conn")

rep_mod <- brodietable %>% filter(group == "raw-pd") %>% 
    rename(item = term, Value = estimate, Std.Error = std.error) %>% 
    select(item, Value, Std.Error) %>% 
    mutate(group = "Rep")
rep_mod$item <- c("(Intercept)", "forest_structure", "access_log10.z", "HDI.z", "PA")

items <- c("PA", "awf_ptg.z", "PA:awf_ptg.z", "forest_structure", 
           "access_log10.z", "HDI.z", "(Intercept)")
labels <- c("Protected Area (PA)", "Connectivity", "PA | Connectivity", "Forest Canopy Height",
            "Site Accessibility", "HDI", "Intercept")
dat <- rbind(conn_mod, rep_mod) %>% 
    mutate(item = factor(
        item, levels = rev(items), labels = rev(labels)))


ggplot(dat, aes(Value, item, col=group)) +
    # add in a dotted line at zero
    geom_vline(xintercept = 0, lty = 2) +
    # geom_point() +
    scale_colour_manual(values=c("#009E73", "black"))+
    geom_pointrange(size=0.2, aes(xmin = Value - 2 * Std.Error, 
                                  xmax = Value + 2 * Std.Error), 
                    position = position_dodge(width = .1)) +
    # errorbar has ends
    labs(
        x = "Estimate of effect",
        y = NULL, 
    )+
    theme_tufte()+
    theme(strip.text.x = element_text(angle = 0, hjust = 0),
          panel.border = element_rect(colour='white', fill=NA),
          legend.position = "none")+
    title_theme

ggsave(here("results/figures/Conn100_rep_compare.png"), width = 3, height = 6.5, bg = "white")

# Study area figure
library(rnaturalearth)
pas <- st_read("data/derived/public/clean_pas.geojson") %>% rmapshaper::ms_simplify()
pts <- st_read("data/derived/public/pts_bird.geojson") %>% 
    left_join(dat_clean %>% select(station, awf_ptg.z), by = "station")

bry <- ne_countries(scale = 10, continent = "asia", returnclass = "sf") %>% select(name) %>% 
    st_transform(st_crs(pas)) %>% 
    st_intersection(st_as_sfc(st_bbox(pas)))

highlight <- st_as_sf(data.frame(x = 107.5, y = 12.5), 
                          coords = c("x", "y"),
                          crs = 4326) %>% 
    st_transform(st_crs(pas)) %>% st_buffer(2 * 1e5)

pas_hl <- st_intersection(pas, highlight)
pts_hl <- st_intersection(pts, highlight)

g_highlight <- ggplotGrob(ggplot() + 
    geom_sf(data = highlight, color = "#009E73", fill = "cornsilk", linewidth = 0.6) +
    geom_sf(data = pas_hl, color = "grey20", fill = "transparent", linewidth = 0.5) +
    geom_sf(data = pts_hl, aes(fill = awf_ptg.z), size = 2, shape = 21, color = "white") +
    ggsci::scale_fill_material(name = 'Connectivity', "deep-orange", limit = c(-1, 3)) +
    coord_sf() + theme_map() + theme(legend.position = "none"))

ggplot() + 
    geom_sf(data = bry, color = "cornsilk", fill = "cornsilk") +
    geom_sf(data = pas, color = "grey20", fill = "transparent", linewidth = 0.5) +
    geom_sf(data = pts, aes(fill = awf_ptg.z), size = 1, 
            shape = 21, color = "white", alpha = 0.3) +
    geom_sf(data = highlight, color = '#009E73', 
            fill = "transparent", linewidth = 0.6) +
    xlim(st_bbox(pas)[1], st_bbox(pas)[3]) +
    ylim(st_bbox(pas)[2], st_bbox(pas)[4]) +
    ggsci::scale_fill_material(name = 'Connectivity', "deep-orange", limit = c(-1, 3)) +
    annotation_custom(grob = g_highlight,
                      xmin = 2173181, xmax = 3602136,
                      ymin = 906978, ymax = 2704781) +
    coord_sf() + theme_linedraw() +
    theme(legend.position = "none",
          text = element_text(size = 12))

ggsave(here("results/figures/study_area.png"), width = 6.5, height = 6.5, bg = "white")
