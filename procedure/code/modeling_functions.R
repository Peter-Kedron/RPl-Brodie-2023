## -------------------------------------------------------------------
## Script name: modeling_functions
## Purpose of script: A few functions for modeling.

## Usage
## In your script:
## source(path/to/modeling_functions.R)
## -------------------------------------------------------------------

# load packages and dependencies as of 2023-12-05
library(groundhog)
pkgs <- c("tidyverse", "cowplot", "here", "dagitty", "ggdag", "Hmisc", 
          "MatchIt", "modelsummary", "optmatch", "nlme")
groundhog.library(pkgs, "2023-12-05")

prepare_data <- function(taxon, sta_dir, conn_dir){
    # Load the data provided by Brodie et al. (2023)
    dat_brodie <- data.frame(read.csv(
        file.path(sta_dir, sprintf("training/%s_data_230326.csv", taxon)), 
        header = T))
    
    # Simplify the variable names of site identifier and geographic coordinates
    names(dat_brodie)[names(dat_brodie) == "site"] <- "station"
    names(dat_brodie)[names(dat_brodie) == "lat_wgs84"] <- "lat"
    names(dat_brodie)[names(dat_brodie) == "long_wgs84"] <- "long"
    
    # Assign stations the HDI value of its country
    # Reference: Human Development Report 2020: The Next Frontier—Human Development 
    # and the Anthropocene (United Nations Development Programme, 2020).
    # https://hdr.undp.org/data-center/human-development-index#/indicies/HDI
    # https://hdr.undp.org/sites/default/files/2021-22_HDR/HDR21-22_Statistical_Annex_HDI_Table.xlsx
    
    dat_HDI <- data.frame(
        country = c("Cambodia", "China", "Indonesia", "Laos", "Malaysia", 
                    "Singapore", "Thailand", "Vietnam", "Brunei"), 
        HDI = c(0.593, 0.768, 0.705, 0.607, 0.803, 0.939, 0.800, 
                0.703, 0.829))
    dat_brodie <- left_join(dat_brodie, dat_HDI, by = "country")
    
    # Create dataframe containing the subset of variable used in the analysis
    if (taxon == "bird"){
        dat <- dat_brodie %>% 
            select(station, country, PA, utm_east, utm_north, 
                   Hansen_recentloss,access_log10, HDI, dist_to_PA, 
                   PA_size_km2, rh_95_a0.pred, pavd_0_5.pred, 
                   pai_a0.pred, fhd_pai_1m_a0.pred, cover_a0.pred, 
                   agbd_a0.pred, asymptPD, maxFRic, SR.mean)
    } else {
        dat <- dat_brodie %>% 
            select(station, study_area, country, PA, utm_east, utm_north, 
                   Hansen_recentloss,access_log10, HDI, dist_to_PA, 
                   PA_size_km2, rh_95_a0.pred, pavd_0_5.pred, 
                   pai_a0.pred, fhd_pai_1m_a0.pred, cover_a0.pred, 
                   agbd_a0.pred, asymptPD, maxFRic, SR.mean)
    }
    
    # Add the connectivity variables for each station calculated with 
    # calc_conn_metrics.R
    dat_conn_metrics <- data.frame(read.csv(
        file.path(conn_dir, sprintf("conn_flux_%s_10_150.csv", taxon)), 
        header = T))
    dat <- left_join(dat, dat_conn_metrics, by = "station")
    
    # Scale subset of continuous variables in dat
    # Peter - Need to add the connectivity measures to the scaling list when they 
    # are introduced
    dat_scale <- subset(
        dat, select = c("utm_east", "utm_north", 
                        "HDI", "access_log10", 
                        "PA_size_km2", 
                        "dist_to_PA", 
                        "rh_95_a0.pred", 
                        "pavd_0_5.pred", 
                        "pai_a0.pred", 
                        "fhd_pai_1m_a0.pred", 
                        "cover_a0.pred", 
                        "agbd_a0.pred",
                        "awf_ptg")) %>% 
            scale(. , center = TRUE, scale = TRUE) %>% data.frame()
    
    # Append scaled variables to data with .z suffixes
    dat[paste0(names(dat_scale), '.z')] <- dat_scale
    
    # Rename scaled, predicted relative canopy height at 95% (rh_95_a0.pred.z) 
    # as forest_structure to match DAG 
    names(dat)[names(dat) == "rh_95_a0.pred.z"] <- "forest_structure"
    
    # Rename scaled, predicted plant area volume density between 0m and 5m 
    # (pavd_0_5.pred.z) as understory_density to match DAG
    names(dat)[names(dat) == "pavd_0_5.pred.z"] <- "understory_density"
    
    # Exclude stations that underwent recent forest loss as defined by 
    # Hansen et al. (2013)
    dat_clean <- subset(dat, Hansen_recentloss == 0)
    
    return(dat_clean)
}

modeling <- function(dat, 
                     independent_variable = "asymptPD",
                     med_dispersal_dist = 100, 
                     outliers = NULL, 
                     mod_type = "connec",
                     taxon = "bird"){
    dat <- subset(dat, med_dist == med_dispersal_dist)
    # Outliers:
    # stations to remove
    # mod_type: model type
    # ["brodie", "connec"]
    
    # Remove high-leverage outliers identified by Brodie et al. 
    ### Brodie et al. identified outlier using the hatvalue function. We should run 
    ### the analysis with their outlier set removed, but also run the hatvalues 
    ### analysis to identify the outliers for our particular specification thereby 
    ### mirroring their procedure.
    if (!is.null(outliers)){
        dat <- dat[!dat$station %in% outliers, ]
    }
    
    # Select variables for analysis and restrict to rows with complete values 
    # Peter - We need to add the eventual connectivity measures to this selection 
    # so they are in place for our extended analysis
    dat <- dat %>% rename(y = all_of(independent_variable))
    
    if ("study_area" %in% names(dat)){
        dat <- dat %>% 
            select(y, PA, study_area, country, utm_east, 
                   utm_north, utm_east.z, utm_north.z, forest_structure, 
                   access_log10.z, HDI.z, awf_ptg.z)
    } else {
        dat <- dat %>% 
            select(y, PA, country, utm_east, 
                   utm_north, utm_east.z, utm_north.z, forest_structure, 
                   access_log10.z, HDI.z, awf_ptg.z)
    }
    
    dat <- dat[complete.cases(dat), ]
    
    # Perform propensity score matching following the DAG developed in the 
    # structural causal modeling  and retrieve the matched dataset
    match_mod <- matchit(PA ~ utm_east.z + utm_north.z + forest_structure + 
                             access_log10.z + HDI.z,	
                         data = dat, method = "full", 
                         distance = "glm", link = "probit", replace = F)
    dat_matched <- match.data(match_mod)
    
    if (mod_type == "brodie"){
        # Run original Brodie linear mixed effects model with exponential spatial 
        # correlation structure for the residuals 
        if (taxon == "bird"){
            mod_efficacy <- lme(
                y ~ forest_structure + access_log10.z 
                + HDI.z + PA, random = list(~1 | country), 
                data = dat_matched, weights = ~I(1/weights), 
                correlation = corExp(form = ~utm_east + utm_north, 
                                     nugget = TRUE))
        } else if (taxon == "mammal"){
            mod_efficacy <- lme(
                y ~ forest_structure + access_log10.z 
                + HDI.z + PA, random = list(~1 | country, ~1 | study_area), 
                data = dat_matched, weights = ~I(1/weights), 
                correlation = corExp(form = ~utm_east + utm_north, 
                                     nugget = TRUE))
        } else {
            stop("No such taxon, mammal or bird.")
        }
        
    } else if (mod_type == "connec") {
        # Peter - we may want to introduce a spatial autocorrelation check of the model 
        # residuals
        
        # Run linear mixed effect model with the addition of connectivity moderator 
        # + conn + conn:PA
        if (taxon == "bird"){
            mod_efficacy <- lme(
                y ~ forest_structure + access_log10.z 
                + HDI.z + PA + awf_ptg.z + awf_ptg.z:PA, 
                random = list(~1 | country), data = dat_matched, 
                weights = ~I(1/weights), 
                correlation = corExp(form = ~utm_east + utm_north, 
                                     nugget = TRUE))
        } else if (taxon == "mammal"){
            mod_efficacy <- lme(
                y ~ forest_structure + access_log10.z 
                + HDI.z + PA + awf_ptg.z + awf_ptg.z:PA, 
                random = list(~1 | country, ~1 | study_area), 
                data = dat_matched, weights = ~I(1/weights), 
                correlation = corExp(form = ~utm_east + utm_north, 
                                     nugget = TRUE))
        } else {
            stop("No such taxon, mammal or bird.")
        }
    } else {
        stop("No such model type, brodie or connec.")
    }
    
    # Return the model
    return(mod_efficacy)
}
