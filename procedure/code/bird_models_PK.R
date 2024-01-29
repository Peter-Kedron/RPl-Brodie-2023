# Adapted from Brodie et al. (2023) 
# Code downloaded from: https://doi.org/10.5281/zenodo.7796347
# Data downloaded from: https://figshare.com/articles/dataset/Data/22527295
#
# Modified by Lei Song (lsong@ucsb.edu), Peter Kedron (peterkedron@ucsb.edu)
# Last updated: Dec 13, 2023
#
# Comments on original code and modifications:----------------------------------
# 
# Issue 1: The author provided data file did not include several variables used 
#          in the analysis. We corrected the data file following the description
#          of data construction provided in the published article. 
#
# Issue 2: While measures of the response variables PD, FR, and SR were provided 
#          in the data file, how these variables were constructed is not 
#          described. Moreover, the original data used to produce these 
#          variables was not available, which limited their evaluation. 
#
# Issue 3: Procedures to construct the predictor variables forest structure and
#          understory density were not provided. Only the final values were 
#          provided. Original GEDI can be tracked down online, but cannot be 
#          used to construct the predictor variables.
# 
# Issue 4: The authors did not base the standardization of their predictors in 
#          the spillover analyses on the subsets of data used in those analyses. 
#          Rather the standardization from the full data set was carried into 
#          these analyses.


# Load Packages and Set Working Directory --------------------------------------

# load packages and dependencies as of 2023-12-05
library(groundhog)
pkgs <- c("tidyverse", "cowplot", "here", "dagitty", "ggdag", "Hmisc", 
          "MatchIt", "modelsummary", "optmatch", "nlme")
groundhog.library(pkgs, "2023-12-05")

# Set folder hierarchy to data folder
here('data/')

# ------------------------------------------------------------------------------
# Structural Causal Modelling  
# ------------------------------------------------------------------------------

# Create the directed acyclic graph (DAG) of potential causal pathways with the 
# addition of connectivity as a moderator of the PA effect. Assign PA as 
# the exposure and diversity as the outcome.

dagBrodie <- dagitty("dag {
  PA -> Diversity
  Connectivity -> Diversity
  ForestStructure -> PA
  SiteAccessibility -> PA
  Bioclimate -> ForestStructure -> PA -> Diversity
  Bioclimate -> UnderstoryDensity -> Diversity
  Bioclimate -> Diversity
  Elevation -> Bioclimate
  Elevation -> ForestStructure
  Elevation -> UnderstoryDensity
  Elevation -> SiteAccessibility
  Topography -> UnderstoryDensity
  Topography -> ForestStructure
  Topography -> SiteAccessibility
  Topography -> Diversity
  HDI -> Diversity
  PA [exposure]
  Diversity [outcome]
               }"
                     )

# Organize the data into a visual hierarchy
coordinates( dagBrodie ) <-  list(x = c(Diversity = 3, 
                                        UnderstoryDensity = 1, 
                                        ForestStructure = 2, 
                                        PA = 3, 
                                        SiteAccessibility = 4, 
                                        HDI = 5, 
                                        Bioclimate = 2, 
                                        Elevation = 3, 
                                        Topography = 4,
                                        Connectivity = 6),
                                  y = c(Diversity = 3, 
                                        UnderstoryDensity = 2, 
                                        ForestStructure = 2, 
                                        PA = 2, 
                                        SiteAccessibility = 2, 
                                        HDI = 2, 
                                        Bioclimate = 1, 
                                        Elevation = 1, 
                                        Topography = 1, 
                                        Connectivity = 2))

# Plot the DAG to confirm the visual structure and exposure
ggdag_status(dagBrodie) + theme_dag()

# Identify the set of adjustment variables needed to identify 
# the effect of PA on biodiversity and Test whether connectivity fulfills the 
# adjustment criterion
adjustmentSets(dagBrodie)
isAdjustmentSet(dagBrodie, c("Connectivity"))

# Plot the alternative adjustment sets
ggdag_adjustment_set(dagBrodie, 
                     node_size = 20, 
                     text_col = "black"
                     ) + theme(legend.position = "bottom")


# ------------------------------------------------------------------------------
# Process Data and Construct Variables
# ------------------------------------------------------------------------------
# Deviation 1: The data file shared by the original authors did not include the
# human development index (HDI) measure. Following the in-text description we 
# assigned each station the HDI of its country.

# Clarification 1: Brodie et al. identified two GEDI-based measures of forest 
# structure as the best fit for their diversity data. We follow their procedure, 
# but clarify the mapping of measure to concepts by matching variable names to 
# those used in the structural causal modelling

# Extension 1: We introduce 2 measures of station to PA connectivity, calculated
# following the procedures defined in calc_conn_metrics.R.


# Clear the workspace
rm(list = ls())

# Load the data provided by Brodie et al. (2023)
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
dat <- dat_brodie %>% select(station, country, PA, utm_east, utm_north, 
                             Hansen_recentloss,access_log10, HDI, dist_to_PA, 
                             PA_size_km2, rh_95_a0.pred, pavd_0_5.pred, 
                             pai_a0.pred, fhd_pai_1m_a0.pred, cover_a0.pred, 
                             agbd_a0.pred,asymptPD, maxFRic, SR.mean)

# Add the connectivity variables for each station calculated with 
# calc_conn_metrics.R
dat_conn_metrics <- data.frame(read.csv(
                                here("data/derived/public/conn_flux_bird_10_150.csv"), 
                                header = T))
dat <- left_join(dat, dat_conn_metrics, by = "station")

# Scale subset of continuous variables in dat
# Peter - Need to add the connectivity measures to the scaling list when they 
# are introduced
dat_scale <- data.frame(scale(subset(dat, select = c("utm_east", "utm_north", 
                                                     "HDI", "access_log10", 
                                                     "PA_size_km2", 
                                                     "dist_to_PA", 
                                                     "rh_95_a0.pred", 
                                                     "pavd_0_5.pred", 
                                                     "pai_a0.pred", 
                                                     "fhd_pai_1m_a0.pred", 
                                                     "cover_a0.pred", 
                                                     "agbd_a0.pred",
                                                     "awf_rst_ptp2"), 
                                     ), 
                              center = TRUE, scale = TRUE))

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


# ------------------------------------------------------------------------------
# Linear Mixed-Effects Modelling and Propensity Score Matching
# ------------------------------------------------------------------------------

# ----------------- Analysis of phylogenetic diversity (PD) --------------------

# Linear Mixed Effects Model of PA Efficacy ------------------------------------

# Extension 1: We introduce an interaction term to assess the moderating effect
# of PA connectivity on the effect of PA status on biodiversity

# Modification: We may attempt to match without the lat long coordinates 
# included in the probit model.

# Replicate data frame for PD sub-analysis focus analysis on 100m dispersal 
# distance
dat_PD_efficacy <- subset(dat_clean, med_dist == 80)

# Remove high-leverage outliers identified by Brodie et al. 
### Brodie et al. identified outlier using the hatvalue function. We should run 
### the analysis with their outlier set removed, but also run the hatvalues 
### analysis to identify the outliers for our particular specification thereby 
### mirroring their procedure.

PD_efficacy_outliers <- c("L2422371", "L3776738", "L2521761", "L6127181", 
                          "L3865754")
dat_PD_efficacy <- dat_PD_efficacy[! dat_PD_efficacy$station %in% 
                                       PD_efficacy_outliers, ]

# Select variables for analysis and restrict to rows with complete values 
# Peter - We need to add the eventual connectivity measures to this selection 
# so they are in place for our extended analysis

dat_PD_efficacy <- dat_PD_efficacy %>% select(asymptPD, PA, country, utm_east, 
                          utm_north, utm_east.z, utm_north.z, forest_structure, 
                          access_log10.z, HDI.z, awf_rst_ptp2.z)
dat_PD_efficacy <- dat_PD_efficacy[complete.cases(dat_PD_efficacy), ]

# Perform propensity score matching following the DAG developed in the 
# structural causal modeling  and retrieve the matched dataset
match_PD <- matchit(PA ~ utm_east.z + utm_north.z + forest_structure + 
                        access_log10.z + HDI.z,	
                    data = dat_PD_efficacy, method = "full", 
                    distance = "glm", link = "probit", replace = F)
dat_matched_PD <- match.data(match_PD)

# Run original Brodie linear mixed effects model with exponential spatial 
# correlation structure for the residuals 
mod_PD_efficacy <- lme(asymptPD ~ forest_structure + access_log10.z 
                       + HDI.z + PA, random = list(~1 | country), 
                       data = dat_matched_PD, weights = ~I(1/weights), 
                       correlation = corExp(form = ~utm_east + utm_north, 
                                            nugget = TRUE))
summary(mod_PD_efficacy) 

# Peter - we may want to introduce a spatial autocorrelation check of the model 
# residuals

# Run linear mixed effect model with the addition of connectivity moderator 
# + conn + conn:PA
mod_CN_efficacy <- lme(asymptPD ~ forest_structure + access_log10.z 
                       + HDI.z + PA + awf_rst_ptp2.z + awf_rst_ptp2.z:PA, 
                       random = list(~1 | country), data = dat_matched_PD, 
                       weights = ~I(1/weights), 
                       correlation = corExp(form = ~utm_east + utm_north, 
                       nugget = TRUE))
summary(mod_CN_efficacy)

# Summarize the two model outputs in a table
msummary(list(mod_PD_efficacy, mod_CN_efficacy))


# Linear Mixed Effects Model of PA Size Spillovers -----------------------------

# Modification: We will re-standardize the predictor variables after subset for
# each of the two spillover analyses, rather than just carrying through the 
# universal standardization.

# Subset dataset to stations outside of PA
dat_PA_size <- subset(dat_clean, PA==0)

# Label nearest PA that are above a 500km2 size threshold as large PA
PA_size_threshold.z <- (500 - mean(dat_clean$PA_size_km2)) / 
                        sd(dat_clean$PA_size_km2)
dat_PA_size$BigPA <- ifelse(dat_PA_size$PA_size_km2.z < PA_size_threshold.z, 
                            0, 1)

# Remove high-leverage outliers identified by Brodie et al. 
dat_PD_size_outliers <- c("L1084299", "L4225511", "L3846512", "L2129865", 
                          "L3267752")
dat_PD_spill <- dat_PA_size[! dat_PA_size$station %in% 
                             dat_PD_size_outliers, ]

# Select variables for analysis and restrict to rows with complete values 
# Peter - We need to add the eventual connectivity measures to this selection 
# so they are in place for our extended analysis
dat_PD_spill <- dat_PD_spill %>% select(asymptPD, BigPA, dist_to_PA.z, 
                                      country, utm_east, 
                                      utm_north, utm_east.z, utm_north.z, 
                                      forest_structure, access_log10.z, HDI.z)
dat_PD_spill <- dat_PD_spill[complete.cases(dat_PD_spill), ]

# Perform propensity score matching again following the DAG
match_size <- matchit(BigPA ~ utm_north.z + utm_east.z + forest_structure + 
                              access_log10.z + HDI.z + dist_to_PA.z, 
                      data = dat_PD_spill, method = "full", distance = "glm", 
                      link = "probit", replace = F)
data_matched_size <- match.data(match_size)

# Run original Brodie linear mixed effects model with exponential spatial 
# correlation structure for the residuals 
mod_PD_size <- lme(asymptPD ~ forest_structure + access_log10.z + HDI.z + 
                              dist_to_PA.z + BigPA, 
                   random = list(~1 | country), data = data_matched_size, 
                   weights = ~I(1/weights), 
                   correlation = corExp(form = ~utm_east + utm_north, 
                                        nugget = TRUE))
summary(mod_PD_size)

# Peter - repeat the model above with the connectivity metrics included if we 
# think this is an analysis we wish to pursue


# Linear Mixed Effects Model of PA Distance Spillovers -----------------------------

# Modification: We will re-standardize the predictor variables after subset for
# each of the two spillover analyses, rather than just carrying through the 
# universal standardization.

# Label nearest PA that are > 2km  size threshold as large PA
dat_PA_size$CloseToPA <- ifelse(dat_PA_size$dist_to_PA > 2, 0, 1)

# Select variables for analysis and restrict to rows with complete values 

# Peter - We need to add the eventual connectivity measures to this selection 
# so they are in place for our extended analysis.
dat_PD_spill <- dat_PA_size[! dat_PA_size$station %in% 
                                dat_PD_size_outliers, ]

dat_PD_spill <- dat_PD_spill %>% select(asymptPD, CloseToPA, PA_size_km2.z, 
                                       country, utm_east, utm_north, 
                                       utm_east.z, utm_north.z, 
                                       forest_structure, access_log10.z, HDI.z)
dat_PD_spill <- dat_PD_spill[complete.cases(dat_PD_spill), ]

# Perform propensity score matching again following the DAG
match_dist <- matchit(CloseToPA ~ utm_north.z + utm_east.z + forest_structure 
                                  + access_log10.z + HDI.z + PA_size_km2.z, 
                      data = dat_PD_spill, method = "full", distance = "glm", 
                      link = "probit", replace = F)
dat_matched_dist <- match.data(match_dist)

# Run original Brodie linear mixed effects model with exponential spatial 
# correlation structure for the residuals 
mod_PD_dist <- lme(asymptPD ~ forest_structure + access_log10.z + HDI.z + 
                              PA_size_km2.z + CloseToPA, 
                   random = list(~1 | country), data = dat_matched_dist, 
                   weights = ~I(1/weights), 
                   correlation = corExp(form = ~utm_east + utm_north, 
                                        nugget = TRUE))
summary(mod_PD_dist)

# Peter - repeat the model above with the connectivity metrics included if we 
# think this is an analysis we wish to pursue

# ----------------- Analysis of Functional Richness (FR) -----------------------

# Linear Mixed Effects Model of PA Efficacy ------------------------------------

# Extension 1: We introduce an interaction term to assess the moderating effect
# of PA connectivity on the effect of PA status on biodiversity

# Modification: We may attempt to match without the lat long coordinates 
# included in the probit model.

# Start over, replicate data frame for FR sub-analysis
dat_FR_efficacy <- subset(dat_clean, med_dist == 100)

# Remove high-leverage outliers identified by Brodie et al. 
### Brodie et al. identified outlier using the hatvalue function. We should run 
### the analysis with their outlier set removed, but also run the hatvalues 
### analysis to identify the outliers for our particular specification thereby 
### mirroring their procedure.

FR_efficacy_outliers <- c("L921125", "L2422371", "L4331944", "L13465594")
dat_FR_efficacy <- dat_FR_efficacy[! dat_FR_efficacy$station %in% 
                                       FR_efficacy_outliers, ]

# Select variables for analysis and restrict to rows with complete values 
# Peter - We need to add the eventual connectivity measures to this selection 
# so they are in place for our extended analysis

# No variable named maxFR in the csv, only maxFRic, not fully sure if they are the same
dat_FR_efficacy <- dat_FR_efficacy %>% 
    select(maxFRic, PA, country, utm_east, 
           utm_north, utm_east.z, utm_north.z, forest_structure, 
           access_log10.z, HDI.z, awf_rst_ptp2.z)
dat_FR_efficacy <- dat_FR_efficacy[complete.cases(dat_FR_efficacy), ]

# Perform propensity score matching following the DAG developed in the 
# structural causal modeling  and retrieve the matched dataset
match_FR <- matchit(PA ~ utm_east.z + utm_north.z + forest_structure + 
                        access_log10.z + HDI.z,	
                    data = dat_FR_efficacy, method = "full", 
                    distance = "glm", link = "probit", replace = F)
dat_matched_FR <- match.data(match_FR)

# Run original Brodie linear mixed effects model with exponential spatial 
# correlation structure for the residuals 
mod_FR_efficacy <- lme(maxFRic ~ forest_structure + access_log10.z 
                       + HDI.z + PA, random = list(~1 | country), 
                       data = dat_matched_FR, weights = ~I(1/weights), 
                       correlation = corExp(form = ~utm_east + utm_north, 
                                            nugget = TRUE))
summary(mod_FR_efficacy)

# Peter - we may want to introduce a spatial autocorrelation check of the model 
# residuals

# Run linear mixed effect model with the addition of connectivity moderator 
# + conn + conn:PA
mod_FR_CN_efficacy <- lme(maxFRic ~ forest_structure + access_log10.z 
                       + HDI.z + PA + awf_rst_ptp2.z + awf_rst_ptp2.z:PA, 
                       random = list(~1 | country), data = dat_matched_FR, 
                       weights = ~I(1/weights), 
                       correlation = corExp(form = ~utm_east + utm_north, 
                                            nugget = TRUE))
summary(mod_FR_CN_efficacy)

# Summarize the two model outputs in a table
msummary(list(mod_FR_efficacy, mod_FR_CN_efficacy))


# Linear Mixed Effects Model of PA Size Spillovers -----------------------------

# Modification: We will re-standardize the predictor variables after subset for
# each of the two spillover analyses, rather than just carrying through the 
# universal standardization.

# Subset dataset to stations outside of PA
dat_PA_size <- subset(dat_clean, PA==0)

# Label nearest PA that are above a 500km2 size threshold as large PA
PA_size_threshold.z <- (500 - mean(dat_clean$PA_size_km2)) / 
    sd(dat_clean$PA_size_km2)
dat_PA_size$BigPA <- ifelse(dat_PA_size$PA_size_km2.z < PA_size_threshold.z, 
                            0, 1)

# Remove high-leverage outliers identified by Brodie et al. 
dat_FR_size_outliers <- c("L4225511", "L5969878", "L3267752", "L4331944", 
                          "L13465594", "L1084299")
dat_FR_spill <- dat_PA_size[! dat_PA_size$station %in% 
                                dat_FR_size_outliers, ]

# Select variables for analysis and restrict to rows with complete values 
# Peter - We need to add the eventual connectivity measures to this selection 
# so they are in place for our extended analysis
dat_FR_spill <- dat_FR_spill %>% select(maxFRic, BigPA, dist_to_PA.z, 
                                        country, utm_east, 
                                        utm_north, utm_east.z, utm_north.z, 
                                        forest_structure, access_log10.z, HDI.z)
dat_FR_spill <- dat_FR_spill[complete.cases(dat_FR_spill), ]

# Perform propensity score matching again following the DAG
match_size <- matchit(BigPA ~ utm_north.z + utm_east.z + forest_structure + 
                          access_log10.z + HDI.z + dist_to_PA.z, 
                      data = dat_FR_spill, method = "full", distance = "glm", 
                      link = "probit", replace = F)
data_matched_size <- match.data(match_size)

# Run original Brodie linear mixed effects model with exponential spatial 
# correlation structure for the residuals 
mod_FR_size <- lme(maxFRic ~ forest_structure + access_log10.z + HDI.z + 
                       dist_to_PA.z + BigPA, 
                   random = list(~1 | country), data = data_matched_size, 
                   weights = ~I(1/weights), 
                   correlation = corExp(form = ~utm_east + utm_north, 
                                        nugget = TRUE))
summary(mod_FR_size)

# Peter - repeat the model above with the connectivity metrics included if we 
# think this is an analysis we wish to pursue


# Linear Mixed Effects Model of PA Distance Spillovers -----------------------------

# Modification: We will re-standardize the predictor variables after subset for
# each of the two spillover analyses, rather than just carrying through the 
# universal standardization.

# Label nearest PA that are > 2km  size threshold as large PA
dat_PA_size$CloseToPA <- ifelse(dat_PA_size$dist_to_PA > 2, 0, 1)

# Select variables for analysis and restrict to rows with complete values 

# Peter - We need to add the eventual connectivity measures to this selection 
# so they are in place for our extended analysis.
dat_FR_spill <- dat_PA_size[! dat_PA_size$station %in% 
                                dat_FR_size_outliers, ]

dat_FR_spill <- dat_FR_spill %>% 
    select(maxFRic, CloseToPA, PA_size_km2.z, 
           country, utm_east, utm_north, 
           utm_east.z, utm_north.z, 
           forest_structure, access_log10.z, HDI.z)
dat_FR_spill <- dat_FR_spill[complete.cases(dat_FR_spill), ]

# Perform propensity score matching again following the DAG
match_dist <- matchit(CloseToPA ~ utm_north.z + utm_east.z + forest_structure 
                      + access_log10.z + HDI.z + PA_size_km2.z, 
                      data = dat_FR_spill, method = "full", distance = "glm", 
                      link = "probit", replace = F)
dat_matched_dist <- match.data(match_dist)

# Run original Brodie linear mixed effects model with exponential spatial 
# correlation structure for the residuals 
mod_FR_dist <- lme(maxFRic ~ forest_structure + access_log10.z + HDI.z + 
                       PA_size_km2.z + CloseToPA, 
                   random = list(~1 | country), data = dat_matched_dist, 
                   weights = ~I(1/weights), 
                   correlation = corExp(form = ~utm_east + utm_north, 
                                        nugget = TRUE))
summary(mod_FR_dist)

# ------------------- Analysis of species richness (SR) ------------------------

# Linear Mixed Effects Model of PA Efficacy ------------------------------------

# Extension 1: We introduce an interaction term to assess the moderating effect
# of PA connectivity on the effect of PA status on biodiversity

# Modification: We may attempt to match without the lat long coordinates 
# included in the probit model.

# Replicate data frame for SR sub-analysis
dat_SR_efficacy <- subset(dat_clean, med_dist == 100)

# Remove high-leverage outliers identified by Brodie et al. 
### Brodie et al. identified outlier using the hatvalue function. We should run 
### the analysis with their outlier set removed, but also run the hatvalues 
### analysis to identify the outliers for our particular specification thereby 
### mirroring their procedure.

SR_efficacy_outliers <- c("L4789498", "L921125", "L1122096", 
                          "L7010824", "L3865754", "L3776738")
dat_SR_efficacy <- dat_SR_efficacy[! dat_SR_efficacy$station %in% 
                                       SR_efficacy_outliers, ]

# Select variables for analysis and restrict to rows with complete values 
# Peter - We need to add the eventual connectivity measures to this selection 
# so they are in place for our extended analysis

dat_SR_efficacy <- dat_SR_efficacy %>% 
    select(SR.mean, PA, country, utm_east, 
           utm_north, utm_east.z, utm_north.z, forest_structure, 
           access_log10.z, HDI.z, awf_rst_ptp2.z)
dat_SR_efficacy <- dat_SR_efficacy[complete.cases(dat_SR_efficacy), ]

# Perform propensity score matching following the DAG developed in the 
# structural causal modeling  and retrieve the matched dataset
match_SR <- matchit(PA ~ utm_east.z + utm_north.z + forest_structure + 
                        access_log10.z + HDI.z,	
                    data = dat_SR_efficacy, method = "full", 
                    distance = "glm", link = "probit", replace = F)
dat_matched_SR <- match.data(match_SR)

# Run original Brodie linear mixed effects model with exponential spatial 
# correlation structure for the residuals 
mod_SR_efficacy <- lme(SR.mean ~ forest_structure + access_log10.z 
                       + HDI.z + PA, random = list(~1 | country), 
                       data = dat_matched_SR, weights = ~I(1/weights), 
                       correlation = corExp(form = ~utm_east + utm_north, 
                                            nugget = TRUE))
summary(mod_SR_efficacy) 

# Peter - we may want to introduce a spatial autocorrelation check of the model 
# residuals

# Run linear mixed effect model with the addition of connectivity moderator 
# + conn + conn:PA
mod_SR_CN_efficacy <- lme(SR.mean ~ forest_structure + access_log10.z 
                          + HDI.z + PA + awf_rst_ptp2.z + awf_rst_ptp2.z:PA, 
                          random = list(~1 | country), data = dat_matched_SR, 
                          weights = ~I(1/weights), 
                          correlation = corExp(form = ~utm_east + utm_north, 
                                               nugget = TRUE))
summary(mod_SR_CN_efficacy)

# Summarize the two model outputs in a table
msummary(list(mod_FR_efficacy, mod_SR_CN_efficacy))


# Linear Mixed Effects Model of PA Size Spillovers -----------------------------

# Modification: We will re-standardize the predictor variables after subset for
# each of the two spillover analyses, rather than just carrying through the 
# universal standardization.

# Subset dataset to stations outside of PA
dat_PA_size <- subset(dat_clean, PA==0)

# Label nearest PA that are above a 500km2 size threshold as large PA
PA_size_threshold.z <- (500 - mean(dat_clean$PA_size_km2)) / 
    sd(dat_clean$PA_size_km2)
dat_PA_size$BigPA <- ifelse(dat_PA_size$PA_size_km2.z < PA_size_threshold.z, 
                            0, 1)

# Remove high-leverage outliers identified by Brodie et al. 
dat_SR_size_outliers <- c("L4225511", "L5624588", "L3321319", "L14087870")
dat_SR_spill <- dat_PA_size[! dat_PA_size$station %in% 
                                dat_SR_size_outliers, ]

# Select variables for analysis and restrict to rows with complete values 
# Peter - We need to add the eventual connectivity measures to this selection 
# so they are in place for our extended analysis
dat_SR_spill <- dat_SR_spill %>% 
    select(SR.mean, BigPA, dist_to_PA.z, 
           country, utm_east, 
           utm_north, utm_east.z, utm_north.z, 
           forest_structure, access_log10.z, HDI.z)
dat_SR_spill <- dat_SR_spill[complete.cases(dat_SR_spill), ]

# Perform propensity score matching again following the DAG
match_size <- matchit(BigPA ~ utm_north.z + utm_east.z + forest_structure + 
                          access_log10.z + HDI.z + dist_to_PA.z, 
                      data = dat_SR_spill, method = "full", distance = "glm", 
                      link = "probit", replace = F)
data_matched_size <- match.data(match_size)

# Run original Brodie linear mixed effects model with exponential spatial 
# correlation structure for the residuals 
mod_SR_size <- lme(SR.mean ~ forest_structure + access_log10.z + HDI.z + 
                       dist_to_PA.z + BigPA, 
                   random = list(~1 | country), data = data_matched_size, 
                   weights = ~I(1/weights), 
                   correlation = corExp(form = ~utm_east + utm_north, 
                                        nugget = TRUE))
summary(mod_SR_size)

# Peter - repeat the model above with the connectivity metrics included if we 
# think this is an analysis we wish to pursue


# Linear Mixed Effects Model of PA Distance Spillovers -----------------------------

# Modification: We will re-standardize the predictor variables after subset for
# each of the two spillover analyses, rather than just carrying through the 
# universal standardization.

# Label nearest PA that are > 2km  size threshold as large PA
dat_PA_size$CloseToPA <- ifelse(dat_PA_size$dist_to_PA > 2, 0, 1)

# Select variables for analysis and restrict to rows with complete values 

# Peter - We need to add the eventual connectivity measures to this selection 
# so they are in place for our extended analysis.
dat_SR_spill <- dat_PA_size[! dat_PA_size$station %in% 
                                dat_SR_size_outliers, ]

dat_SR_spill <- dat_SR_spill %>% 
    select(asymptPD, CloseToPA, PA_size_km2.z, 
           country, utm_east, utm_north, 
           utm_east.z, utm_north.z, 
           forest_structure, access_log10.z, HDI.z)
dat_SR_spill <- dat_SR_spill[complete.cases(dat_SR_spill), ]

# Perform propensity score matching again following the DAG
match_dist <- matchit(CloseToPA ~ utm_north.z + utm_east.z + forest_structure 
                      + access_log10.z + HDI.z + PA_size_km2.z, 
                      data = dat_SR_spill, method = "full", distance = "glm", 
                      link = "probit", replace = F)
dat_matched_dist <- match.data(match_dist)

# Run original Brodie linear mixed effects model with exponential spatial 
# correlation structure for the residuals 
mod_SR_dist <- lme(SR.mean ~ forest_structure + access_log10.z + HDI.z + 
                       PA_size_km2.z + CloseToPA, 
                   random = list(~1 | country), data = dat_matched_dist, 
                   weights = ~I(1/weights), 
                   correlation = corExp(form = ~utm_east + utm_north, 
                                        nugget = TRUE))
summary(mod_SR_dist)
