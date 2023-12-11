# Adapted from code download from https://doi.org/10.5281/zenodo.7796347
# Modified by Lei Song (lsong@ucsb.edu), Peter Kedron (peter.kedron@ucsb.edu)
# Time: Dec 5, 2023

# csv for the code:
# https://figshare.com/articles/dataset/Data/22527295

# Author Comments:
# - Species observations, trait data and phylogeny construction
##  Data source and pre-cleaning are described clear in the main text.
##  PD calculation is described very vague and they did not share any code.

# - Variables
##  I think we could have all independent variables.

# - Diversity estimation
##  They did not share the code, but they share the final training dataset used
##  for the code

# - Structural causal modelling
##  I think they just use this as a conceptual tool to select variables
##  So maybe no code at all.

# - Linear mixed-effects modelling and propensity score matching
##  The code they shared

# Load Packages ----------------------------------------------------------------
library(tidyverse)
library(cowplot)
library(here)
# who would name a package with their own name?? # good point
# functions for data cleaning, processing, analysis
library(Hmisc) 
# https://kosukeimai.github.io/MatchIt/
# matching for covariate balance in observational studies
library(MatchIt)
# https://markmfredrickson.github.io/optmatch/
# optimal full matching algorithm, for propensity score matching
library(optmatch)
# https://www.jaredknowles.com/journal/2013/11/25/getting-started-with-mixed-effect-models-in-r
# non-/Linear mixed-effect model
library(lme4) # not used
# another relatively older package to do mixed-effect model
library(nlme)
# for tests in Linear Mixed Effects Models
library(lmerTest) # not used
library(dplyr)
library(dagitty) # for drawing the directed acyclic graph (DAG)
library(ggdag) # for drawing dag

# Set working directory
setwd(here('data/raw/public/training')) 


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
ggdag_adjustment_set(tidy_dagBrodie, 
                     node_size = 20, 
                     text_col = "black"
                     ) + theme(legend.position = "bottom")


# ------------------------------------------------------------------------------
# Process Data and Construct Variables
# ------------------------------------------------------------------------------
# Lei - Hmm... The file they shared not the same as the one they used here.

# Deviation 1: The data file shared by the original authors did not include a
# human development index (HDI) measure. Following the in-text description we 
# assigned each station the HDI of it country.

# Clarification 1: Brodie et al. identified two GEDI-based measures of forest 
# structure as the best fit for their diversity data. We follow their procedure, 
# but clarify the mapping of measure to concepts by matching variable names to 
# those used in the structural causal modelling

# Extension 1: We introduce X measures of station to PA connectivity, calculated
# following the procedures defined in calc_conn_metrics.R.


# Clear the workspace
rm(list = ls())

# Load the data provided by Brodie et al. (2023)
dat_brodie <- data.frame(read.csv("bird_data_230326.csv", header = T))

# Simplify the variable names of site identifier and geographic coordinates
names(dat_brodie)[names(dat_brodie) == "site"] <- "station"
names(dat_brodie)[names(dat_brodie) == "lat_wgs84"] <- "lat"
names(dat_brodie)[names(dat_brodie) == "long_wgs84"] <- "long"

# Search for HDI in the column names 
grep("hdi", names(dat_brodie), value = TRUE)
grep("HDI", names(dat_brodie), value = TRUE)

# Assign stations the HDI value of its country
dat_HDI <- data.frame(country = unique(dat_brodie$country), 
                      HDI = c(0.593, 0.768, 0.705, 0.607, 0.803, 0.939, 0.800, 
                              0.703, 0.829)
                      )
dat_brodie <- left_join(dat_brodie, dat_HDI, by = "country")

# Create dataframe containing the subset of variable used in the analysis
dat <- dat_brodie %>% select(station, country, PA, utm_east, utm_north, 
                             Hansen_recentloss,access_log10, HDI, dist_to_PA, 
                             PA_size_km2, rh_95_a0.pred, pavd_0_5.pred, 
                             pai_a0.pred, fhd_pai_1m_a0.pred, cover_a0.pred, 
                             agbd_a0.pred,asymptPD, maxFRic, SR.mean)

# Add the connectivity variables for each station calculated with 
# calc_conn_metrics.R
# dat_conn_metrics <- data.frame(read.csv("conn_metrics.csv", header = T))
# dat <- left_join(dat, data_conn_metrics, by = "station")

# Scale subset of continuous variables in dat
# Lei - Drive me nut. The commented out variables are missing in the csv. 
### Peter deleted because these variables are also not used in the subsequent 
### analysis.

# Peter - Need to add the connectivity measure to the scaling when it is 
# introduced

dat_scale <- data.frame(scale(subset(dat, select = c("utm_east", "utm_north", 
                                                     "HDI", "access_log10", 
                                                     "PA_size_km2", 
                                                     "dist_to_PA", 
                                                     "rh_95_a0.pred", 
                                                     "pavd_0_5.pred", 
                                                     "pai_a0.pred", 
                                                     "fhd_pai_1m_a0.pred", 
                                                     "cover_a0.pred", 
                                                     "agbd_a0.pred"), 
                                     center = TRUE, scale = TRUE)))

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
# Peter changed naming so data_clean == d1 in old code
# ----------------- Analysis of phylogenetic diversity (PD) --------------------

# Linear Mixed Effects Model of PA Efficacy ------------------------------------

# Lei - Hmm...Remove, based on what justification?
### Peter - There is no statistical justification in the paper, but high 
### leverage outliers are typically removed because they exert undue control on 
### The regression estimates. Brodie et al. identified outlier using the 
### hatvalue function. We should run the analysis with their outlier set remove, 
### but also # run the hatvalues analysis to identify the outliers for our 
### particular specification thereby mirroring their procedure.

dat_clean$y <- dat_clean$asymptPD

# Replicate data frame for PD sub-analysis
dat_PD_PA_efficacy <- dat_clean

# Remove high-leverage outliers identified by Brodie et al. 
PD_PA_efficacy_outliers <- c("L2422371", "L3776738", "L2521761", "L6127181", 
                             "L3865754")
dat_PD_PA_efficacy <- dat_PD_PA_efficacy[! dat_PD_PA_efficacy$station %in% 
                                           PD_PA_efficacy_outliers, ]

# Peter changed naming so d1b == dat_PD_PA_efficacy
# Question: what is the rationality to include spatial coordinate? 
d1b <- d1b %>% select(y, PA, country, utm_east, utm_north, utm_east.z, utm_north.z, gedi.m, access_log10.z, HDI.z)
d1b <- d1b[complete.cases(d1b), ]
matchb <- MatchIt::matchit(PA ~ utm_east.z + utm_north.z + gedi.m + access_log10.z + HDI.z, 
	data = d1b, method = "full", distance = "glm", link = "probit", replace = F)
dmatchb <- match.data(matchb)

mm02b <- nlme::lme(y ~ gedi.m + access_log10.z + HDI.z + PA, 
                   random = list(~1 | country), 
                   data = dmatchb, weights = ~I(1/weights), 
	correlation = corExp(form = ~utm_east + utm_north, nugget = TRUE))
summary(mm02b) 

#----- Outside PAs - 'PA size' effect - Matched -----
d3 <- subset(d1, PA == 0)
d3$PA_size_km2 <- (d3$PA_size_km2.z * s44) + m44
PA_size_threshold.z <- (500 - m44) / s44
d3$BigPA <- ifelse(d3$PA_size_km2.z < PA_size_threshold.z, 0, 1)


#--- Matched dataset
d3b <- d3
# Remove high-leverage outliers
d3b <- subset(d3b, station != "L1084299")
d3b <- subset(d3b, station != "L4225511")
d3b <- subset(d3b, station != "L3846512")
d3b <- subset(d3b, station != "L2129865")
d3b <- subset(d3b, station != "L3267752")

d3b <- d3b %>% select(y, BigPA, dist_to_PA.z, PA_size_km2.z, country, utm_east, utm_north, utm_east.z, utm_north.z, gedi.m, access_log10.z, HDI.z)
d3b <- d3b[complete.cases(d3b), ]

matchb <- MatchIt::matchit(BigPA ~ utm_north.z + utm_east.z + gedi.m + access_log10.z + HDI.z + 
	dist_to_PA.z, data = d3b, method = "full", distance = "glm",
	link = "probit", replace = F)
dmatchb <- match.data(matchb)
mm02b <- nlme::lme(y ~ gedi.m + access_log10.z + HDI.z + dist_to_PA.z + BigPA, random = list(~1 | country),
	data = dmatchb, weights = ~I(1/weights), correlation = corExp(form = ~utm_east + utm_north, nugget = TRUE))
summary(mm02b)

#----- Outside PAs - 'Distance to PA' effect - Matched -----
d3b$dist_to_PA <- (d3b$dist_to_PA.z * s45) + m45
PA_dist_threshold <- 2 # threshold distance (km)
d3b$CloseToPA <- ifelse(d3b$dist_to_PA > PA_dist_threshold, 0, 1)

d3b <- d3b %>% select(y, CloseToPA, PA_size_km2.z, country, utm_east, utm_north, utm_east.z, utm_north.z, gedi.m, access_log10.z, HDI.z)
d3b <- d3b[complete.cases(d3b), ]

match2 <- MatchIt::matchit(CloseToPA ~ utm_north.z + utm_east.z + gedi.m + access_log10.z + HDI.z + 
	PA_size_km2.z, data = d3b, method = "full", distance = "glm",
	link = "probit", replace = F)
summary(match2)
dmatch <- match.data(match2)

mm02b <- nlme::lme(y ~ gedi.m + access_log10.z + HDI.z + PA_size_km2.z + CloseToPA, random = list(~1 | country),
	data = dmatch, weights = ~I(1/weights), correlation = corExp(form = ~utm_east + utm_north, nugget = TRUE))
summary(mm02b)

########## Below are repeat code for FR and SR ##################

############################################################################################################
#--------------------------- FR ----------------------------------------------------------------------------
d1$y <- d1$maxFR # no such variable in the csv, only maxFRic, not fully sure if they are the same

#----- All sites - PA effect - Matched -----
d1b <- d1
# Remove high-leverage outliers
d1b <- subset(d1b, station != "L921125")
d1b <- subset(d1b, station != "L2422371")
d1b <- subset(d1b, station != "L4331944")
d1b <- subset(d1b, station != "L13465594")
d1b <- d1b %>% select(y, PA, country, utm_east, utm_north, utm_east.z, utm_north.z, gedi.m, access_log10.z, HDI.z)
d1b <- d1b[complete.cases(d1b), ]
matchb <- MatchIt::matchit(PA ~ utm_east.z + utm_north.z + gedi.m + access_log10.z + HDI.z, 
	data = d1b, method = "full", distance = "glm", link = "probit", replace = F)
dmatchb <- match.data(matchb)
mm02b <- nlme::lme(y ~ gedi.m + access_log10.z + HDI.z + PA, random = list(~1 | country), data = dmatchb, weights = ~I(1/weights), 
	correlation = corExp(form = ~utm_east + utm_north, nugget = TRUE))
summary(mm02b) 

#----- Outside PAs - 'PA size' effect - Matched -----
d3 <- subset(d1, PA == 0)
d3$PA_size_km2 <- (d3$PA_size_km2.z * s44) + m44
PA_size_threshold.z <- (500 - m44) / s44
d3$BigPA <- ifelse(d3$PA_size_km2.z < PA_size_threshold.z, 0, 1)

#--- Matched dataset
d3b <- d3
# Remove high-leverage outliers
d3b <- subset(d3b, station != "L4225511")
d3b <- subset(d3b, station != "L5969878")
d3b <- subset(d3b, station != "L3267752")
d3b <- subset(d3b, station != "L4331944")
d3b <- subset(d3b, station != "L13465594")
d3b <- subset(d3b, station != "L1084299")
d3b <- d3b %>% select(y, BigPA, dist_to_PA.z, PA_size_km2.z, country, utm_east, utm_north, utm_east.z, utm_north.z, gedi.m, access_log10.z, HDI.z)
d3b <- d3b[complete.cases(d3b), ]
matchb <- MatchIt::matchit(BigPA ~ utm_north.z + utm_east.z + gedi.m + access_log10.z + HDI.z + 
	dist_to_PA.z, data = d3b, method = "full", distance = "glm",
	link = "probit", replace = F)
dmatchb <- match.data(matchb)
mm02b <- nlme::lme(y ~ gedi.m + access_log10.z + HDI.z + dist_to_PA.z + BigPA, random = list(~1 | country),
	data = dmatchb, weights = ~I(1/weights), correlation = corExp(form = ~utm_east + utm_north, nugget = TRUE))
summary(mm02b) 

#----- Outside PAs - 'Distance to PA' effect - Matched -----
d3b$dist_to_PA <- (d3b$dist_to_PA.z * s45) + m45
PA_dist_threshold <- 2 # threshold distance (km)
d3b$CloseToPA <- ifelse(d3b$dist_to_PA > PA_dist_threshold, 0, 1)

#--- Matched dataset
d3b <- d3b %>% select(y, CloseToPA, PA_size_km2.z, country, utm_east, utm_north, utm_east.z, utm_north.z, gedi.m, access_log10.z, HDI.z)
d3b <- d3b[complete.cases(d3b), ]
match2 <- MatchIt::matchit(CloseToPA ~ utm_north.z + utm_east.z + gedi.m + access_log10.z + HDI.z + 
	PA_size_km2.z, data = d3b, method = "full", distance = "glm",
	link = "probit", replace = F)
summary(match2)
dmatch <- match.data(match2)

#--- Models - Binary ---
mm02 <- nlme::lme(y ~ gedi.m + access_log10.z + HDI.z + PA_size_km2.z + CloseToPA, random = list(~1 | country),
	data = dmatch, weights = ~I(1/weights), correlation = corExp(form = ~utm_east + utm_north, nugget = TRUE))
summary(mm02)

############################################################################################################
#--------------------------- SR ----------------------------------------------------------------------------
d1$y <- d1$SR.mean

#----- All sites - PA effect - Matched -----
d1b <- d1
# Remove high-leverage outliers
d1b <- subset(d1b, station != "L4789498")
d1b <- subset(d1b, station != "L921125")
d1b <- subset(d1b, station != "L1122096")
d1b <- subset(d1b, station != "L7010824")
d1b <- subset(d1b, station != "L3865754")
d1b <- subset(d1b, station != "L3776738")
d1b <- d1b %>% select(y, PA, country, utm_east, utm_north, utm_east.z, utm_north.z, gedi.m, access_log10.z, HDI.z)
d1b <- d1b[complete.cases(d1b), ]
matchb <- MatchIt::matchit(PA ~ utm_east.z + utm_north.z + gedi.m + access_log10.z + HDI.z, 
	data = d1b, method = "full", distance = "glm", link = "probit", replace = F)
dmatchb <- match.data(matchb)
mm02b <- nlme::lme(y ~ gedi.m + access_log10.z + HDI.z + PA, random = list(~1 | country), data = dmatchb, weights = ~I(1/weights), 
	correlation = corExp(form = ~utm_east + utm_north, nugget = TRUE))
summary(mm02b) 

#----- Outside PAs - 'PA size' effect - Matched -----
d3 <- subset(d1, PA == 0)
d3$PA_size_km2 <- (d3$PA_size_km2.z * s44) + m44
PA_size_threshold.z <- (500 - m44) / s44
d3$BigPA <- ifelse(d3$PA_size_km2.z < PA_size_threshold.z, 0, 1)

#--- Matched dataset
d3b <- d3
# Remove high-leverage outliers
d3b <- subset(d3b, station != "L4225511")
d3b <- subset(d3b, station != "L5624588")
d3b <- subset(d3b, station != "L3321319")
d3b <- subset(d3b, station != "L14087870")
d3b <- d3b %>% select(y, BigPA, dist_to_PA.z, PA_size_km2.z, country, utm_east, utm_north, utm_east.z, utm_north.z, gedi.m, access_log10.z, HDI.z)
d3b <- d3b[complete.cases(d3b), ]
matchb <- MatchIt::matchit(BigPA ~ utm_north.z + utm_east.z + gedi.m + access_log10.z + HDI.z + 
	dist_to_PA.z, data = d3b, method = "full", distance = "glm",
	link = "probit", replace = F)
dmatchb <- match.data(matchb)
mm02b <- nlme::lme(y ~ gedi.m + access_log10.z + HDI.z + dist_to_PA.z + BigPA, random = list(~1 | country),
	data = dmatchb, weights = ~I(1/weights), correlation = corExp(form = ~utm_east + utm_north, nugget = TRUE))
summary(mm02b) 

#----- Outside PAs - 'Distance to PA' effect - Matched -----
d3b$dist_to_PA <- (d3b$dist_to_PA.z * s45) + m45
PA_dist_threshold <- 2 # threshold distance (km)
d3b$CloseToPA <- ifelse(d3b$dist_to_PA > PA_dist_threshold, 0, 1)

#--- Matched dataset
match2 <- MatchIt::matchit(CloseToPA ~ utm_north.z + utm_east.z + gedi.m + access_log10.z + HDI.z + 
	PA_size_km2.z, data = d3b, method = "full", distance = "glm",
	link = "probit", replace = F)
summary(match2)
dmatch <- match.data(match2)

#--- Models - Binary ---
mm02 <- nlme::lme(y ~ gedi.m + access_log10.z + HDI.z + PA_size_km2.z + CloseToPA, random = list(~1 | country),
	data = dmatch, weights = ~I(1/weights), correlation = corExp(form = ~utm_east + utm_north, nugget = TRUE))
summary(mm02)
