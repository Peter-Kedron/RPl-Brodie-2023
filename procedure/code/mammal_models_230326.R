# Code download from https://doi.org/10.5281/zenodo.7796347
# Modified by Lei Song (lsong@ucsb.edu)
# Time: Nov 25, 2023

# Libraries
library(tidyverse)
library(cowplot)
library(dplyr)
library(here)
library(Hmisc)
library(MatchIt)
library(optmatch)
library(nlme)

setwd('data/raw/public/training')

#--------------------------- LOAD DATA ---------------------------------------------------
rm(list = ls())
dat <- data.frame(read.csv("mammal_data_230326.csv", header = T))
names(dat)[names(dat) == "lat_wgs84"] <- "lat"
names(dat)[names(dat) == "long_wgs84"] <- "long"

# HDI, the same for birds
# Reference:
# Human Development Report 2020: The Next Frontier—Human Development and the 
# Anthropocene (United Nations Development Programme, 2020).
# website: https://hdr.undp.org/data-center/human-development-index#/indicies/HDI
# https://hdr.undp.org/sites/default/files/2021-22_HDR/HDR21-22_Statistical_Annex_HDI_Table.xlsx
HDI <- data.frame(
    country = unique(dat$country), 
    HDI = c(0.803, 0.768, 0.800, 0.705, 0.939,  0.703))
dat <- left_join(dat, HDI, by = "country")

dat <- dat %>% select(
    station, study_area, country, PA, utm_east, utm_north, Hansen_recentloss,
    access_log10, HDI, dist_to_PA, PA_size_km2, rh_95_a0.pred, 
    pavd_0_5.pred, pai_a0.pred, fhd_pai_1m_a0.pred, cover_a0.pred, agbd_a0.pred,
    asymptPD, maxFRic, SR.mean)

#--------------------------- CLEAN DATA ---------------------------------------------------
# Standardize continuous variables
m01<-mean(dat$utm_east,na.rm=T); s01<-sd(dat$utm_east,na.rm=T); dat$utm_east.z<-(dat$utm_east-m01)/s01
m02<-mean(dat$utm_north,na.rm=T); s02<-sd(dat$utm_north,na.rm=T); dat$utm_north.z<-(dat$utm_north-m02)/s02
m03<-mean(dat$HDI,na.rm=T); s03<-sd(dat$HDI,na.rm=T); dat$HDI.z<-(dat$HDI-m03)/s03
MalaysiaHDI <- (0.81 - m03)/s03 # have no idea this value is used where?
m08<-mean(dat$access_log10,na.rm=T); s08<-sd(dat$access_log10,na.rm=T); dat$access_log10.z<-(dat$access_log10-m08)/s08
m44<-mean(dat$PA_size_km2,na.rm=T); s44<-sd(dat$PA_size_km2,na.rm=T); dat$PA_size_km2.z<-(dat$PA_size_km2-m44)/s44
m45<-mean(dat$dist_to_PA,na.rm=T); s45<-sd(dat$dist_to_PA,na.rm=T); dat$dist_to_PA.z<-(dat$dist_to_PA-m45)/s45
m47<-mean(dat$rh_95_a0.pred,na.rm=T); s47<-sd(dat$rh_95_a0.pred,na.rm=T); dat$rh_95_a0.pred.z<-(dat$rh_95_a0.pred-m47)/s47
m49<-mean(dat$pavd_0_5.pred,na.rm=T); s49<-sd(dat$pavd_0_5.pred,na.rm=T); dat$pavd_0_5.pred.z<-(dat$pavd_0_5.pred-m49)/s49
m51<-mean(dat$pai_a0.pred,na.rm=T); s51<-sd(dat$pai_a0.pred,na.rm=T); dat$pai_a0.pred.z<-(dat$pai_a0.pred-m51)/s51
m53<-mean(dat$fhd_pai_1m_a0.pred,na.rm=T); s53<-sd(dat$fhd_pai_1m_a0.pred,na.rm=T); dat$fhd_pai_1m_a0.pred.z<-(dat$fhd_pai_1m_a0.pred-m53)/s53
m55<-mean(dat$cover_a0.pred,na.rm=T); s55<-sd(dat$cover_a0.pred,na.rm=T); dat$cover_a0.pred.z<-(dat$cover_a0.pred-m55)/s55
m60<-mean(dat$agbd_a0.pred,na.rm=T); s60<-sd(dat$agbd_a0.pred,na.rm=T); dat$agbd_a0.pred.z<-(dat$agbd_a0.pred-m60)/s60

######################################################################################
# It is not wise to remove NA with all these non-used columns, so ignore these lines #
# Subset of data
# d1 <- dat
# d1 <- subset(d1, select = -c(PA_cat, dist_to_PA.z, PA_edge_effect.z))
# d1 <- d1[complete.cases(d1), ]
######################################################################################

######################### Exploratory analysis ###############################
#----- Find the best GEDI variable
# t1 <- lm(asymptPD ~ rh_95_a0.pred.z, data = d1)
# t2 <- lm(asymptPD ~ pavd_0_5.pred.z, data = d1)
# t3 <- lm(asymptPD ~ pai_a0.pred.z, data = d1)
# t4 <- lm(asymptPD ~ fhd_pai_1m_a0.pred.z, data = d1)
# t5 <- lm(asymptPD ~ cover_a0.pred.z, data = d1)
# t6 <- lm(asymptPD ~ agbd_a0.pred.z, data = d1)
# AIC(t1)
# AIC(t2)
# AIC(t3)
# AIC(t4) 
# AIC(t5) 
# AIC(t6)
##############################################################################
d1 <- dat
d1$gedi.m <- d1$rh_95_a0.pred.z 
d1$understory.m <- d1$pavd_0_5.pred.z

#----- Recently deforested sites
# dim(d1)
# tmp <- subset(d1, Hansen_recentloss == 1)
# dim(tmp)
d1 <- subset(d1, Hansen_recentloss == 0)

############################################################################################################
#--------------------------- PD ----------------------------------------------------------------------------
d1$y <- d1$asymptPD

#----- All sites - PA effect - Matched -----
d1b <- d1 %>% select(y, PA, study_area, country, utm_east, utm_north, utm_east.z, utm_north.z, gedi.m, access_log10.z, HDI.z)
d1b <- d1b[complete.cases(d1b), ]
match2 <- MatchIt::matchit(PA ~ utm_east.z + utm_north.z + gedi.m + access_log10.z + HDI.z, 
	data = d1b, method = "full", distance = "glm", link = "probit", replace = F)
dmatch <- match.data(match2)

#--- Analysis using matched data
mm02 <- nlme::lme(y ~ gedi.m + access_log10.z + HDI.z + PA, 
	random = list(~1 | country, ~1 | study_area), data = dmatch, weights = ~I(1/weights), 
	correlation = corExp(form = ~utm_east + utm_north, nugget = TRUE))
summary(mm02)

#----- Outside PAs - 'PA size' effect - Matched -----
d3 <- subset(d1, PA == 0)
d3$PA_size_km2 <- (d3$PA_size_km2.z * s44) + m44
PA_size_threshold.z <- (500 - m44) / s44
d3$BigPA <- ifelse(d3$PA_size_km2.z < PA_size_threshold.z, 0, 1)

#--- Models - Binary ---
d3b <- d3
# Remove high-leverage outliers
d3b <- subset(d3b, station != "WM-OP009")
d3b <- subset(d3b, station != "WM-HCV003")
d3b <- subset(d3b, station != "C24A25")
d3b <- subset(d3b, station != "C1A09")
d3b <- subset(d3b, station != "C1B12")
d3b <- d3b %>% select(y, BigPA, dist_to_PA.z, PA_size_km2.z, study_area, country, utm_east, utm_north, utm_east.z, utm_north.z, gedi.m, access_log10.z, HDI.z)
d3b <- d3b[complete.cases(d3b), ]
matchb <- MatchIt::matchit(BigPA ~ utm_north.z + utm_east.z + gedi.m + access_log10.z + HDI.z + 
	dist_to_PA.z, data = d3b, method = "full", distance = "glm",
	link = "probit", replace = F)
dmatchb <- match.data(matchb)
mm02b <- nlme::lme(y ~ gedi.m + access_log10.z + HDI.z + dist_to_PA.z + BigPA, random = list(~1 | country, ~1 | study_area),
	data = dmatchb, weights = ~I(1/weights), correlation = corExp(form = ~utm_east + utm_north, nugget = TRUE))
summary(mm02b) 

#####################################################################################################
#--------------------------- FR ---------------------------------------------------------------------
d1$y <- d1$maxFR

#----- All sites - PA effect - Matched -----
d1b <- d1 %>% select(y, PA, study_area, country, utm_east, utm_north, utm_east.z, utm_north.z, gedi.m, access_log10.z, HDI.z)
d1b <- d1b[complete.cases(d1b), ]
match2 <- MatchIt::matchit(PA ~ utm_east.z + utm_north.z + gedi.m + access_log10.z + HDI.z, 
	data = d1b, method = "full", distance = "glm", link = "probit", replace = F)
dmatch <- match.data(match2)

mm02 <- nlme::lme(y ~ gedi.m + access_log10.z + HDI.z + PA, 
	random = list(~1 | country, ~1 | study_area), data = dmatch, weights = ~I(1/weights), 
	correlation = corExp(form = ~utm_east + utm_north, nugget = TRUE))
summary(mm02)

#----- Outside PAs - 'PA size' effect - Matched -----
d3 <- subset(d1, PA == 0)
d3$PA_size_km2 <- (d3$PA_size_km2.z * s44) + m44
PA_size_threshold.z <- (500 - m44) / s44
d3$BigPA <- ifelse(d3$PA_size_km2.z < PA_size_threshold.z, 0, 1)

#--- Matched dataset
d3b <- d3
# Remove high-leverage outliers
d3b <- subset(d3b, station != "Bal013a")
d3b <- subset(d3b, station != "Bal017a")
d3b <- subset(d3b, station != "C1CT21")
d3b <- d3b %>% select(y, BigPA, dist_to_PA.z, PA_size_km2.z, study_area, country, utm_east, utm_north, utm_east.z, utm_north.z, gedi.m, access_log10.z, HDI.z)
d3b <- d3b[complete.cases(d3b), ]
matchb <- MatchIt::matchit(BigPA ~ utm_north.z + utm_east.z + gedi.m + access_log10.z + HDI.z + 
	dist_to_PA.z, data = d3b, method = "full", distance = "glm",
	link = "probit", replace = F)
dmatchb <- match.data(matchb)
mm02b <- nlme::lme(y ~ gedi.m + access_log10.z + HDI.z + dist_to_PA.z + BigPA, random = list(~1 | country, ~1 | study_area),
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
dmatch <- match.data(match2)

#--- Models - Binary ---
mm02 <- nlme::lme(y ~ gedi.m + access_log10.z + HDI.z + PA_size_km2.z + CloseToPA, random = list(~1 | country, ~1 | study_area),
	data = dmatch, weights = ~I(1/weights), correlation = corExp(form = ~utm_east + utm_north, nugget = TRUE))
summary(mm02)
MuMIn::r.squaredGLMM(mm02) # R2c is conditional R2, explaining var in fixed+random effects


####################################################################################################
#--------------------------- SR --------------------------------------------------------------------
d1$y <- d1$SR.mean

#----- All sites - PA effect - Matched -----
d1b <- d1 %>% select(y, PA, study_area, country, utm_east, utm_north, utm_east.z, utm_north.z, gedi.m, access_log10.z, HDI.z)
d1b <- d1b[complete.cases(d1b), ]
match2 <- MatchIt::matchit(PA ~ utm_east.z + utm_north.z + gedi.m + access_log10.z + HDI.z, 
	data = d1b, method = "full", distance = "glm", link = "probit", replace = F)
dmatch <- match.data(match2)

#--- Analysis using matched data
mm02 <- nlme::lme(y ~ gedi.m + access_log10.z + HDI.z + PA, 
	random = list(~1 | country, ~1 | study_area), data = dmatch, weights = ~I(1/weights), 
	correlation = corExp(form = ~utm_east + utm_north, nugget = TRUE))
summary(mm02)

#----- Outside PAs - 'PA size' effect - Matched -----
d3 <- subset(d1, PA == 0)
d3$PA_size_km2 <- (d3$PA_size_km2.z * s44) + m44
PA_size_threshold.z <- (500 - m44) / s44
d3$BigPA <- ifelse(d3$PA_size_km2.z < PA_size_threshold.z, 0, 1)

#--- Matched dataset
d3b <- d3
d3b <- subset(d3b, station != "Bal011")
d3b <- subset(d3b, station != "C1CT50")
d3b <- subset(d3b, station != "C24A25")
d3b <- d3b %>% select(y, BigPA, dist_to_PA.z, PA_size_km2.z, study_area, country, utm_east, utm_north, utm_east.z, utm_north.z, gedi.m, access_log10.z, HDI.z)
d3b <- d3b[complete.cases(d3b), ]
matchb <- MatchIt::matchit(BigPA ~ utm_north.z + utm_east.z + gedi.m + access_log10.z + HDI.z + 
	dist_to_PA.z, data = d3b, method = "full", distance = "glm",
	link = "probit", replace = F)
dmatchb <- match.data(matchb)
mm02b <- nlme::lme(y ~ gedi.m + access_log10.z + HDI.z + dist_to_PA.z + BigPA, random = list(~1 | country, ~1 | study_area),
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
dmatch <- match.data(match2)

#--- Models - Binary ---
mm02 <- nlme::lme(y ~ gedi.m + access_log10.z + HDI.z + PA_size_km2.z + CloseToPA, random = list(~1 | country, ~1 | study_area),
	data = dmatch, weights = ~I(1/weights), correlation = corExp(form = ~utm_east + utm_north, nugget = TRUE))
summary(mm02)
