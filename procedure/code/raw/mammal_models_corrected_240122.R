


library(tidyverse)
library(Hmisc)
library(optmatch)
library(MatchIt)
library(lme4)
library(nlme)
library(lmerTest)



setwd('')



#--------------------------- LOAD DATA ---------------------------------------------------
rm(list = ls())
dat <- data.frame(read.csv("mammal_data_corrected_240122.csv", header = T))



#--------------------------- CLEAN DATA ---------------------------------------------------
m01<-mean(dat$utm_east,na.rm=T); s01<-sd(dat$utm_east,na.rm=T); dat$utm_east.z<-(dat$utm_east-m01)/s01
m02<-mean(dat$utm_north,na.rm=T); s02<-sd(dat$utm_north,na.rm=T); dat$utm_north.z<-(dat$utm_north-m02)/s02
m03<-mean(dat$HDI,na.rm=T); s03<-sd(dat$HDI,na.rm=T); dat$HDI.z<-(dat$HDI-m03)/s03
m06<-mean(dat$elev,na.rm=T); s06<-sd(dat$elev,na.rm=T); dat$elev.z<-(dat$elev-m06)/s06
m06f<-mean(dat$TPI,na.rm=T); s06f<-sd(dat$TPI,na.rm=T); dat$TPI.z<-(dat$TPI-m06f)/s06f
m08<-mean(dat$access_log10,na.rm=T); s08<-sd(dat$access_log10,na.rm=T); dat$access_log10.z<-(dat$access_log10-m08)/s08
m09<-mean(dat$bio01,na.rm=T); s09<-sd(dat$bio01,na.rm=T); dat$bio01.z<-(dat$bio01-m09)/s09
m44<-mean(dat$PA_size_km2,na.rm=T); s44<-sd(dat$PA_size_km2,na.rm=T); dat$PA_size_km2.z<-(dat$PA_size_km2-m44)/s44
m45<-mean(dat$dist_to_PA,na.rm=T); s45<-sd(dat$dist_to_PA,na.rm=T); dat$dist_to_PA.z<-(dat$dist_to_PA-m45)/s45
m47<-mean(dat$rh_95_a0.pred,na.rm=T); s47<-sd(dat$rh_95_a0.pred,na.rm=T); dat$rh_95_a0.pred.z<-(dat$rh_95_a0.pred-m47)/s47
m51<-mean(dat$pai_a0.pred,na.rm=T); s51<-sd(dat$pai_a0.pred,na.rm=T); dat$pai_a0.pred.z<-(dat$pai_a0.pred-m51)/s51
m53<-mean(dat$fhd_pai_1m_a0.pred,na.rm=T); s53<-sd(dat$fhd_pai_1m_a0.pred,na.rm=T); dat$fhd_pai_1m_a0.pred.z<-(dat$fhd_pai_1m_a0.pred-m53)/s53
m55<-mean(dat$cover_a0.pred,na.rm=T); s55<-sd(dat$cover_a0.pred,na.rm=T); dat$cover_a0.pred.z<-(dat$cover_a0.pred-m55)/s55
m60<-mean(dat$agbd_a0.pred,na.rm=T); s60<-sd(dat$agbd_a0.pred,na.rm=T); dat$agbd_a0.pred.z<-(dat$agbd_a0.pred-m60)/s60
dat1 <- subset(dat, select = c(study_area, station, country, utm_east, utm_north, utm_east.z, utm_north.z, 
	elev.z, TPI.z, HDI.z, access_log10.z, bio01.z, 
	PA, PA_size_km2.z, dist_to_PA.z, Hansen_recentloss, 
	rh_95_a0.pred.z, pai_a0.pred.z, fhd_pai_1m_a0.pred.z, cover_a0.pred.z, agbd_a0.pred.z,  
	maxFRic, asymptPD, SR.mean)) 
dat1 <- subset(dat1, !is.na(dat1$bio01.z))
d1 <- subset(dat1, select = -c(dist_to_PA.z))
d1 <- d1[complete.cases(d1), ]
d1$gedi.m <- d1$rh_95_a0.pred.z 
d1 <- subset(d1, Hansen_recentloss == 0)






############################################################################################################
#--------------------------- PD ----------------------------------------------------------------------------
d1$y <- d1$asymptPD



#----- All sites - PA effect - Matched -----
match2 <- MatchIt::matchit(PA ~ utm_east.z + utm_north.z + gedi.m + access_log10.z + HDI.z, 
	data = d1, method = "full", distance = "glm", link = "probit", replace = F)
summary(match2)
dmatch <- match.data(match2)
mm02 <- nlme::lme(y ~ gedi.m + access_log10.z + HDI.z + PA, 
	random = list(~1 | country, ~1 | study_area), data = dmatch, weights = ~I(1/weights), 
	correlation = corExp(form = ~utm_east + utm_north, nugget = TRUE))


#----- Outside PAs - 'PA size' effect - Matched -----
d3 <- subset(d1, PA == 0)
tmp <- subset(dat1, select = c(station, dist_to_PA.z))
d3 <- left_join(d3, tmp, by = "station")
d3$PA_size_km2 <- (d3$PA_size_km2.z * s44) + m44
PA_size_threshold.z <- (500 - m44) / s44
d3$BigPA <- ifelse(d3$PA_size_km2.z < PA_size_threshold.z, 0, 1)
d3b <- subset(d3, station != "WM-OP009")
d3b <- subset(d3b, station != "WM-HCV003")
d3b <- subset(d3b, station != "C24A25")
d3b <- subset(d3b, station != "C1A09")
d3b <- subset(d3b, station != "C1B12")
matchb <- MatchIt::matchit(BigPA ~ utm_north.z + utm_east.z + gedi.m + access_log10.z + HDI.z + 
	dist_to_PA.z, data = d3b, method = "full", distance = "glm",
	link = "probit", replace = F)
dmatchb <- match.data(matchb)
mm02b <- nlme::lme(y ~ gedi.m + access_log10.z + HDI.z + dist_to_PA.z + BigPA, random = list(~1 | country, ~1 | study_area),
	data = dmatchb, weights = ~I(1/weights), correlation = corExp(form = ~utm_east + utm_north, nugget = TRUE))



#----- Outside PAs - 'Distance to PA' effect - Matched -----
d3b$dist_to_PA <- (d3b$dist_to_PA.z * s45) + m45
PA_dist_threshold <- 2 
d3b$CloseToPA <- ifelse(d3b$dist_to_PA > PA_dist_threshold, 0, 1)
match2 <- MatchIt::matchit(CloseToPA ~ utm_north.z + utm_east.z + gedi.m + access_log10.z + HDI.z + 
	PA_size_km2.z, data = d3b, method = "full", distance = "glm",
	link = "probit", replace = F)
dmatch <- match.data(match2)
mm02 <- nlme::lme(y ~ gedi.m + access_log10.z + HDI.z + PA_size_km2.z + CloseToPA, random = list(~1 | study_area),
	data = dmatch, weights = ~I(1/weights), correlation = corExp(form = ~utm_east + utm_north, nugget = TRUE))






#####################################################################################################
#--------------------------- FR ---------------------------------------------------------------------
d1$y <- d1$maxFR



#----- All sites - PA effect - Matched -----
match2 <- MatchIt::matchit(PA ~ utm_east.z + utm_north.z + gedi.m + access_log10.z + HDI.z, 
	data = d1, method = "full", distance = "glm", link = "probit", replace = F)
dmatch <- match.data(match2)
mm02 <- nlme::lme(y ~ gedi.m + access_log10.z + HDI.z + PA, 
	random = list(~1 | country, ~1 | study_area), data = dmatch, weights = ~I(1/weights), 
	correlation = corExp(form = ~utm_east + utm_north, nugget = TRUE))



#----- Outside PAs - 'PA size' effect - Matched -----
d3 <- subset(d1, PA == 0)
tmp <- subset(dat1, select = c(station, dist_to_PA.z))
d3 <- left_join(d3, tmp, by = "station")
d3$PA_size_km2 <- (d3$PA_size_km2.z * s44) + m44
PA_size_threshold.z <- (500 - m44) / s44
d3$BigPA <- ifelse(d3$PA_size_km2.z < PA_size_threshold.z, 0, 1)
d3b <- subset(d3, station != "Bal013a")
d3b <- subset(d3b, station != "Bal017a")
d3b <- subset(d3b, station != "C1CT21")
matchb <- MatchIt::matchit(BigPA ~ utm_north.z + utm_east.z + gedi.m + access_log10.z + HDI.z + 
	dist_to_PA.z, data = d3b, method = "full", distance = "glm",
	link = "probit", replace = F)
dmatchb <- match.data(matchb)
mm02b <- nlme::lme(y ~ gedi.m + access_log10.z + HDI.z + dist_to_PA.z + BigPA, random = list(~1 | country, ~1 | study_area),
	data = dmatchb, weights = ~I(1/weights), correlation = corExp(form = ~utm_east + utm_north, nugget = TRUE))



#----- Outside PAs - 'Distance to PA' effect - Matched -----
d3b$dist_to_PA <- (d3b$dist_to_PA.z * s45) + m45
PA_dist_threshold <- 2 
d3b$CloseToPA <- ifelse(d3b$dist_to_PA > PA_dist_threshold, 0, 1)
match2 <- MatchIt::matchit(CloseToPA ~ utm_north.z + utm_east.z + gedi.m + access_log10.z + HDI.z + 
	PA_size_km2.z, data = d3b, method = "full", distance = "glm",
	link = "probit", replace = F)
dmatch <- match.data(match2)
mm02 <- nlme::lme(y ~ gedi.m + access_log10.z + HDI.z + PA_size_km2.z + CloseToPA, random = list(~1 | country, ~1 | study_area),
	data = dmatch, weights = ~I(1/weights), correlation = corExp(form = ~utm_east + utm_north, nugget = TRUE))










####################################################################################################
#--------------------------- SR --------------------------------------------------------------------
d1$y <- d1$SR.mean



#----- All sites - PA effect - Matched -----
match2 <- MatchIt::matchit(PA ~ utm_east.z + utm_north.z + gedi.m + access_log10.z + HDI.z, 
	data = d1, method = "full", distance = "glm", link = "probit", replace = F)
dmatch <- match.data(match2)
mm02 <- nlme::lme(y ~ gedi.m + access_log10.z + HDI.z + PA, 
	random = list(~1 | country, ~1 | study_area), data = dmatch, weights = ~I(1/weights), 
	correlation = corExp(form = ~utm_east + utm_north, nugget = TRUE))



#----- Outside PAs - 'PA size' effect - Matched -----
d3 <- subset(d1, PA == 0)
tmp <- subset(dat1, select = c(station, dist_to_PA.z))
d3 <- left_join(d3, tmp, by = "station")
d3$PA_size_km2 <- (d3$PA_size_km2.z * s44) + m44
PA_size_threshold.z <- (500 - m44) / s44
d3$BigPA <- ifelse(d3$PA_size_km2.z < PA_size_threshold.z, 0, 1)
d3b <- subset(d3, station != "Bal011")
d3b <- subset(d3b, station != "C1CT50")
d3b <- subset(d3b, station != "C24A25")
matchb <- MatchIt::matchit(BigPA ~ utm_north.z + utm_east.z + gedi.m + access_log10.z + HDI.z + 
	dist_to_PA.z, data = d3b, method = "full", distance = "glm",
	link = "probit", replace = F)
dmatchb <- match.data(matchb)
mm02b <- nlme::lme(y ~ gedi.m + access_log10.z + HDI.z + dist_to_PA.z + BigPA, random = list(~1 | country, ~1 | study_area),
	data = dmatchb, weights = ~I(1/weights), correlation = corExp(form = ~utm_east + utm_north, nugget = TRUE))



#----- Outside PAs - 'Distance to PA' effect - Matched -----
d3b$dist_to_PA <- (d3b$dist_to_PA.z * s45) + m45
PA_dist_threshold <- 2 
d3b$CloseToPA <- ifelse(d3b$dist_to_PA > PA_dist_threshold, 0, 1)
match2 <- MatchIt::matchit(CloseToPA ~ utm_north.z + utm_east.z + gedi.m + access_log10.z + HDI.z + 
	PA_size_km2.z, data = d3b, method = "full", distance = "glm",
	link = "probit", replace = F)
dmatch <- match.data(match2)
mm02 <- nlme::lme(y ~ gedi.m + access_log10.z + HDI.z + PA_size_km2.z + CloseToPA, random = list(~1 | country, ~1 | study_area),
	data = dmatch, weights = ~I(1/weights), correlation = corExp(form = ~utm_east + utm_north, nugget = TRUE))















