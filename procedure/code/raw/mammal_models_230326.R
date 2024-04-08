


library(tidyverse)
library(Hmisc)
library(MatchIt)
library(optmatch)
library(lme4)
library(nlme)
library(lmerTest)
library(cowplot)





setwd('~Analysis/Mammals')







#--------------------------- LOAD DATA ---------------------------------------------------
rm(list = ls())
dat <- data.frame(read.csv("mammal_data_230326.csv", header = T))
names(dat)[names(dat) == "lat_wgs84"] <- "lat"
names(dat)[names(dat) == "long_wgs84"] <- "long"







#--------------------------- CLEAN DATA ---------------------------------------------------
# Standardize continuous variables
m01<-mean(dat$utm_east,na.rm=T); s01<-sd(dat$utm_east,na.rm=T); dat$utm_east.z<-(dat$utm_east-m01)/s01
m02<-mean(dat$utm_north,na.rm=T); s02<-sd(dat$utm_north,na.rm=T); dat$utm_north.z<-(dat$utm_north-m02)/s02
m03<-mean(dat$HDI,na.rm=T); s03<-sd(dat$HDI,na.rm=T); dat$HDI.z<-(dat$HDI-m03)/s03
MalaysiaHDI <- (0.81 - m03)/s03
m04<-mean(dat$IHDI,na.rm=T); s04<-sd(dat$IHDI,na.rm=T); dat$IHDI.z<-(dat$IHDI-m04)/s04
m05<-mean(dat$Gini,na.rm=T); s05<-sd(dat$Gini,na.rm=T); dat$Gini.z<-(dat$Gini-m05)/s05
m06<-mean(dat$elev,na.rm=T); s06<-sd(dat$elev,na.rm=T); dat$elev.z<-(dat$elev-m06)/s06
m06e<-mean(dat$slope,na.rm=T); s06e<-sd(dat$slope,na.rm=T); dat$slope.z<-(dat$slope-m06e)/s06e
m06f<-mean(dat$TPI,na.rm=T); s06f<-sd(dat$TPI,na.rm=T); dat$TPI.z<-(dat$TPI-m06f)/s06f
m06a<-mean(dat$humfoot,na.rm=T); s06a<-sd(dat$humfoot,na.rm=T); dat$humfoot.z<-(dat$humfoot-m06a)/s06a
m06b<-mean(dat$defaunind,na.rm=T); s06b<-sd(dat$defaunind,na.rm=T); dat$defaunind.z<-(dat$defaunind-m06b)/s06b
m06c<-mean(dat$globio4,na.rm=T); s06c<-sd(dat$globio4,na.rm=T); dat$globio4.z<-(dat$globio4-m06c)/s06c
m06d<-mean(dat$land_hunt,na.rm=T); s06d<-sd(dat$land_hunt,na.rm=T); dat$land_hunt.z<-(dat$land_hunt-m06d)/s06d
m07<-mean(dat$access,na.rm=T); s07<-sd(dat$access,na.rm=T); dat$access.z<-(dat$access-m07)/s07
m08<-mean(dat$access_log10,na.rm=T); s08<-sd(dat$access_log10,na.rm=T); dat$access_log10.z<-(dat$access_log10-m08)/s08
m09<-mean(dat$bio01,na.rm=T); s09<-sd(dat$bio01,na.rm=T); dat$bio01.z<-(dat$bio01-m09)/s09
m10<-mean(dat$bio02,na.rm=T); s10<-sd(dat$bio02,na.rm=T); dat$bio02.z<-(dat$bio02-m10)/s10
m11<-mean(dat$bio03,na.rm=T); s11<-sd(dat$bio03,na.rm=T); dat$bio03.z<-(dat$bio03-m11)/s11
m12<-mean(dat$bio04,na.rm=T); s12<-sd(dat$bio04,na.rm=T); dat$bio04.z<-(dat$bio04-m12)/s12
m13<-mean(dat$bio05,na.rm=T); s13<-sd(dat$bio05,na.rm=T); dat$bio05.z<-(dat$bio05-m13)/s13
m14<-mean(dat$bio06,na.rm=T); s14<-sd(dat$bio06,na.rm=T); dat$bio06.z<-(dat$bio06-m14)/s14
m15<-mean(dat$bio07,na.rm=T); s15<-sd(dat$bio07,na.rm=T); dat$bio07.z<-(dat$bio07-m15)/s15
m16<-mean(dat$bio08,na.rm=T); s16<-sd(dat$bio08,na.rm=T); dat$bio08.z<-(dat$bio08-m16)/s16
m17<-mean(dat$bio09,na.rm=T); s17<-sd(dat$bio09,na.rm=T); dat$bio09.z<-(dat$bio09-m17)/s17
m18<-mean(dat$bio10,na.rm=T); s18<-sd(dat$bio10,na.rm=T); dat$bio10.z<-(dat$bio10-m18)/s18
m19<-mean(dat$bio11,na.rm=T); s19<-sd(dat$bio11,na.rm=T); dat$bio11.z<-(dat$bio11-m19)/s19
m20<-mean(dat$bio12,na.rm=T); s20<-sd(dat$bio12,na.rm=T); dat$bio12.z<-(dat$bio12-m20)/s20
m21<-mean(dat$bio13,na.rm=T); s21<-sd(dat$bio13,na.rm=T); dat$bio13.z<-(dat$bio13-m21)/s21
m22<-mean(dat$bio14,na.rm=T); s22<-sd(dat$bio14,na.rm=T); dat$bio14.z<-(dat$bio14-m22)/s22
m23<-mean(dat$bio15,na.rm=T); s23<-sd(dat$bio15,na.rm=T); dat$bio15.z<-(dat$bio15-m23)/s23
m24<-mean(dat$bio16,na.rm=T); s24<-sd(dat$bio16,na.rm=T); dat$bio16.z<-(dat$bio16-m24)/s24
m25<-mean(dat$bio17,na.rm=T); s25<-sd(dat$bio17,na.rm=T); dat$bio17.z<-(dat$bio17-m25)/s25
m26<-mean(dat$bio18,na.rm=T); s26<-sd(dat$bio18,na.rm=T); dat$bio18.z<-(dat$bio18-m26)/s26
m27<-mean(dat$bio19,na.rm=T); s27<-sd(dat$bio19,na.rm=T); dat$bio19.z<-(dat$bio19-m27)/s27
m28<-mean(dat$bio20,na.rm=T); s28<-sd(dat$bio20,na.rm=T); dat$bio20.z<-(dat$bio20-m28)/s28
m29<-mean(dat$bio21,na.rm=T); s29<-sd(dat$bio21,na.rm=T); dat$bio21.z<-(dat$bio21-m29)/s29
m30<-mean(dat$bio22,na.rm=T); s30<-sd(dat$bio22,na.rm=T); dat$bio22.z<-(dat$bio22-m30)/s30
m31<-mean(dat$bio23,na.rm=T); s31<-sd(dat$bio23,na.rm=T); dat$bio23.z<-(dat$bio23-m31)/s31
m32<-mean(dat$bio24,na.rm=T); s32<-sd(dat$bio24,na.rm=T); dat$bio24.z<-(dat$bio24-m32)/s32
m33<-mean(dat$bio25,na.rm=T); s33<-sd(dat$bio25,na.rm=T); dat$bio25.z<-(dat$bio25-m33)/s33
m34<-mean(dat$bio26,na.rm=T); s34<-sd(dat$bio26,na.rm=T); dat$bio26.z<-(dat$bio26-m34)/s34
m35<-mean(dat$bio27,na.rm=T); s35<-sd(dat$bio27,na.rm=T); dat$bio27.z<-(dat$bio27-m35)/s35
m36<-mean(dat$bio28,na.rm=T); s36<-sd(dat$bio28,na.rm=T); dat$bio28.z<-(dat$bio28-m36)/s36
m37<-mean(dat$bio29,na.rm=T); s37<-sd(dat$bio29,na.rm=T); dat$bio29.z<-(dat$bio29-m37)/s37
m38<-mean(dat$bio30,na.rm=T); s38<-sd(dat$bio30,na.rm=T); dat$bio30.z<-(dat$bio30-m38)/s38
m39<-mean(dat$bio31,na.rm=T); s39<-sd(dat$bio31,na.rm=T); dat$bio31.z<-(dat$bio31-m39)/s39
m40<-mean(dat$bio32,na.rm=T); s40<-sd(dat$bio32,na.rm=T); dat$bio32.z<-(dat$bio32-m40)/s40
m41<-mean(dat$bio33,na.rm=T); s41<-sd(dat$bio33,na.rm=T); dat$bio33.z<-(dat$bio33-m41)/s41
m42<-mean(dat$bio34,na.rm=T); s42<-sd(dat$bio34,na.rm=T); dat$bio34.z<-(dat$bio34-m42)/s42
m43<-mean(dat$bio35,na.rm=T); s43<-sd(dat$bio35,na.rm=T); dat$bio35.z<-(dat$bio35-m43)/s43
m44<-mean(dat$PA_size_km2,na.rm=T); s44<-sd(dat$PA_size_km2,na.rm=T); dat$PA_size_km2.z<-(dat$PA_size_km2-m44)/s44
m45<-mean(dat$dist_to_PA,na.rm=T); s45<-sd(dat$dist_to_PA,na.rm=T); dat$dist_to_PA.z<-(dat$dist_to_PA-m45)/s45
m46<-mean(dat$PA_edge_effect,na.rm=T); s46<-sd(dat$PA_edge_effect,na.rm=T); dat$PA_edge_effect.z<-(dat$PA_edge_effect-m46)/s46
m47<-mean(dat$rh_95_a0.pred,na.rm=T); s47<-sd(dat$rh_95_a0.pred,na.rm=T); dat$rh_95_a0.pred.z<-(dat$rh_95_a0.pred-m47)/s47
m48<-mean(dat$rh_95_a0.var,na.rm=T); s48<-sd(dat$rh_95_a0.var,na.rm=T); dat$rh_95_a0.var.z<-(dat$rh_95_a0.var-m48)/s48
m49<-mean(dat$pavd_0_5.pred,na.rm=T); s49<-sd(dat$pavd_0_5.pred,na.rm=T); dat$pavd_0_5.pred.z<-(dat$pavd_0_5.pred-m49)/s49
m50<-mean(dat$pavd_0_5.var,na.rm=T); s50<-sd(dat$pavd_0_5.var,na.rm=T); dat$pavd_0_5.var.z<-(dat$pavd_0_5.var-m50)/s50
m51<-mean(dat$pai_a0.pred,na.rm=T); s51<-sd(dat$pai_a0.pred,na.rm=T); dat$pai_a0.pred.z<-(dat$pai_a0.pred-m51)/s51
m52<-mean(dat$pai_a0.var,na.rm=T); s52<-sd(dat$pai_a0.var,na.rm=T); dat$pai_a0.var.z<-(dat$pai_a0.var-m52)/s52
m53<-mean(dat$fhd_pai_1m_a0.pred,na.rm=T); s53<-sd(dat$fhd_pai_1m_a0.pred,na.rm=T); dat$fhd_pai_1m_a0.pred.z<-(dat$fhd_pai_1m_a0.pred-m53)/s53
m54<-mean(dat$fhd_pai_1m_a0.var,na.rm=T); s54<-sd(dat$fhd_pai_1m_a0.var,na.rm=T); dat$fhd_pai_1m_a0.var.z<-(dat$fhd_pai_1m_a0.var-m54)/s54
m55<-mean(dat$cover_a0.pred,na.rm=T); s55<-sd(dat$cover_a0.pred,na.rm=T); dat$cover_a0.pred.z<-(dat$cover_a0.pred-m55)/s55
m56<-mean(dat$cover_a0.var,na.rm=T); s56<-sd(dat$cover_a0.var,na.rm=T); dat$cover_a0.var.z<-(dat$cover_a0.var-m56)/s56
#m57<-mean(dat$Hansen_recentloss,na.rm=T); s57<-sd(dat$Hansen_recentloss,na.rm=T); dat$Hansen_recentloss.z<-(dat$Hansen_recentloss-m57)/s57
#m58<-mean(dat$Hansen_recentgain,na.rm=T); s58<-sd(dat$Hansen_recentgain,na.rm=T); dat$Hansen_recentgain.z<-(dat$Hansen_recentgain-m58)/s58
#m59<-mean(dat$Hansen_2019forestcover,na.rm=T); s59<-sd(dat$Hansen_2019forestcover,na.rm=T); dat$Hansen_2019forestcover.z<-(dat$Hansen_2019forestcover-m59)/s59
m60<-mean(dat$agbd_a0.pred,na.rm=T); s60<-sd(dat$agbd_a0.pred,na.rm=T); dat$agbd_a0.pred.z<-(dat$agbd_a0.pred-m60)/s60
m61<-mean(dat$agbd_a0.var,na.rm=T); s61<-sd(dat$agbd_a0.var,na.rm=T); dat$agbd_a0.var.z<-(dat$agbd_a0.var-m61)/s61


# Add some polynomial terms
dat$elev2.z <- dat$elev.z^2
dat$slope2.z <- dat$slope.z^2
dat$utm_north2.z <- dat$utm_north.z^2
dat$utm_east2.z <- dat$utm_east.z^2
#dat$utm_north3.z <- dat$utm_north.z^3 # these 3rd-order spatial terms actually increase model error
#dat$utm_east3.z <- dat$utm_east.z^3


# Choose a full set of predictor and response variables
dat1 <- subset(dat, select = c(study_area, station, country, utm_east, utm_north, utm_east.z, utm_north.z, 
	elev.z, elev2.z, slope.z, slope2.z, TPI.z, 
	HDI.z, access_log10.z, 
	bio01.z, bio02.z, bio03.z, bio04.z, bio05.z, bio06.z, bio07.z, bio08.z, bio09.z, bio10.z, bio11.z, bio12.z, bio13.z, bio14.z, 
	bio15.z, bio16.z, bio17.z, bio18.z, bio19.z, bio20.z, bio21.z, bio22.z, bio23.z, bio24.z, bio25.z, bio26.z, bio27.z, bio28.z,
	bio29.z, bio30.z, bio31.z, bio32.z, bio33.z, bio34.z, bio35.z, 
	PA, PA_name, PA_year, PA_cat, PA_size_km2.z, dist_to_PA.z, PA_edge_effect.z, Hansen_recentloss, 
	rh_95_a0.pred.z, pavd_0_5.pred.z, pai_a0.pred.z, fhd_pai_1m_a0.pred.z, cover_a0.pred.z, agbd_a0.pred.z,  
	rh_95_a0.var.z, pavd_0_5.var.z, pai_a0.var.z, fhd_pai_1m_a0.var.z, cover_a0.var.z, agbd_a0.var.z, 
	asymptSR, maxFRic, asymptPD, median.body.size, SR.mean)) 
dat1 <- subset(dat1, !is.na(dat1$bio01.z))


# PCA of bioclim vars
bioclim <- subset(dat1, select = c(
	bio01.z, bio02.z, bio03.z, bio04.z, bio05.z, bio06.z, bio07.z, bio08.z, bio09.z, bio10.z,
	bio11.z, bio12.z, bio13.z, bio14.z, bio15.z, bio16.z, bio17.z, bio18.z, bio19.z, bio20.z,
	bio21.z, bio22.z, bio23.z, bio24.z, bio25.z, bio26.z, bio27.z, bio28.z, bio29.z, bio30.z,
	bio31.z, bio32.z, bio33.z, bio34.z, bio35.z))
pc <- prcomp(bioclim, scale = F)
#summary(pc) # first 2 axes represent 76.0% of variance in data
trg <- data.frame(predict(pc, bioclim))
dat1$bioclim.pc1 <- trg$PC1
dat1$bioclim.pc2 <- trg$PC2
dat1$bioclim.pc3 <- trg$PC3


# Subset of data
d1 <- dat1
d1 <- subset(d1, select = -c(PA_cat, dist_to_PA.z, PA_edge_effect.z))
d1 <- d1[complete.cases(d1), ]


#----- Find the best GEDI variable
t1 <- lm(asymptPD ~ rh_95_a0.pred.z, data = d1)
t2 <- lm(asymptPD ~ pavd_0_5.pred.z, data = d1)
t3 <- lm(asymptPD ~ pai_a0.pred.z, data = d1)
t4 <- lm(asymptPD ~ fhd_pai_1m_a0.pred.z, data = d1)
t5 <- lm(asymptPD ~ cover_a0.pred.z, data = d1)
t6 <- lm(asymptPD ~ agbd_a0.pred.z, data = d1)
AIC(t1)
AIC(t2)
AIC(t3)
AIC(t4) 
AIC(t5) 
AIC(t6) 

d1$gedi.m <- d1$rh_95_a0.pred.z 
d1$understory.m <- d1$pavd_0_5.pred.z


#----- Recently deforested sites
dim(d1)
tmp <- subset(d1, Hansen_recentloss == 1)
dim(tmp)
d1 <- subset(d1, Hansen_recentloss == 0)







############################################################################################################
#--------------------------- PD ----------------------------------------------------------------------------
d1$y <- d1$asymptPD





#----- All sites - PA effect - Matched -----
match2 <- MatchIt::matchit(PA ~ utm_east.z + utm_north.z + gedi.m + access_log10.z + HDI.z, 
	data = d1, method = "full", distance = "glm", link = "probit", replace = F)
dmatch <- match.data(match2)



#--- Analysis using matched data
mm02 <- nlme::lme(y ~ gedi.m + access_log10.z + HDI.z + PA, 
	random = list(~1 | country, ~1 | study_area), data = dmatch, weights = ~I(1/weights), 
	correlation = corExp(form = ~utm_east + utm_north, nugget = TRUE))
summary(mm02)



#----- Outside PAs - 'PA size' effect - Matched -----
d3 <- subset(d1, PA == 0)
tmp <- subset(dat1, select = c(station, dist_to_PA.z))
d3 <- left_join(d3, tmp, by = "station")
d3$PA_size_km2 <- (d3$PA_size_km2.z * s44) + m44
PA_size_threshold.z <- (500 - m44) / s44
d3$BigPA <- ifelse(d3$PA_size_km2.z < PA_size_threshold.z, 0, 1)

#--- Matched dataset
match2 <- MatchIt::matchit(BigPA ~ utm_north.z + utm_east.z + gedi.m + access_log10.z + HDI.z + 
	dist_to_PA.z, data = d3, method = "full", distance = "glm",
	link = "probit", replace = F)
dmatch <- match.data(match2)

#--- Models - Binary ---
d3b <- d3
# Remove high-leverage outliers
d3b <- subset(d3b, station != "WM-OP009")
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
summary(mm02b) 










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
summary(mm02)



#----- Outside PAs - 'PA size' effect - Matched -----
d3 <- subset(d1, PA == 0)
tmp <- subset(dat1, select = c(station, dist_to_PA.z))
d3 <- left_join(d3, tmp, by = "station")
d3$PA_size_km2 <- (d3$PA_size_km2.z * s44) + m44
PA_size_threshold.z <- (500 - m44) / s44
d3$BigPA <- ifelse(d3$PA_size_km2.z < PA_size_threshold.z, 0, 1)

#--- Matched dataset
match2 <- MatchIt::matchit(BigPA ~ utm_north.z + utm_east.z + gedi.m + access_log10.z + HDI.z + 
	dist_to_PA.z, data = d3, method = "full", distance = "glm",
	link = "probit", replace = F)
dmatch <- match.data(match2)

#--- Models - Binary ---
d3b <- d3
# Remove high-leverage outliers
d3b <- subset(d3b, station != "Bal013a")
d3b <- subset(d3b, station != "Bal017a")
d3b <- subset(d3b, station != "C1CT21")
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
match2 <- MatchIt::matchit(PA ~ utm_east.z + utm_north.z + gedi.m + access_log10.z + HDI.z, 
	data = d1, method = "full", distance = "glm", link = "probit", replace = F)
dmatch <- match.data(match2)

#--- Analysis using matched data
mm02 <- nlme::lme(y ~ gedi.m + access_log10.z + HDI.z + PA, 
	random = list(~1 | country, ~1 | study_area), data = dmatch, weights = ~I(1/weights), 
	correlation = corExp(form = ~utm_east + utm_north, nugget = TRUE))
summary(mm02)




#----- Outside PAs - 'PA size' effect - Matched -----
d3 <- subset(d1, PA == 0)
tmp <- subset(dat1, select = c(station, dist_to_PA.z))
d3 <- left_join(d3, tmp, by = "station")
d3$PA_size_km2 <- (d3$PA_size_km2.z * s44) + m44
PA_size_threshold.z <- (500 - m44) / s44
d3$BigPA <- ifelse(d3$PA_size_km2.z < PA_size_threshold.z, 0, 1)

#--- Matched dataset
match2 <- MatchIt::matchit(BigPA ~ utm_north.z + utm_east.z + gedi.m + access_log10.z + HDI.z + 
	dist_to_PA.z, data = d3, method = "full", distance = "glm",
	link = "probit", replace = F)
dmatch <- match.data(match2)

#--- Models - Binary ---
d3b <- d3
d3b <- subset(d3b, station != "Bal011")
d3b <- subset(d3b, station != "C1CT50")
d3b <- subset(d3b, station != "C24A25")
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































