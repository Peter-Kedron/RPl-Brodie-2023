## -----------------------------------------------------------------------------
## model_pa_spillover
## -----------------------------------------------------------------------------

# Purpose of script: -----------------------------------------------------------
# Create model for PA spillover for species richness,functional richness, and 
# phylogenetic diversity.
#
# Author:       Lei Song, Peter Kedron
# Date Created: 2024-08-13
# Last Update:  2024-08-25
# Email:        lsong@ucsb.edu
#
# Import from package: tidyverse, Hmisc, MatchIt, optmatch, lme4, nlme, 
#                      lmerTest, cowplot

# Inputs: ---------------------------------------------------------------------
# dat (data.frame): The data.frame of modeling data. The user should preprocess
#                   it (e.g. outlier removal) before feeding it into the model.
#
# mod_type (character): The model type. It is either brodie for running model
#                       defined in brodie paper, or connec for running model 
#                       with connectivity included.
#
# taxon (character): The taxon name. Either mammal or bird.
#
# binary_var (character): The binary exposure variable to use. Either BigPA for 
#                         PA size spillover or CloseToPA for distance to PA 
#                         spillover.
#
# independent_variable (character): The independent variable to read. It is one
#                                   of ["asymptPD", "maxFRic", "SR.mean"] 
#                                   for [phylogenetic diversity, Functional 
#                                   Richness, species richness]. 
#
#                                   The default is "asymptPD".
#
# outliers (character or vector): The user could pass a vector of outliers to
#                                 remove or set "auto" to call function 
#                                 identify_outliers to automatically detect 
#                                 outliers and remove them.
#
# Outputs: --------------------------------------------------------------------
## The lme model object.

# Usage example: ---------------------------------------------------------------
# 
# taxon <- "bird"
# conn_metrics <- 'awf_ptg'
# src_dir <- "data/raw/public"
# conn_dir <- "data/derived/public"
# dst_dir <- "data/derived/public"
# dat_clean <- clean_data(taxon, conn_metrics, src_dir, conn_dir, dst_dir)
## outliers <- NULL
#
## Create model for phylogenetic diversity of bird with connectivity calculated
## with the median dispersal distance of 100 km.
# mod <- model_pa_spillover(dat_clean, "connec", "bird", "BigPA", 
#                           "asymptPD", outliers)
#
## -----------------------------------------------------------------------------


model_pa_spillover <- function(dat, # leave outliers removal outside of function
                               mod_type = "connec", # 3A
                               taxon = "bird", # 2T
                               binary_var = "BigPA", # 2M 
                               independent_variable = "asymptPD", # 3RV
                               outliers){
    # Check inputs
    if (!taxon %in% c("bird", "mammal")){
        stop("Taxon must be either bird or mammal.")
    }
    if (!mod_type %in% c("brodie", "connec", "connec+")){
        stop("mod_type must be one of (brodie, connec, connec+).")
    }
    if (!binary_var %in% c("BigPA", "CloseToPA")){
        stop("binary_var must be either BigPA or CloseToPA.")
    }
    if (!independent_variable %in% c("asymptPD", "maxFRic", "SR.mean")){
        stop("independent_variable must be one of (asymptPD, maxFRic, SR.mean).")
    }
    
    # Subset the modeling data
    dat <- subset(dat, PA == 0)
    
    # Rename the dependent variable to have a common name
    dat <- dat %>% rename(y = all_of(independent_variable))
    
    ## select the variables
    nms <- c("station", 'y', 'country', 'utm_east', 'utm_north', 
             'utm_east.z', 'utm_north.z', 'forest_structure', 
             'access_log10.z', 'HDI.z', "connectivity.z")
    if("study_area" %in% names(dat)){
        nms <- append(nms, "study_area", after = 2)} else {nms <- nms}
    if (binary_var == "BigPA"){
        nms <- append(nms, c(binary_var, "dist_to_PA.z"), after = 1)
    } else{
        nms <- append(nms, c(binary_var, "PA_size_km2.z"), after = 1)
    }
    
    dat <- dat %>% select(all_of(nms))
    
    # Remove NAs
    dat <- dat[complete.cases(dat), ]
    
    # Remove outliers
    if (length(outliers) == 1 & "auto" %in% outliers){
        outliers <- identify_outliers()
    }
    dat <- dat %>% filter(!station %in% outliers)
    dat <- dat %>% select(-station)
    
    # Internal note: 2M
    if (binary_var == "BigPA"){
        # Create mix linear model accordingly
        ## Internal note: 3A
        if (mod_type == "brodie"){
            # Perform propensity score matching following the DAG developed in the 
            # structural causal modeling and retrieve the matched dataset
            match_mod <- matchit(BigPA ~ utm_east.z + utm_north.z + forest_structure + 
                                     access_log10.z + HDI.z + dist_to_PA.z,
                                 data = dat, method = "full", 
                                 distance = "glm", link = "probit", replace = F)
            dat_matched <- match.data(match_mod)
            
            # Run original Brodie linear mixed effects model with exponential spatial 
            # correlation structure for the residuals
            ## Internal note: 2T
            if (taxon == "bird"){
                mod_spillover <- lme(
                    y ~ forest_structure + access_log10.z 
                    + HDI.z + dist_to_PA.z + BigPA,  
                    random = list(~1 | country), 
                    data = dat_matched, weights = ~I(1/weights), 
                    correlation = corExp(form = ~utm_east + utm_north, 
                                         nugget = TRUE))
            } else if (taxon == "mammal"){
                mod_spillover <- lme(
                    y ~ forest_structure + access_log10.z 
                    + HDI.z + dist_to_PA.z + BigPA, 
                    random = list(~1 | country, ~1 | study_area), 
                    data = dat_matched, weights = ~I(1/weights), 
                    correlation = corExp(form = ~utm_east + utm_north, 
                                         nugget = TRUE))
            }
        } else {
            # Perform propensity score matching following the DAG developed in the 
            # structural causal modeling and retrieve the matched dataset
            match_mod <- matchit(
                BigPA ~ utm_east.z + utm_north.z + forest_structure + 
                    access_log10.z + HDI.z + dist_to_PA.z + connectivity.z,
                data = dat, method = "full", 
                distance = "glm", link = "probit", replace = F)
            dat_matched <- match.data(match_mod)
            
            if (mod_type == "connec"){
                # Run original Brodie linear mixed effects model with exponential spatial 
                # correlation structure for the residuals
                ## Internal note: 2T
                if (taxon == "bird"){
                    mod_spillover <- lme(
                        y ~ forest_structure + access_log10.z 
                        + HDI.z + dist_to_PA.z + BigPA + connectivity.z,  
                        random = list(~1 | country), 
                        data = dat_matched, weights = ~I(1/weights), 
                        correlation = corExp(form = ~utm_east + utm_north, 
                                             nugget = TRUE),
                        control =list(msMaxIter = 1000, msMaxEval = 1000))
                } else if (taxon == "mammal"){
                    mod_spillover <- lme(
                        y ~ forest_structure + access_log10.z 
                        + HDI.z + dist_to_PA.z + BigPA + connectivity.z, 
                        random = list(~1 | country, ~1 | study_area), 
                        data = dat_matched, weights = ~I(1/weights), 
                        correlation = corExp(form = ~utm_east + utm_north, 
                                             nugget = TRUE),
                        control =list(msMaxIter = 1000, msMaxEval = 1000))
                }
            } else {
                # Run original Brodie linear mixed effects model with exponential spatial 
                # correlation structure for the residuals
                ## Internal note: 2T
                if (taxon == "bird"){
                    mod_spillover <- lme(
                        y ~ forest_structure + access_log10.z 
                        + HDI.z + dist_to_PA.z + BigPA +
                            connectivity.z + connectivity.z:BigPA,  
                        random = list(~1 | country), 
                        data = dat_matched, weights = ~I(1/weights), 
                        correlation = corExp(form = ~utm_east + utm_north, 
                                             nugget = TRUE),
                        control =list(msMaxIter = 1000, msMaxEval = 1000))
                } else if (taxon == "mammal"){
                    mod_spillover <- lme(
                        y ~ forest_structure + access_log10.z 
                        + HDI.z + dist_to_PA.z + BigPA +
                            connectivity.z + connectivity.z:BigPA, 
                        random = list(~1 | country, ~1 | study_area), 
                        data = dat_matched, weights = ~I(1/weights), 
                        correlation = corExp(form = ~utm_east + utm_north, 
                                             nugget = TRUE),
                        control =list(msMaxIter = 1000, msMaxEval = 1000))
                }
            }
        }
    } else { # CloseToPA
        # Create mix linear model accordingly
        if (mod_type == "brodie"){
            # Perform propensity score matching following the DAG developed in the 
            # structural causal modeling and retrieve the matched dataset
            match_mod <- matchit(CloseToPA ~ utm_east.z + utm_north.z + 
                                     forest_structure + access_log10.z + HDI.z +
                                     PA_size_km2.z,
                                 data = dat, method = "full", 
                                 distance = "glm", link = "probit", replace = F)
            dat_matched <- match.data(match_mod)
            
            # Run original Brodie linear mixed effects model with exponential spatial 
            # correlation structure for the residuals
            ## Internal note: 2T
            if (taxon == "bird"){
                mod_spillover <- lme(
                    y ~ forest_structure + access_log10.z 
                    + HDI.z + PA_size_km2.z + CloseToPA, 
                    random = list(~1 | country), 
                    data = dat_matched, weights = ~I(1/weights), 
                    correlation = corExp(form = ~utm_east + utm_north, 
                                         nugget = TRUE))
            } else if (taxon == "mammal"){
                mod_spillover <- lme(
                    y ~ forest_structure + access_log10.z 
                    + HDI.z + PA_size_km2.z + CloseToPA, 
                    random = list(~1 | country, ~1 | study_area), 
                    data = dat_matched, weights = ~I(1/weights), 
                    correlation = corExp(form = ~utm_east + utm_north, 
                                         nugget = TRUE))
            }
        } else {
            # Perform propensity score matching following the DAG developed in the 
            # structural causal modeling and retrieve the matched dataset
            match_mod <- matchit(CloseToPA ~ utm_east.z + utm_north.z + 
                                     forest_structure + access_log10.z + HDI.z +
                                     PA_size_km2.z + connectivity.z,
                                 data = dat, method = "full", 
                                 distance = "glm", link = "probit", replace = F)
            dat_matched <- match.data(match_mod)
            
            if (mod_type == "connec"){
                # Run original Brodie linear mixed effects model with exponential spatial 
                # correlation structure for the residuals
                ## Internal note: 2T
                if (taxon == "bird"){
                    mod_spillover <- lme(
                        y ~ forest_structure + access_log10.z 
                        + HDI.z + PA_size_km2.z + CloseToPA + connectivity.z, 
                        random = list(~1 | country), 
                        data = dat_matched, weights = ~I(1/weights), 
                        correlation = corExp(form = ~utm_east + utm_north, 
                                             nugget = TRUE),
                        control =list(msMaxIter = 1000, msMaxEval = 1000))
                } else if (taxon == "mammal"){
                    mod_spillover <- lme(
                        y ~ forest_structure + access_log10.z 
                        + HDI.z + PA_size_km2.z + CloseToPA + connectivity.z, 
                        random = list(~1 | country, ~1 | study_area), 
                        data = dat_matched, weights = ~I(1/weights), 
                        correlation = corExp(form = ~utm_east + utm_north, 
                                             nugget = TRUE),
                        control =list(msMaxIter = 1000, msMaxEval = 1000))
                }
            } else {
                # Run original Brodie linear mixed effects model with exponential spatial 
                # correlation structure for the residuals
                ## Internal note: 2T
                if (taxon == "bird"){
                    mod_spillover <- lme(
                        y ~ forest_structure + access_log10.z 
                        + HDI.z + PA_size_km2.z + CloseToPA + connectivity.z +
                            connectivity.z:CloseToPA, 
                        random = list(~1 | country), 
                        data = dat_matched, weights = ~I(1/weights), 
                        correlation = corExp(form = ~utm_east + utm_north, 
                                             nugget = TRUE),
                        control =list(msMaxIter = 1000, msMaxEval = 1000))
                } else if (taxon == "mammal"){
                    mod_spillover <- lme(
                        y ~ forest_structure + access_log10.z 
                        + HDI.z + PA_size_km2.z + CloseToPA + connectivity.z +
                            connectivity.z:CloseToPA, 
                        random = list(~1 | country, ~1 | study_area), 
                        data = dat_matched, weights = ~I(1/weights), 
                        correlation = corExp(form = ~utm_east + utm_north, 
                                             nugget = TRUE),
                        control =list(msMaxIter = 1000, msMaxEval = 1000))
                }
            }
        }
    }
    
    # Return the model
    return(mod_spillover)
}
