## -----------------------------------------------------------------------------
## model_pa_efficacy
## -----------------------------------------------------------------------------
#
# Purpose of script: -----------------------------------------------------------
# Create main model for PA efficacy for species richness, functional richness, 
# and phylogenetic diversity.
#
# Author:       Lei Song, Peter Kedron
# Date Created: 2024-08-13
# Last Update:  2024-08-25
# Email:        lsong@ucsb.edu

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
#
# Usage example: ---------------------------------------------------------------
# taxon <- "bird"
# conn_metrics <- 'awf_ptg'
# src_dir <- "data/raw/public"
# conn_dir <- "data/derived/public"
# dst_dir <- "data/derived/public"
# dat_clean <- clean_data(taxon, conn_metrics, src_dir, conn_dir, dst_dir)
## outliers <- "auto"
#
## Create model for phylogenetic diversity of bird with connectivity calculated
## with the median dispersal distance of 100 km.
# mod <- model_pa_efficacy(dat_clean, "connec", "bird", "asymptPD", outliers)
#
# -------------------------------------------------------------------

model_pa_efficacy <- function(dat, # leave outliers removal outside of function
                              mod_type = "connec",
                              taxon = "bird",
                              independent_variable = "asymptPD",
                              outliers){
    # Check inputs
    if (!taxon %in% c("bird", "mammal")){
        stop("Taxon must be either bird or mammal.")
    }
    if (!mod_type %in% c("brodie", "connec")){
        stop("mod_type must be either brodie or connec.")
    }
    if (!independent_variable %in% c("asymptPD", "maxFRic", "SR.mean")){
        stop("independent_variable must be one of (asymptPD, maxFRic, SR.mean).")
    }
    
    # Rename the dependent variable to have a common name
    dat <- dat %>% rename(y = all_of(independent_variable))
    
    ## Select the variables
    nms <- c("station", 'y', 'PA', 'country', 'utm_east', 'utm_north', 
             'utm_east.z', 'utm_north.z', 'forest_structure', 
             'access_log10.z', 'HDI.z', "connectivity.z")
    if("study_area" %in% names(dat)){
        nms <- append(nms, "study_area", after = 2)} else {nms <- nms}
    dat <- dat %>% select(all_of(nms))
    
    # Remove NAs
    dat <- dat[complete.cases(dat), ]
    
    # Remove outliers
    if (length(outliers) == 1 & "auto" %in% outliers){
        outliers <- identify_outliers()
    }
    dat <- dat %>% filter(!station %in% outliers)
    dat <- dat %>% select(-station)
    
    if (mod_type == "brodie"){
        # Perform propensity score matching following the DAG developed in the 
        # structural causal modeling and retrieve the matched dataset
        match_mod <- matchit(PA ~ utm_east.z + utm_north.z + forest_structure + 
                                 access_log10.z + HDI.z,
                             data = dat, method = "full", 
                             distance = "glm", link = "probit", replace = F)
        dat_matched <- match.data(match_mod)
        
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
        }
        
    } else {
        # Perform propensity score matching following the DAG developed in the 
        # structural causal modeling and retrieve the matched dataset
        match_mod <- matchit(PA ~ utm_east.z + utm_north.z + forest_structure + 
                                 access_log10.z + HDI.z + connectivity.z,
                             data = dat, method = "full", 
                             distance = "glm", link = "probit", replace = F)
        dat_matched <- match.data(match_mod)
        
        # Run linear mixed effect model with the addition of connectivity moderator 
        # + conn + conn:PA
        if (taxon == "bird"){
            mod_efficacy <- lme(
                y ~ forest_structure + access_log10.z 
                + HDI.z + PA + connectivity.z + connectivity.z:PA, 
                random = list(~1 | country), data = dat_matched, 
                weights = ~I(1/weights), 
                correlation = corExp(form = ~utm_east + utm_north, 
                                     nugget = TRUE))
        } else if (taxon == "mammal"){
            mod_efficacy <- lme(
                y ~ forest_structure + access_log10.z 
                + HDI.z + PA + connectivity.z + connectivity.z:PA, 
                random = list(~1 | country, ~1 | study_area), 
                data = dat_matched, weights = ~I(1/weights), 
                correlation = corExp(form = ~utm_east + utm_north, 
                                     nugget = TRUE))
        }
    }
    
    # Return the model
    return(mod_efficacy)
}
