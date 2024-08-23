## -------------------------------------------------------------------
## Script name: identify_outliers
## Purpose of script: A tool function to identify and optionally remove
## outliers in the cleaned data from function clean_data.
## Author: Wenxin Yang | August 22nd, 2024
## Inputs:
## dat (data.frame): The data.frame of the processed data to use. dat must have
## two columns: y for independent variable and station for station index. It
## also should include multiple feature columns to define outliers.
## keep (logical): Keep the outliers in the input data or not. TRUE to keep,
## and FALSE to remove. The default is FALSE.

## Outputs:
## outliers (list): A list of outliers with station ids.
## -------------------------------------------------------------------

# Because outlier detection depends on the modeling objective, it is easier
# to deal with the data.frame directly, rather than loading csv files.
identify_outliers <- function(dat, mod_type = "connec", taxon = "bird",
                              independent_variable = "maxFRic",
                              leverage_threshold){
    
    # because as.formula (a following step) does not take random effects 
    # of string type, we have to convert it to many dummy variables
    createAllDummy <- function(a_dat_frame, li_country,
                            li_studyarea = NULL){

        # 1. create country dummies
        for (i in 1:length(li_country)){
            country_name <- li_country[i]
            varname <- paste0("if_", country_name)
            a_dat_frame[, varname] <- 0
            a_dat_frame[a_dat_frame$country == country_name, varname] <- 1
        }
        
        # 2. if mammals, we add study area dummies as well
        if(!is.null(li_studyarea)){
            for(i in 1:length(li_studyarea)){
                area_name <- li_studyarea[i]
                varname <- paste0("if_", area_name)
                a_dat_frame[, varname] <- 0
                a_dat_frame[a_dat_frame$study_area == area_name, varname] <- 1
            }
        }

        return(a_dat_frame)
    }
    
    # another function to prepare for an input of a following step
    getFormula <- function(fixed_string, a_list){
        rand_formula_string <- paste0(a_list, sep = " + ", collapse = "")
        rand_formula_string <- substr(rand_formula_string, 1, 
                                      nchar(rand_formula_string)-3)
        
        formula_mod_string <- paste0(fixed_string, rand_formula_string, ")")
        formula_mod <- as.formula(formula_mod_string)
        
        return(formula_mod)
    }
    
    # ====== 1. run a model with raw data ======
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
    
    # Remove NAs, skip outlier removal and keep "station"
    dat <- dat[complete.cases(dat), ]
    
    if (mod_type == "brodie"){
        # Perform propensity score matching following the DAG developed in the 
        # structural causal modeling and retrieve the matched dataset
        match_mod <- matchit(PA ~ utm_east.z + utm_north.z + forest_structure + 
                                 access_log10.z + HDI.z,
                             data = dat, method = "full", 
                             distance = "glm", link = "probit", replace = F)
        dat_matched <- match.data(match_mod)
        # li_country, li_studyarea, and fixed_formula_string are necessary
        # for computing hat values
        li_country <- list(unique(dat_matched$country))[[1]]
        li_studyarea <- list(unique(dat_matched$study_area))[[1]]
        fixed_formula_string <- "y ~ forest_structure + access_log10.z + HDI.z + PA + (1|"
        # Run original Brodie linear mixed effects model with exponential spatial 
        # correlation structure for the residuals 
        if (taxon == "bird"){
            li_studyarea <- NULL
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
        # li_country, li_studyarea, and fixed_formula_string are necessary
        # for computing hat values
        li_country <- list(unique(dat_matched$country))[[1]]
        li_studyarea <- list(unique(dat_matched$study_area))[[1]]
        fixed_formula_string <- "y ~ forest_structure + access_log10.z + HDI.z + PA + connectivity.z + connectivity.z:PA + (1|"
        # Run linear mixed effect model with the addition of connectivity moderator 
        # + conn + conn:PA
        if (taxon == "bird"){
            li_studyarea <- NULL
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

    # ====== 2. prep data to calculate hat values ======
    # some strings might contain "-", this would cause a problem when writing
    # the formula, so we substitute it with an underscore
    li_country <- sub("-", "_", li_country)
    li_studyarea <- sub("-", "_", li_studyarea)
    dat_hat <- createAllDummy(dat_matched, li_country, li_studyarea)
    li_country_varnames <- paste0('if_', li_country)
    if(!is.null(li_studyarea)){
        li_studyarea_varnames <- paste0('if_', li_studyarea)
    } else{
        li_studyarea_varnames <- NULL
    }
    li_all_varnames <- c(li_country_varnames, li_studyarea_varnames)
    formula_mod <- getFormula(fixed_formula_string, li_all_varnames)
    # ====== 3. run hat values ======
    mod_frame <- model.frame(formula_mod, dat_hat)
    mod_matrix <- model.matrix(mod_frame, dat_hat)
    mod_hatvals <- hat(mod_matrix)
    # ===== 4. use threshold to get outliers =====
    # length(mod_hatvals) == nrow(dat_hat)
    dat_hatvals <- as.data.frame(mod_hatvals)
    dat_hatvals$if_highlev <- 0
    dat_hatvals[dat_hatvals$mod_hatvals>=leverage_threshold, "if_highlev"] <- 1
    
    print(max(dat_hatvals$mod_hatvals))
    
    dat_matched$if_highlev <- dat_hatvals$if_highlev
    outliers <- list(dat_matched[dat_matched$if_highlev==1, "station"])[[1]]
    
    # Return
  return(outliers)
}

## tidyverse, Hmisc, MatchIt, optmatch, lme4, nlme, lmerTest, cowplot
library(tidyverse)
library(Hmisc)
library(MatchIt)
library(optmatch)
library(lme4)
library(nlme)
library(lmerTest)
library(cowplot)
library(here)

med_dist <- 100
taxon = "bird"
dat_raw <- data.frame(read.csv("~/Desktop/GitHub/RPl-Brodie-2023/data/derived/public/dat_analysis_bird.csv"), header = T)
dat_test <- dat_raw[dat_raw$med_dist == med_dist,]
leverage_threshold <- 0.02

outliers <- identify_outliers(dat = dat_test, mod_type = "connec", taxon = taxon, 
                  independent_variable = "SR.mean", # maxFRic, SR.mean
                  leverage_threshold = leverage_threshold)

