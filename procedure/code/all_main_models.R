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
source(here("procedure/code/modeling_functions.R"))

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

# ---------------------------------- Bird --------------------------------------

# Clear the work space
rm(list = ls())
dst_dir <- here("data/derived/public")

# Load and process dataset
dat_clean <- prepare_data("bird", here("data/raw/public"),
                          here("data/derived/public"))

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
mod_orig <- modeling(dat_clean, "asymptPD", 10, 
                c("L2422371", "L3776738", "L2521761", "L6127181", "L3865754"), 
                "brodie", "bird")
bird_pd_models <- lapply(seq(10, 150, 10), function(dist){
    modeling(dat_clean, "asymptPD", dist, 
             c("L2422371", "L3776738", "L2521761", "L6127181", "L3865754"), 
             "connec", "bird")
})

bird_pd_models <- c(list(mod_orig), bird_pd_models)
names(bird_pd_models) <- c("brodie", paste("connec", seq(10, 150, 10), sep = "_"))

save(bird_pd_models, file = file.path(dst_dir, "bird_pd_models.rda"))

# Peter - repeat the model above with the connectivity metrics included if we 
# think this is an analysis we wish to pursue

# ----------------- Analysis of Functional Richness (FR) -----------------------

# Linear Mixed Effects Model of PA Efficacy ------------------------------------

# Extension 1: We introduce an interaction term to assess the moderating effect
# of PA connectivity on the effect of PA status on biodiversity

# Modification: We may attempt to match without the lat long coordinates 
# included in the probit model.
mod_orig <- modeling(dat_clean, "maxFRic", 10, 
                     c("L921125", "L2422371", "L4331944", "L13465594"), 
                     "brodie", "bird")
bird_fr_models <- lapply(seq(10, 150, 10), function(dist){
    modeling(dat_clean, "maxFRic", dist, 
             c("L921125", "L2422371", "L4331944", "L13465594"), 
             "connec", "bird")
})

bird_fr_models <- c(list(mod_orig), bird_fr_models)
names(bird_fr_models) <- c("brodie", paste("connec", seq(10, 150, 10), sep = "_"))

save(bird_fr_models, file = file.path(dst_dir, "bird_fr_models.rda"))

# ------------------- Analysis of species richness (SR) ------------------------

# Linear Mixed Effects Model of PA Efficacy ------------------------------------

# Extension 1: We introduce an interaction term to assess the moderating effect
# of PA connectivity on the effect of PA status on biodiversity

# Modification: We may attempt to match without the lat long coordinates 
# included in the probit model.

mod_orig <- modeling(dat_clean, "SR.mean", 10, 
                     c("L4789498", "L921125", "L1122096", 
                       "L7010824", "L3865754", "L3776738"), 
                     "brodie", "bird")
bird_sr_models <- lapply(seq(10, 150, 10), function(dist){
    modeling(dat_clean, "SR.mean", dist, 
             c("L4789498", "L921125", "L1122096", 
               "L7010824", "L3865754", "L3776738"), 
             "connec", "bird")
})

bird_sr_models <- c(list(mod_orig), bird_sr_models)
names(bird_sr_models) <- c("brodie", paste("connec", seq(10, 150, 10), sep = "_"))

save(bird_sr_models, file = file.path(dst_dir, "bird_sr_models.rda")) 

# --------------------------------- Mammal -------------------------------------

# Clear the work space
rm(list = ls())
dst_dir <- here("data/derived/public")

# Load and process dataset
dat_clean <- prepare_data("mammal", here("data/raw/public"),
                          here("data/derived/public"))

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
mod_orig <- modeling(dat_clean, "asymptPD", 10, NULL, "brodie", "mammal")
mammal_pd_models <- lapply(seq(10, 150, 10), function(dist){
    modeling(dat_clean, "asymptPD", dist, NULL, "connec", "mammal")
})

mammal_pd_models <- c(list(mod_orig), mammal_pd_models)
names(mammal_pd_models) <- c("brodie", paste("connec", seq(10, 150, 10), sep = "_"))

save(mammal_pd_models, file = file.path(dst_dir, "mammal_pd_models.rda"))

# Peter - repeat the model above with the connectivity metrics included if we 
# think this is an analysis we wish to pursue

# ----------------- Analysis of Functional Richness (FR) -----------------------

# Linear Mixed Effects Model of PA Efficacy ------------------------------------

# Extension 1: We introduce an interaction term to assess the moderating effect
# of PA connectivity on the effect of PA status on biodiversity

# Modification: We may attempt to match without the lat long coordinates 
# included in the probit model.
mod_orig <- modeling(dat_clean, "maxFRic", 10, NULL, "brodie", "mammal")
mammal_fr_models <- lapply(seq(10, 150, 10), function(dist){
    modeling(dat_clean, "maxFRic", dist, NULL, "connec", "mammal")
})

mammal_fr_models <- c(list(mod_orig), mammal_fr_models)
names(mammal_fr_models) <- c("brodie", paste("connec", seq(10, 150, 10), sep = "_"))

save(mammal_fr_models, file = file.path(dst_dir, "mammal_fr_models.rda"))

# ------------------- Analysis of species richness (SR) ------------------------

# Linear Mixed Effects Model of PA Efficacy ------------------------------------

# Extension 1: We introduce an interaction term to assess the moderating effect
# of PA connectivity on the effect of PA status on biodiversity

# Modification: We may attempt to match without the lat long coordinates 
# included in the probit model.

mod_orig <- modeling(dat_clean, "SR.mean", 10, NULL, "brodie", "mammal")
mammal_sr_models <- lapply(seq(10, 150, 10), function(dist){
    modeling(dat_clean, "SR.mean", dist, NULL, "connec", "mammal")
})

mammal_sr_models <- c(list(mod_orig), mammal_sr_models)
names(mammal_sr_models) <- c("brodie", paste("connec", seq(10, 150, 10), sep = "_"))

save(mammal_sr_models, file = file.path(dst_dir, "mammal_sr_models.rda"))
