## -------------------------------------------------------------------
## Script name: kick_off
## Purpose of script: Setting up working environment before any functions
## Author: Lei Song
## Date Created: 2024-08-15
## Email: lsong@ucsb.edu

## Inputs:
## code_dir (character): The directory for all scripts.

## Outputs:
## No output. Set up the computational environment.
## -------------------------------------------------------------------

kick_off <- function(code_dir){
    # Load libraries
    pkgs <- c("here", "sf", "dplyr", "terra", "wdpar", "countrycode", "stringr", 
              "tidyverse", "cowplot", "here", "dagitty", "ggdag", "Hmisc", 
              "MatchIt", "modelsummary", "optmatch", "nlme", "MuMIn", "pbapply")
    for (pkg in pkgs){
        if (!require(pkg, character.only = TRUE)) {
            install.packages(pkg)
            require(pkg, character.only = TRUE)
        }
    }
    rm(pkgs, pkg)
    
    # Load functions
    functions <- c("clean_pa", "prep_dm", "calc_conn", "clean_data", 
                   "identify_outliers", "model_pa_efficacy", 
                   "model_pa_spillover", "compile_models")
    for (func in file.path(code_dir, sprintf("%s.R", functions))){
        source(func)
    }
}
