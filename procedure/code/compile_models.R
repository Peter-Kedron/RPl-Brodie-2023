## -------------------------------------------------------------------
## Script name: compile_models
## Purpose of script: Compile all models for PA efficacy and spillover.
## Author: Lei Song, Peter Kedron
## Date Created: 2024-08-15
## Email: lsong@ucsb.edu

## Inputs:
## taxon (character): The taxon name. Either mammal or bird.
## outliers (character): The method to detect outliers. auto for automatic
## detection and brodie for using the values reported by brodie. Also check 
## parameter outliers in function model_pa_efficacy and model_pa_spillover.
## med_dist (integer): The median dispersal distance in km. Make sure
## the connectivity values with this dispersal distance have been calculated.
## Otherwise, run function calc_conn to calculate.
## src_dir (character): The directory of data source.
## dst_dir (character): The directory to save files to.

## Outputs:
## A list of 24 lme model object. The names are: 
## pd_efficacy_brodie, pd_efficacy_connec, 
## fr_efficacy_brodie, fr_efficacy_connec,
## sr_efficacy_brodie, sr_efficacy_connec, 
## pd_size_spillover_brodie, pd_size_spillover_connec, pd_size_spillover_connec+,
## fr_size_spillover_brodie, fr_size_spillover_connec, fr_size_spillover_connec+, 
## sr_size_spillover_brodie, sr_size_spillover_connec, sr_size_spillover_connec+, 
## pd_dist_spillover_brodie, pd_dist_spillover_connec, pd_dist_spillover_connec+, 
## fr_dist_spillover_brodie, fr_dist_spillover_connec, fr_dist_spillover_connec+, 
## sr_dist_spillover_brodie, sr_dist_spillover_connec, sr_dist_spillover_connec+
## Details
## [A]_efficacy_[B]: The PA efficacy models. A means independent variables. pd
## for phylogenetic diversity, fr for Functional Richness, and sr for
## species richness. B means model type. brodie for brodie model setting, and
## connec for including connectivity and interaction between connectivity and PA.
## [A]_[B]_spillover_[C]: The PA spillover models. A means independent variables. 
## pd for phylogenetic diversity, fr for Functional Richness, and sr for
## species richness. B means binary variable to check spillover. size for PA size
## and dist for distance to PA. C means means model type. brodie for brodie 
## model setting, connec for only including connectivity, and connec+ for 
## including both connectivity and interaction between connectivity and binary.
## -------------------------------------------------------------------

compile_models <- function(taxon = "bird",
                           outliers = "brodie", # auto or brodie
                           med_dist = 100, 
                           src_dir,
                           dst_dir){
    # Check inputs
    if (!outliers %in% c("auto", "brodie")){
        stop("outliers must be either auto or brodie.")
    }
    if (!taxon %in% c("bird", "mammal")){
        stop("Taxon must be either bird or mammal.")
    }
    
    # Give convenience for calculation but keep the meaningful name for inputs.
    dist <- med_dist
    
    # Load data
    dat <- read.csv(file.path(src_dir, sprintf("dat_analysis_%s.csv", taxon)))
    
    # Run models for efficacy by independent var, model type, and dispersal distance
    efficacy_models <- lapply(c("asymptPD", "maxFRic", "SR.mean"), function(rv){
        # Set outliers for each response variable if using brodie
        if (outliers == "brodie"){
            if (rv == "asymptPD"){
                if (taxon == "bird"){
                    outliers <- c("L2422371", "L3776738", "L2521761", 
                                  "L6127181", "L3865754")
                } else{
                    outliers <- NULL
                }
            } else if (rv == "maxFRic"){
                if (taxon == "bird"){
                    outliers <- c("L921125", "L2422371", "L4331944", "L13465594")
                } else{
                    outliers <- NULL
                }
            } else {
                if (taxon == "bird"){
                    outliers <- c("L4789498", "L921125", "L1122096", 
                                  "L7010824", "L3865754", "L3776738")
                } else{
                    outliers <- NULL
                }
            }
        }
        
        dat_clean <- subset(dat, med_dist == dist)
        mods <- lapply(c("brodie", "connec"), function(mod_type){
            model_pa_efficacy(dat_clean, mod_type, taxon, rv, outliers)
        })
        
        # Set names and return
        nm <- sprintf("%s_efficacy", tolower(str_extract(rv, "[A-Z]{2}")))
        names(mods) <- paste(nm, c("brodie", "connec"), sep = "_")
        mods
    })
    # Concatenate results
    efficacy_models <- do.call(c, efficacy_models)
    
    # Run models for spillover by binary_var, independent var, model type,
    # and dispersal distance
    spillover_models <- lapply(c("BigPA", "CloseToPA"), function(bnr_var){
        mods <- lapply(c("asymptPD", "maxFRic", "SR.mean"), 
               function(rv){
                   # Set outliers for each response variable if using brodie
                   if (outliers == "brodie"){
                       if (rv == "asymptPD"){
                           if (taxon == "bird"){
                               outliers <- c("L1084299", "L4225511", "L3846512", 
                                             "L2129865", "L3267752")
                           } else{
                               outliers <- c("WM-OP009", "WM-HCV003", "C24A25", 
                                             "C1A09", "C1B12")
                           }
                       } else if (rv == "maxFRic"){
                           if (taxon == "bird"){
                               outliers <- c("L4225511", "L5969878", "L3267752", 
                                             "L4331944", "L13465594", "L1084299")
                           } else{
                               outliers <- c("Bal013a", "Bal017a", "C1CT21")
                           }
                       } else {
                           if (taxon == "bird"){
                               outliers <- c("L4225511", "L5624588", "L3321319", 
                                             "L14087870")
                           } else{
                               outliers <- c("Bal011", "C1CT50", "C24A25")
                           }
                       }
                   }
                   
                   dat_clean <- subset(dat, med_dist == dist)
                   mods <- lapply(c("brodie", "connec", "connec+"), 
                                  function(mod_type){
                       model_pa_spillover(
                           dat_clean, mod_type, taxon, bnr_var, rv, outliers)
                   })
                   
                   # Set names and return
                   nm <- sprintf("%s_%s_spillover", 
                                 tolower(str_extract(rv, "[A-Z]{2}")),
                                 ifelse(bnr_var == "BigPA", "size", "dist"))
                   names(mods) <- paste(nm, c("brodie", "connec", "connec+"), sep = "_")
                   mods
               })
        
        mods <- do.call(c, mods)
        mods
    })
    # Concatenate results
    spillover_models <- do.call(c, spillover_models)
    
    # Concatenate models and save out
    models <- c(efficacy_models, spillover_models)
    fname <- file.path(dst_dir, sprintf("models_%s_%s.rda", taxon, med_dist))
    save(models, file = fname)
}
