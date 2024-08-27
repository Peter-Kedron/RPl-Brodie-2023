## -------------------------------------------------------------------
## Script name: tidy_coefficients
## Purpose of script: Gather and clean the coefficients for models to compare.
## Author: Lei Song
## Date Created: 2024-08-26
## Email: lsong@ucsb.edu

## Inputs:
## src_dir (character): The directory of models
## dst_dir (character): The directory to save out the plots.

## Outputs:
## Save out a data.frame with columns: Effect, Variable, [SR|FR|PD].[Birds|Mammals].
## -------------------------------------------------------------------

tidy_coefficients <- function(src_dir,
                              dst_dir,
                              # c(Birds, Mammals)
                              med_dists = c(100, 50)){
    # Define taxon names
    taxons <- c("bird", "mammal")
    
    # Manipulate the model results
    coefs <- lapply(1:2, function(i){
        # Load associated models
        fname <- file.path(
            src_dir, sprintf("models_%s_%s.rda", taxons[i], med_dists[i]))
        load(fname)
        
        # Select relevant models
        nms <- names(models)
        nms <- nms[!str_detect(nms, "brodie")]
        nms <- nms[!str_detect(nms, "spillover_connec$")]
        
        # Define variable convert table
        var_cvt <- data.frame(
            term = c("R2", "(Intercept)", "forest_structure", "access_log10.z", 
                     "HDI.z", "PA", "connectivity.z", "PA:connectivity.z",
                     "dist_to_PA.z", "BigPA", "BigPA:connectivity.z",
                     "PA_size_km2.z", "CloseToPA", "CloseToPA:connectivity.z"),
            Variable = c("R2", "(Intercept)", "Forest canopy height",
                         "Site accessibility", "HDI", "PA", "Connectivity",
                         "PA|Connectivity", "Distance to PA", "PA size (binary)",
                         "PA size|Connectivity", "PA size", 
                         "Distance to PA (binary)", "Distance to PA|Connectivity"))
        
        # Per model
        lapply(nms, function(nm){
            # Extract names for identity
            var_nm <- toupper(strsplit(nm, "_")[[1]][1])
            effect_nm <- strsplit(nm, "_")[[1]][2]
            effect_nm <- ifelse(effect_nm == "efficacy", "All sites",
                                ifelse(effect_nm == 'size', 
                                       "Outside protected areas - 'PA size' effect",
                                       "Outside protected areas - 'Distance to PA' effect"))
            r_square <- r.squaredGLMM(models[[nm]])[[2]]
            
            # Do the tidy work
            broom.mixed::tidy(
                models[[nm]], effects='fixed', conf.int = TRUE) %>% 
                mutate(report_value = sprintf("%.3f (%.3f; %.3f)", 
                                             estimate, std.error, p.value)) %>% 
                select(term, report_value) %>% 
                rbind(c("R2", sprintf("%.3f", r_square)), .) %>% 
                left_join(var_cvt, by = "term") %>% 
                mutate(var_nm = var_nm, Effect = effect_nm) %>% 
                select(Effect, Variable, report_value, var_nm)
        }) %>% bind_rows() %>% 
            mutate(var_nm = factor(var_nm, levels = c("SR", "FR", "PD"))) %>% 
            arrange(var_nm) %>% 
            pivot_wider(names_from = var_nm, values_from = report_value)
    })
    coefs <- left_join(coefs[[1]], coefs[[2]], by = c("Variable", "Effect"),
                       suffix = sprintf(".%s", c("Birds", "Mammals")))
    
    # It is very similar to TABLE 1 in supplementary
    save(coefs, file = file.path(dst_dir, "coefs_table1_compare.rda"))
}