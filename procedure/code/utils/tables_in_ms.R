# The script to make table 1&2 in main text
## Load libraries and functions
library(here)
setwd(here())
source(here("procedure/code/kick_off.R"))
# Load required packages and scripts
kick_off(here('procedure/code'))

## Build the models
type <- "brodie" # connec, brodie
name <- "updatePA" # replicate, updatePA
conn_metrics <- 'awf_ptg'
src_dir <- "data/raw/public"
conn_dir <- "data/derived/public"
dst_dir <- "data/derived/public"
dat_clean_bird <- clean_data("bird", conn_metrics, src_dir, conn_dir, dst_dir)
dat_clean_mammal <- clean_data("mammal", conn_metrics, src_dir, conn_dir, dst_dir)

### Make a common catalog
var_catalog <- data.frame(
    response_variable = c("asymptPD", "maxFRic", "SR.mean"), 
    name = c("PD", "FR", "SR"))

### Reproduce models for birds
dat_clean_bird <- subset(dat_clean_bird, med_dist == 100)

rpl_efficacy_bird <- lapply(1:nrow(var_catalog), function(i){
    model_pa_efficacy(dat_clean_bird, type, "bird",
                      var_catalog[[i, "response_variable"]], "auto") 
}); names(rpl_efficacy_bird) <- sprintf("mod_bird_eff_%s", var_catalog$name)

### Reproduce models for mammals
dat_clean_mammal <- subset(dat_clean_mammal, med_dist == 50)

rpl_efficacy_mammal <- lapply(1:nrow(var_catalog), function(i){
    model_pa_efficacy(dat_clean_mammal, type, "mammal",
                      var_catalog[[i, "response_variable"]], "auto")
}); names(rpl_efficacy_mammal) <- sprintf("mod_mammal_eff_%s", var_catalog$name)

### Replicate spillover models for birds w/ connectivity
rpl_spill_bird <- lapply(1:nrow(var_catalog), function(i){
    mods <- lapply(c("BigPA", "CloseToPA"), function(bnr_var){
        model_pa_spillover(dat_clean_bird, type, "bird", bnr_var, 
                           var_catalog[[i, "response_variable"]], "auto")
    })
    names(mods) <- sprintf("mod_bird_%s_%s", c("size", "dist"), 
                           var_catalog[[i, "name"]])
    mods
}); rpl_spill_bird <- do.call(c, rpl_spill_bird)

### Replicate spillover models for mammals w/ connectivity
rpl_spill_mammal <- lapply(1:nrow(var_catalog), function(i){
    mods <- lapply(c("BigPA", "CloseToPA"), function(bnr_var){
        model_pa_spillover(dat_clean_mammal, type, "mammal", bnr_var, 
                           var_catalog[[i, "response_variable"]], "auto")
    })
    names(mods) <- sprintf("mod_mammal_%s_%s", c("size", "dist"), 
                           var_catalog[[i, "name"]])
    mods
}); rpl_spill_mammal <- do.call(c, rpl_spill_mammal)

# Concatenate the model results
mods <- list("bird" = do.call(c, list(rpl_efficacy_bird, rpl_spill_bird)), 
             "mammal" = do.call(c, list(rpl_efficacy_mammal, rpl_spill_mammal)))
save(mods, file = sprintf('results/models_%s.rda', name))

## Insert code to construct a summary table of efficacy models 
## like Table 1 in Brodie et al.
### Define variable convert table
var_cvt <- data.frame(
    term = c("R2", "(Intercept)", "forest_structure", "access_log10.z", 
             "HDI.z", "PA",  "dist_to_PA.z", "BigPA", 
             "PA_size_km2.z", "CloseToPA", "connectivity.z"),
    Variable = c("R<SUP>2</SUP>", "(Intercept)", "Forest canopy height",
                 "Site accessibility", "HDI", "PA",
                 "Distance to PA", "PA size (binary)", "PA size", 
                 "Distance to PA (binary)", "Connectivity"))

coefs <- lapply(names(mods), function(taxon){
    load(sprintf("results/models_%s_brodie_orig.rda", taxon))
    models_orig <- models; rm(models); names_org <- names(models_orig)
    models <- mods[[taxon]]
    lapply(names(models), function(nm){
        # Extract names for identity
        var_nm <- toupper(strsplit(nm, "_")[[1]][4])
        effect_nm <- strsplit(nm, "_")[[1]][3]
        
        # Get the original model
        mod_to_load <- names_org[grep(tolower(var_nm), names_org)]
        mod_to_load <- mod_to_load[grep(effect_nm, mod_to_load)]
        
        effect_nm <- ifelse(
            effect_nm == "eff", "All sites",
            ifelse(effect_nm == "size", 
                   "Outside protected areas - 'PA size' effect",
                   "Outside protected areas - 'Dist to PA' effect"))
        r_square <- r.squaredGLMM(models[[nm]])[[2]]
        
        # Do the tidy work
        vals <- broom.mixed::tidy(
            models[[nm]], effects='fixed', conf.int = TRUE) %>% 
            mutate(report_value = sprintf("%.3f<br>(%.3f; %.3f)", 
                                          estimate, std.error, p.value)) %>% 
            select(term, report_value) %>% 
            rbind(c("R2", sprintf("%.3f", r_square)), .) %>% 
            left_join(var_cvt, by = "term") %>% 
            mutate(var_nm = var_nm, Effect = effect_nm) %>% 
            select(Effect, Variable, report_value, var_nm)
        
        # Calculate the difference (mean(PA) - mean(non_PA) / mean(non_PA))
        if ("BigPA" %in% names(models[[nm]]$data)){
            pred <- models_orig[[mod_to_load]]$fitted
            pred_pa <- pred[models_orig[[mod_to_load]]$data$BigPA == 1]
            pred_npa <- pred[models_orig[[mod_to_load]]$data$BigPA == 0]
            diff_orig <- (mean(pred_pa) - mean(pred_npa)) / mean(pred_npa) * 100
            
            pred <- models[[nm]]$fitted
            pred_pa <- pred[models[[nm]]$data$BigPA == 1]
            pred_npa <- pred[models[[nm]]$data$BigPA == 0]
            diff <- (mean(pred_pa) - mean(pred_npa)) / mean(pred_npa) * 100
            
        } else if ("CloseToPA" %in% names(models[[nm]]$data)){
            pred <- models_orig[[mod_to_load]]$fitted
            pred_pa <- pred[models_orig[[mod_to_load]]$data$CloseToPA == 1]
            pred_npa <- pred[models_orig[[mod_to_load]]$data$CloseToPA == 0]
            diff_orig <- (mean(pred_pa) - mean(pred_npa)) / mean(pred_npa) * 100
            
            pred <- models[[nm]]$fitted
            pred_pa <- pred[models[[nm]]$data$CloseToPA == 1]
            pred_npa <- pred[models[[nm]]$data$CloseToPA == 0]
            diff <- (mean(pred_pa) - mean(pred_npa)) / mean(pred_npa) * 100
        } else {
            pred <- models_orig[[mod_to_load]]$fitted
            pred_pa <- pred[models_orig[[mod_to_load]]$data$PA == 1]
            pred_npa <- pred[models_orig[[mod_to_load]]$data$PA == 0]
            diff_orig <- (mean(pred_pa) - mean(pred_npa)) / mean(pred_npa) * 100
            
            pred <- models[[nm]]$fitted
            pred_pa <- pred[models[[nm]]$data$PA == 1]
            pred_npa <- pred[models[[nm]]$data$PA == 0]
            diff <- (mean(pred_pa) - mean(pred_npa)) / mean(pred_npa) * 100
        }
        
        if (effect_nm == "All sites"){
            vals <- rbind(vals[1:5, ], vals[7, ], vals[6, ])
        }
        
        diffs <- data.frame(Variable = c("Reproduce", "Replicate"),
                            report_value = sprintf("%.1f%%", (c(diff_orig, diff))),
                            var_nm = var_nm, Effect = effect_nm) %>% 
            select(Effect, Variable, report_value, var_nm)
        
        rbind(vals, diffs)
    }) %>% bind_rows() %>% 
        mutate(var_nm = factor(var_nm, levels = c("SR", "FR", "PD"))) %>% 
        arrange(var_nm) %>% na.omit() %>% 
        pivot_wider(names_from = var_nm, values_from = report_value)
    
})
coefs <- left_join(coefs[[1]], coefs[[2]], by = c("Variable", "Effect"),
                   suffix = sprintf(".%s", c("Birds", "Mammals")))

coefs_efcy <- coefs %>% filter(Effect == "All sites") %>% select(-Effect)

write.csv(coefs_efcy, sprintf('results/tables/coefs_efficacy_%s.csv', name), 
          row.names = FALSE)

coefs_size <- coefs %>% 
    filter(Effect == "Outside protected areas - 'PA size' effect") %>% 
    select(-Effect)

write.csv(coefs_size, sprintf('results/tables/coefs_spillover_size_%s.csv', name), 
          row.names = FALSE)

coefs_dist <- coefs %>% 
    filter(Effect == "Outside protected areas - 'Dist to PA' effect") %>% 
    select(-Effect)

write.csv(
    coefs_dist, sprintf('results/tables/coefs_spillover_distance_%s.csv', name), 
    row.names = FALSE)
