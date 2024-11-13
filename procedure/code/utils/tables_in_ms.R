# The script to make table 1&2 in main text
## Load libraries and functions
library(here)
setwd(here())
source(here("procedure/code/kick_off.R"))
# Load required packages and scripts
kick_off(here('procedure/code'))

## Build the models
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
    model_pa_efficacy(dat_clean_bird, "connec", "bird",
                      var_catalog[[i, "response_variable"]], "auto") 
}); names(rpl_efficacy_bird) <- sprintf("mod_bird_eff_%s", var_catalog$name)

### Reproduce models for mammals
dat_clean_mammal <- subset(dat_clean_mammal, med_dist == 50)

rpl_efficacy_mammal <- lapply(1:nrow(var_catalog), function(i){
    message(var_catalog[[i, "response_variable"]])
    model_pa_efficacy(dat_clean_mammal, "connec", "mammal",
                      var_catalog[[i, "response_variable"]], "auto")
}); names(rpl_efficacy_mammal) <- sprintf("mod_mammal_eff_%s", var_catalog$name)

### Replicate spillover models for birds w/ connectivity
rpl_spill_bird <- lapply(1:nrow(var_catalog), function(i){
    mods <- lapply(c("BigPA", "CloseToPA"), function(bnr_var){
        model_pa_spillover(dat_clean_bird, "connec", "bird", bnr_var, 
                           var_catalog[[i, "response_variable"]], "auto")
    })
    names(mods) <- sprintf("mod_bird_%s_%s", c("size", "dist"), 
                           var_catalog[[i, "name"]])
    mods
}); rpl_spill_bird <- do.call(c, rpl_spill_bird)

### Replicate spillover models for mammals w/ connectivity
rpl_spill_mammal <- lapply(1:nrow(var_catalog), function(i){
    mods <- lapply(c("BigPA", "CloseToPA"), function(bnr_var){
        model_pa_spillover(dat_clean_mammal, "connec", "mammal", bnr_var, 
                           var_catalog[[i, "response_variable"]], "auto")
    })
    names(mods) <- sprintf("mod_mammal_%s_%s", c("size", "dist"), 
                           var_catalog[[i, "name"]])
    mods
}); rpl_spill_mammal <- do.call(c, rpl_spill_mammal)

# Concatenate the model results
mods <- list("bird" = do.call(c, list(rpl_efficacy_bird, rpl_spill_bird)), 
             "mammal" = do.call(c, list(rpl_efficacy_mammal, rpl_spill_mammal)))

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
    lapply(names(models)[c(1:3, 4, 6, 8)], function(nm){
        
        # Extract names for identity
        var_nm <- toupper(strsplit(nm, "_")[[1]][4])
        effect_nm <- strsplit(nm, "_")[[1]][3]
        
        # Get the original model
        mod_to_load <- names_org[grep(tolower(var_nm), names_org)]
        mod_to_load <- mod_to_load[grep(effect_nm, mod_to_load)]
        
        effect_nm <- ifelse(
            effect_nm == "eff", "All sites",
                   "Outside protected areas - 'PA size' effect")
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
        
        diffs <- data.frame(Variable = c("Brodie", "Replicate"),
                            report_value = sprintf("%.1f%%", (c(diff_orig, diff))),
                            var_nm = var_nm, Effect = effect_nm) %>% 
            select(Effect, Variable, report_value, var_nm)
        
        rbind(vals, diffs)
        
        
    }) %>% bind_rows() %>% 
        mutate(var_nm = factor(var_nm, levels = c("SR", "FR", "PD"))) %>% 
        arrange(var_nm) %>% 
        pivot_wider(names_from = var_nm, values_from = report_value)
})
coefs <- left_join(coefs[[1]], coefs[[2]], by = c("Variable", "Effect"),
                   suffix = sprintf(".%s", c("Birds", "Mammals")))

coefs_efcy <- coefs %>% filter(Effect == "All sites") %>% select(-Effect)
coefs_size <- coefs %>% filter(Effect != "All sites") %>% select(-Effect)

## Make the kable
# Bold some cells
row_ids <- which(coefs_efcy$Variable %in% 
                     c("PA", "Connectivity",
                       "PA size (binary)"))
detec_sig <- coefs_efcy[row_ids, -1]
detec_sig <- do.call(cbind, lapply(1:ncol(detec_sig), function(i){
    as.numeric(sapply(detec_sig[[i]], function(x){
        str_extract(x, "(?<=;[[:space:]]).*?(?=\\))")}))
}))
detec_sig <- detec_sig <= 0.05
bold_ids <- matrix(data = FALSE, nrow = nrow(coefs_efcy), ncol = ncol(coefs_efcy))
bold_ids[row_ids, 1:ncol(detec_sig) + 1] <- detec_sig

for (i in 1:nrow(bold_ids)){
    for (j in 1:ncol(bold_ids)){
        if (bold_ids[i, j] == TRUE){
            coefs_efcy[i, j] <- cell_spec(coefs_efcy[i, j], bold = T, escape = FALSE)
        }
    }
}

names(coefs_efcy) <- c("Variable", "SR", "FR", "PD", "SR", "FR", "PD")
coefs_efcy %>%
    kableExtra::kbl(escape = FALSE, full_width = FALSE, booktabs = TRUE) %>% 
    kable_paper(html_font = "Helvetica") %>% 
    kable_styling(latex_options = c("scale_down", "hold_position"), 
                  font_size = 10) %>% 
    add_header_above(c(" " = 1, "Bird" = 3, "Mammal" = 3),
                     align = "c") %>% 
    save_kable("results/figures/coefs_efficacy.pdf")

# Bold some cells
row_ids <- which(coefs_size$Variable %in% 
                     c("PA", "Connectivity",
                       "PA size (binary)"))
detec_sig <- coefs_size[row_ids, -1]
detec_sig <- do.call(cbind, lapply(1:ncol(detec_sig), function(i){
    as.numeric(sapply(detec_sig[[i]], function(x){
        str_extract(x, "(?<=;[[:space:]]).*?(?=\\))")}))
}))
detec_sig <- detec_sig <= 0.05
bold_ids <- matrix(data = FALSE, nrow = nrow(coefs_size), ncol = ncol(coefs_size))
bold_ids[row_ids, 1:ncol(detec_sig) + 1] <- detec_sig

for (i in 1:nrow(bold_ids)){
    for (j in 1:ncol(bold_ids)){
        if (bold_ids[i, j] == TRUE){
            coefs_size[i, j] <- cell_spec(coefs_size[i, j], bold = T, escape = FALSE)
        }
    }
}
names(coefs_size) <- c("Variable", "SR", "FR", "PD", "SR", "FR", "PD")
coefs_size %>%
    kableExtra::kbl(escape = FALSE, full_width = FALSE, booktabs = TRUE) %>% 
    kable_paper(html_font = "Helvetica") %>% 
    kable_styling(latex_options = c("scale_down", "hold_position"), 
                  font_size = 10) %>% 
    add_header_above(c(" " = 1, "Bird" = 3, "Mammal" = 3),
                     align = "c") %>% 
    save_kable("results/figures/coefs_spillover_size.pdf")

compile_models(taxon = "bird", mod_type = "brodie",
               outliers = "brodie", med_dist = 10,
               src_dir = "data/derived/public", 
               dst_dir = "results")
compile_models(taxon = "mammal", mod_type = "brodie",
               outliers = "brodie", med_dist = 10,
               src_dir = "data/derived/public", 
               dst_dir = "results") 

load("results/models_bird_brodie_updated.rda")
models_bird <- models
load("results/models_mammal_brodie_updated.rda")
models_mammal <- models

# Concatenate the model results
mods <- list("bird" = models_bird, "mammal" = models_mammal)

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
    lapply(names(models)[1:6], function(nm){
        
        # Extract names for identity
        var_nm <- toupper(strsplit(nm, "_")[[1]][1])
        effect_nm <- strsplit(nm, "_")[[1]]
        effect_nm <- effect_nm[2]
        
        # Get the original model
        mod_to_load <- names_org[grep(tolower(var_nm), names_org)]
        mod_to_load <- mod_to_load[grep(effect_nm, mod_to_load)]
        
        effect_nm <- ifelse(
            effect_nm == "efficacy", "All sites",
            "Outside protected areas - 'PA size' effect")
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
        
        diffs <- data.frame(Variable = c("Brodie", "Replicate"),
                            report_value = sprintf("%.1f%%", (c(diff_orig, diff))),
                            var_nm = var_nm, Effect = effect_nm) %>% 
            select(Effect, Variable, report_value, var_nm)
        
        rbind(vals, diffs)
        
        
    }) %>% bind_rows() %>% 
        mutate(var_nm = factor(var_nm, levels = c("SR", "FR", "PD"))) %>% 
        arrange(var_nm) %>% 
        pivot_wider(names_from = var_nm, values_from = report_value)
})
coefs <- left_join(coefs[[1]], coefs[[2]], by = c("Variable", "Effect"),
                   suffix = sprintf(".%s", c("Birds", "Mammals")))

coefs_efcy <- coefs %>% dplyr::filter(Effect == "All sites") %>% select(-Effect)
coefs_size <- coefs %>% dplyr::filter(Effect != "All sites") %>% select(-Effect)

## Make the kable
# Bold some cells
row_ids <- which(coefs_efcy$Variable %in% 
                     c("PA", "Connectivity",
                       "PA size (binary)"))
detec_sig <- coefs_efcy[row_ids, -1]
detec_sig <- do.call(cbind, lapply(1:ncol(detec_sig), function(i){
    as.numeric(sapply(detec_sig[[i]], function(x){
        str_extract(x, "(?<=;[[:space:]]).*?(?=\\))")}))
}))
detec_sig <- detec_sig <= 0.05
bold_ids <- matrix(data = FALSE, nrow = nrow(coefs_efcy), ncol = ncol(coefs_efcy))
bold_ids[row_ids, 1:ncol(detec_sig) + 1] <- detec_sig

for (i in 1:nrow(bold_ids)){
    for (j in 1:ncol(bold_ids)){
        if (bold_ids[i, j] == TRUE){
            coefs_efcy[i, j] <- cell_spec(coefs_efcy[i, j], bold = T, escape = FALSE)
        }
    }
}

names(coefs_efcy) <- c("Variable", "SR", "FR", "PD", "SR", "FR", "PD")
coefs_efcy %>%
    kableExtra::kbl(escape = FALSE, full_width = FALSE, booktabs = TRUE) %>% 
    kable_paper(html_font = "Helvetica") %>% 
    kable_styling(latex_options = c("scale_down", "hold_position"), 
                  font_size = 10) %>% 
    add_header_above(c(" " = 1, "Bird" = 3, "Mammal" = 3),
                     align = "c") %>% 
    save_kable("results/figures/coefs_efficacy_brodie_updated.pdf")


# Bold some cells
row_ids <- which(coefs_size$Variable %in% 
                     c("PA", "Connectivity",
                       "PA size (binary)"))
detec_sig <- coefs_size[row_ids, -1]
detec_sig <- do.call(cbind, lapply(1:ncol(detec_sig), function(i){
    as.numeric(sapply(detec_sig[[i]], function(x){
        str_extract(x, "(?<=;[[:space:]]).*?(?=\\))")}))
}))
detec_sig <- detec_sig <= 0.05
bold_ids <- matrix(data = FALSE, nrow = nrow(coefs_size), ncol = ncol(coefs_size))
bold_ids[row_ids, 1:ncol(detec_sig) + 1] <- detec_sig

for (i in 1:nrow(bold_ids)){
    for (j in 1:ncol(bold_ids)){
        if (bold_ids[i, j] == TRUE){
            coefs_size[i, j] <- cell_spec(coefs_size[i, j], bold = T, escape = FALSE)
        }
    }
}
names(coefs_size) <- c("Variable", "SR", "FR", "PD", "SR", "FR", "PD")
coefs_size %>%
    kableExtra::kbl(escape = FALSE, full_width = FALSE, booktabs = TRUE) %>% 
    kable_paper(html_font = "Helvetica") %>% 
    kable_styling(latex_options = c("scale_down", "hold_position"), 
                  font_size = 10) %>% 
    add_header_above(c(" " = 1, "Bird" = 3, "Mammal" = 3),
                     align = "c") %>% 
    save_kable("results/figures/coefs_spillover_size_brodie_updated.pdf")
