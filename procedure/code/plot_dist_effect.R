## -------------------------------------------------------------------
## Script name: plot_dist_effect
## Purpose of script: Plot the change in coefficient of interested variables
## (PA, connectivity, and interaction) with changing dispersal distance to 
## calculate connectivity.
## Author: Lei Song
## Date Created: 2024-08-26
## Email: lsong@ucsb.edu

## Inputs:
## src_dir (character): The directory of models
## dst_dir (character): The directory to save out the plots.

## Outputs:
## Save out a list of gg object. The structure is:
## taxon_dist_plot:
##  - bird
##      - Phylogenetic diversity (PD)
##      - Functional richness (FR)
##      - Species richness (SR)
##  - mammal
##      - Phylogenetic diversity (PD)
##      - Functional richness (FR)
##      - Species richness (SR)
## -------------------------------------------------------------------

plot_dist_effect <- function(src_dir = "results",
                             dst_dir = "results/figures"){
    # Load associated models
    fnames <- list.files(src_dir, pattern = "[0-9]+.rda$")
    taxons <- unique(str_extract(fnames, "mammal|bird")) # Could just hard code
    
    # Tidy up the models
    coefs <- lapply(taxons, function(taxon){
        # Subset fnames
        fnames <- fnames[str_detect(fnames, taxon)]
        fnames <- fnames[str_detect(fnames, "[0-9]+")]
        
        # Tidy up coefficients of models
        lapply(fnames, function(fname){
            # Load the model object
            load(file.path(src_dir, fname))
            
            # Select relevant models
            nms <- names(models)
            nms <- nms[!str_detect(nms, "brodie")]
            
            # Per model
            lapply(nms, function(nm){
                # Extract names for identity
                var_nm <- strsplit(nm, "_")[[1]][1]
                effect_nm <- strsplit(nm, "_")[[1]][2]
                effect_nm <- ifelse(effect_nm == "efficacy", "efficacy",
                                    paste(effect_nm, "spillover", sep = "_"))
                
                # Do the tidy work
                broom.mixed::tidy(
                    models[[nm]], effects='fixed', conf.int = TRUE) %>% 
                    mutate(var_nm = var_nm, effect_nm = effect_nm)
            }) %>% bind_rows() %>% 
                mutate(dist = as.integer(str_extract(fname, "[0-9]+")))
        }) %>% bind_rows() %>% mutate(taxon = taxon)
    }) %>% bind_rows()
    
    # Make figures
    taxon_dist_plot <- lapply(taxons, function(txn){
        plot_list <- lapply(c("pd", "fr", "sr"), function(rsp){
            # Subset 
            dat <- coefs %>% dplyr::filter(taxon == txn) %>% 
                dplyr::filter(var_nm == rsp) %>% 
                dplyr::filter(term %in% c("PA", "connectivity.z",
                                "BigPA", "CloseToPA")) %>% 
                mutate(effect_nm = factor(
                    effect_nm, levels = c("efficacy", "size_spillover", "dist_spillover"),
                    labels = c("bold(PA~efficacy)", "bold(PA~size~effect)", "bold(Distance~to~PA~effect)"))) %>% 
                mutate(term = ifelse(term == "connectivity.z", "Connectivity", term)) %>% 
                mutate(term = factor(
                    term, levels = c("PA", "BigPA", "CloseToPA", "Connectivity"),
                    labels = c("PA", "PA~Size", "Distance~to~PA", "Connectivity")))
            
            # Plot
            ggplot(dat, 
                   aes(x = dist, y = estimate)) +
                geom_pointrange(
                    size = 0.2,
                    aes(ymin = conf.low, ymax = conf.high),
                    position = position_dodge(width = 2)) +
                geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
                facet_nested_wrap(.~effect_nm + term, scales = "free",
                                  nest_line = element_line(colour = "darkgrey",
                                                           linetype = "dotdash"),
                                  labeller = label_parsed, ncol = 2) +
                scale_x_continuous(breaks = seq(10, 150, by = 30))+
                xlab('Dispersal distance (km)') +
                ylab('Eestimated effect') +
                theme_tufte()+
                theme(panel.border = element_rect(colour = 'white', fill = NA),
                      panel.spacing.y = unit(2, "lines"),
                      strip.text = element_text(size = 12))
        })
        names(plot_list) <- c("Phylogenetic diveristy (PD)", 
                              "Functional richness (FR)",
                              "Species richness (SR)")
        plot_list
    })
    
    names(taxon_dist_plot) <- taxons
    
    save(taxon_dist_plot, file = file.path(dst_dir, "plot_dist_effect_taxons.rda"))
}
