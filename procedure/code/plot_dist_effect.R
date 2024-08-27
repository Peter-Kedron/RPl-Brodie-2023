## -------------------------------------------------------------------
## Script name: plot_dist_effect
## Purpose of script: Plot the change in coefficient of interested variables
## (PA, connectivity, and interaction) with changing dispersal distance to 
## calculate connectivty.
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
##      - Phylogenetic diveristy (PD)
##      - Functional richness (FR)
##      - Species richness (SR)
##  - mammal
##      - Phylogenetic diveristy (PD)
##      - Functional richness (FR)
##      - Species richness (SR)
## -------------------------------------------------------------------

plot_dist_effect <- function(src_dir = "results",
                             dst_dir = "results/figures"){
    # Load associated models
    fnames <- list.files(src_dir, pattern = ".rda$")
    taxons <- unique(str_extract(fnames, "mammal|bird")) # Could just hardcode
    
    # Tidy up the models
    coefs <- lapply(taxons, function(taxon){
        # Subset fames
        fnames <- fnames[str_detect(fnames, taxon)]
        fnames <- fnames[str_detect(fnames, "[0-9]+")]
        
        # Tidy up coefficients of models
        lapply(fnames, function(fname){
            # Load the model object
            load(file.path(src_dir, fname))
            
            # Selete relevant models
            nms <- names(models)
            nms <- nms[!str_detect(nms, "brodie")]
            nms <- nms[!str_detect(nms, "spillover_connec$")]
            
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
            dat <- coefs %>% filter(taxon == txn) %>% 
                filter(var_nm == rsp) %>% 
                filter(term %in% c("PA", "connectivity.z", "PA:connectivity.z",
                                "BigPA", "BigPA:connectivity.z", "CloseToPA",
                                "CloseToPA:connectivity.z")) %>% 
                mutate(effect_nm = factor(
                    effect_nm, levels = c("efficacy", "size_spillover", "dist_spillover"),
                    labels = c("PA efficacy", "PA size effect", "Distance to PA effect"))) %>% 
                mutate(term = ifelse(term == "connectivity.z", "Connectivity",
                                     ifelse(str_detect(term, ":connectivity.z"), "Interaction", term))) %>% 
                mutate(term = factor(
                    term, levels = c("PA", "BigPA", "CloseToPA", "Connectivity", "Interaction"),
                    labels = c("PA", "PA Size", "Distance to PA", "Connectivity", "Interaction")))
            
            # Plot 
            ggplot(dat, 
                   aes(x = dist, y = estimate)) +
                geom_pointrange(
                    size = 0.2,
                    aes(ymin = conf.low, ymax = conf.high),
                    position = position_dodge(width = 2)) +
                facet_wrap(.~effect_nm + term, scales = "free_y") + 
                scale_x_continuous(breaks = seq(10, 150, by = 10))+
                xlab('Dispersal distance (km)') +
                ylab('Eestimated effect') +
                # geom_hline(yintercept=0.38, linetype='dashed', color = "red")+
                # ylab(paste(var, 'effect', sep=' '))+
                # facet_wrap(~ind+varname, scales = 'free_y')+
                theme_tufte()+
                theme(
                    #strip.text.x = element_text(angle = 0, vjust = 0),
                    panel.border = element_rect(colour='white', fill=NA))
        })
        names(plot_list) <- c("Phylogenetic diveristy (PD)", 
                              "Functional richness (FR)",
                              "Species richness (SR)")
        plot_list
    })
    
    names(taxon_dist_plot) <- taxons
    
    save(taxon_dist_plot, file = file.path(dst_dir, "plot_dist_effect_taxons.rda"))
}