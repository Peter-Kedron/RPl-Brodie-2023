pkgs <- c("tidyverse", "hrbrthemes", "viridis", "plotly", "stats",
          "scales", "ggplot2", "here")
lapply(pkgs, library, character.only=TRUE) # use this line if groundhog returns error
here('')


saveBoxPlot <- function(model, modtype, exp1, exp2, y1, y2, yadjust,
                        ylabel, outpath){
    # model: the lme model
    # exp1: the expression that goes with 0, e.g., unprotected areas
    # exp2: the other expression
    # outpath: output path, including the name of the plot
    df <- as.data.frame(model$fitted)
    df$fitted <- fitted(model)
    
    df$iftrtmt <- model$data[modtype]
    df$orig <- model$data$y
    
    df <- df %>% mutate(
        trtmt = ifelse(iftrtmt==0, exp1, exp2)
    )
    df$trtmt <- factor(df$trtmt, levels = c(exp1, exp2))
    
    df1 <- df[df$trtmt==exp1, ]
    df2 <- df[df$trtmt==exp2, ]
    
    l1 <- df[df$trtmt==exp1, ]$fitted
    l2 <- df[df$trtmt==exp2, ]$fitted
    
    pval <- t.test(l1, l2)$p.value
    if(pval < 0.001){
        str_pval = "< 0.001"
    } else{
        str_pval = paste(" =", as.character(round(pval, digits = 3)))
    }
    
    #dltxb <- 100*abs((mean(l1)-mean(l2))*2/(mean(l1)+mean(l2)))
    dltxb <- 100*abs(mean(l2)-mean(l1))/mean(l2)
    str_dltxb <- paste(" = ", as.character(round(dltxb, digits = 1)), "%", sep = "")
    
    lab_pval <- bquote(italic("P") ~ .(str_pval))
    lab_dltxb <- bquote(Delta*italic(bar(x)) ~ .(str_dltxb))
    
    p <- df %>%
        ggplot(aes(x=trtmt, y=fixed, fill=trtmt, group=trtmt)) +
        geom_boxplot(outlier.shape = NA, fatten = NULL, aes(color = trtmt),
                     width=0.3) +
        theme(
            legend.position="none",
            #plot.title = element_text(size=11),
            panel.border = element_blank(),
            axis.line.x = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
            axis.line.y = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
            axis.text.x = element_text(family = "Times New Roman", size = 11, colour = "black"),
            axis.text.y = element_text(family = "Times New Roman", size = 11, colour = "black"),
            axis.title.y.left = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0),
                                             family = 'Times New Roman', size = 11),
            plot.margin = margin(1, 2, 1, 2, "cm"),
            panel.grid.major.x = element_blank(),
            panel.grid.major.y = element_blank(),
            panel.grid.minor.x = element_blank(),
            panel.grid.minor.y = element_blank(),
            axis.line = element_line(size = 0.5),
            axis.ticks = element_line(size = 0.5),
            panel.background = element_rect(colour = 'white', fill = 'white')
        ) +
        scale_color_manual(values=c("#ff702f", "#67c2f3"))+
        scale_fill_manual(values=c("#ff702f", "#67c2f3")) +
        xlab("")+
        ylab(ylabel)+
        geom_text(label = deparse(lab_pval), x = 1, y = y1, size = 4, family='Times New Roman', parse = TRUE)+
        geom_text(label = deparse(lab_dltxb), x = 1, y = y2, size = 4, parse=TRUE, family='Times New Roman') #+
    
    if(length(yadjust==2)){
        p <- p + ylim(yadjust[1], yadjust[2])
    }
        #ylim(1, 3) # only when needed
    p
    
    ggsave(outpath, width = 5, height = 5)
}

# ----------------------- load all models -------------------------------
rm(list = ls())
bird_mods_path <- 'data/derived/models_bird_100.rda'
mammal_mods_path <- 'data/derived/models_mammal_100.rda'
bird_brodie_path <- 'data/derived/models_bird_brodie_orig.rda'
mammal_brodie_path <- 'data/derived/models_mammal_brodie_orig.rda'

load(bird_mods_path)
bird_models <- models
rm(models)

load(mammal_mods_path)
mammal_models <- models
rm(models)


load(bird_brodie_path)
brodie_bird_models <- models
rm(models)

load(mammal_brodie_path)
brodie_mammal_models <- models
rm(models)
# ----------------------- brodie efficacy models -------------------------------
eff_mod <- "PA"
eff1 <- "Unprotected areas"
eff2 <- "PAs"

0.90*max(as.data.frame(bird_models$pd_efficacy_brodie$fitted)$fixed) # 0.94, 0.92
saveBoxPlot(bird_models$pd_efficacy_brodie, eff_mod, eff1, eff2, 
            4.2, 3.92, c(1, 4.5), "PD",
            "results/figures/boxplot/brodie_birds_pd.png")
saveBoxPlot(bird_models$fr_efficacy_brodie, eff_mod, eff1, eff2, 
            365, 340, 0, "FR",
            "results/figures/boxplot/brodie_birds_fr.png")
saveBoxPlot(bird_models$sr_efficacy_brodie, eff_mod, eff1, eff2, 
            215, 200, 0, "SR",
            "results/figures/boxplot/brodie_birds_sr.png")

saveBoxPlot(mammal_models$pd_efficacy_brodie, eff_mod, eff1, eff2, 
            3.3, 3.1, c(1, 3.5), "PD",
            "results/figures/boxplot/brodie_mammals_pd.png")
saveBoxPlot(mammal_models$fr_efficacy_brodie, eff_mod, eff1, eff2,
            18.5, 17.2, c(2, 20), "FR",
            "results/figures/boxplot/brodie_mammals_fr.png")
saveBoxPlot(mammal_models$sr_efficacy_brodie, eff_mod, eff1, eff2, 
            12.2, 11.8, c(7.5, 12.7), "SR",
            "results/figures/boxplot/brodie_mammals_sr.png")

# ----------------------- connec efficacy models -------------------------------
0.90*max(as.data.frame(bird_models$sr_efficacy_connec$fitted)$fixed) # 0.94, 0.92
saveBoxPlot(bird_models$pd_efficacy_connec, eff_mod, eff1, eff2, 
            4.0, 3.85, c(2, 4.2), "PD",
            "results/figures/boxplot/connec_birds_pd.png")
saveBoxPlot(bird_models$fr_efficacy_connec, eff_mod, eff1, eff2, 
            355, 335, c(50, 380), "FR",
            "results/figures/boxplot/connec_birds_fr.png")
saveBoxPlot(bird_models$sr_efficacy_connec, eff_mod, eff1, eff2, 
            205, 193, c(30, 220), "SR",
            "results/figures/boxplot/connec_birds_sr.png")

0.90*max(as.data.frame(mammal_models$sr_efficacy_connec$fitted)$fixed) # 0.94, 0.92
saveBoxPlot(mammal_models$pd_efficacy_connec, eff_mod, eff1, eff2, 
            3.1, 2.9, c(0.5, 3.5), "PD",
            "results/figures/boxplot/connec_mammals_pd.png")
saveBoxPlot(mammal_models$fr_efficacy_connec, eff_mod, eff1, eff2,
            19, 17.8, 0, "FR",
            "results/figures/boxplot/connec_mammals_fr.png")
saveBoxPlot(mammal_models$sr_efficacy_connec, eff_mod, eff1, eff2, 
            13.7, 12.9, c(3, 15), "SR",
            "results/figures/boxplot/connec_mammals_sr.png")

# --------------------------- spillover models ---------------------------------
# ------------------------------ size models -----------------------------------
# bird_size_models <- list(PD_size = mm02b_birds_size_pd, FR_size = mm02b_birds_size_fr, SR_size = mm02b_birds_size_sr)
size_spill <- "BigPA"
size1 <- "Near small PA"
size2 <- "Near large PA"
# brodie
0.90*max(as.data.frame(bird_models$pd_size_spillover_brodie$fitted)$fixed) # 0.94, 0.92
saveBoxPlot(bird_models$pd_size_spillover_brodie, size_spill, size1, size2, 
            4.7, 4.5, c(2, 5), "PD",
            "results/figures/boxplot/brodie_size_birds_pd.png")
saveBoxPlot(bird_models$fr_size_spillover_brodie, size_spill, size1, size2,
            330, 305, 0, "FR",
            "results/figures/boxplot/brodie_size_birds_fr.png")
saveBoxPlot(bird_models$sr_size_spillover_brodie, size_spill, size1, size2,
            235, 220, c(25, 260), "SR",
            "results/figures/boxplot/brodie_size_birds_sr.png")

# connec
0.90*max(as.data.frame(bird_models$sr_size_spillover_connec$fitted)$fixed) # 0.94, 0.92
saveBoxPlot(bird_models$pd_size_spillover_connec, size_spill, size1, size2, 
            4.2, 4, c(2, 4.5), "PD",
            "results/figures/boxplot/connec_size_birds_pd.png")
saveBoxPlot(bird_models$fr_size_spillover_connec, size_spill, size1, size2,
            300, 280, c(100, 320), "FR",
            "results/figures/boxplot/connec_size_birds_fr.png")
saveBoxPlot(bird_models$sr_size_spillover_connec, size_spill, size1, size2,
            210, 190, 0, "SR",
            "results/figures/boxplot/connec_size_birds_sr.png")

# mammal size models
# mammals_size_models <- list(PD_size = mm02b_mammals_size_pd, FR_size = mm02b_mammals_size_fr, SR_size = mm02b_mammals_size_sr)
# brodie
0.90*max(as.data.frame(mammal_models$pd_size_spillover_brodie$fitted)$fixed) # 0.94, 0.92
saveBoxPlot(mammal_models$pd_size_spillover_brodie, size_spill, size1, size2, 
            3.7, 3.5, c(1,4), "PD",
            "results/figures/boxplot/brodie_size_mammals_pd.png")
saveBoxPlot(mammal_models$fr_size_spillover_brodie, size_spill, size1, size2, 
            25, 22.5, c(-8, 28), "FR",
            "results/figures/boxplot/brodie_size_mammals_fr.png")
saveBoxPlot(mammal_models$sr_size_spillover_brodie, size_spill, size1, size2,
            12.3, 11.7, c(5, 12.8), "SR",
            "results/figures/boxplot/brodie_size_mammals_sr.png")

# connec
0.90*max(as.data.frame(mammal_models$sr_size_spillover_connec$fitted)$fixed) # 0.94, 0.92
saveBoxPlot(mammal_models$pd_size_spillover_connec, size_spill, size1, size2, 
            2.8, 2.65, c(1, 3), "PD",
            "results/figures/boxplot/connec_size_mammals_pd.png")
saveBoxPlot(mammal_models$fr_size_spillover_connec, size_spill, size1, size2, 
            18.2, 16.9, 0, "FR",
            "results/figures/boxplot/connec_size_mammals_fr.png")
saveBoxPlot(mammal_models$sr_size_spillover_connec, size_spill, size1, size2,
            12, 11.1, 0, "SR",
            "results/figures/boxplot/connec_size_mammals_sr.png")

# ------------------------------ dist models -----------------------------------
# bird dist models
#bird_dist_models <- list(PD_dist = mm02d_birds_dist_pd, FR_dist = mm02d_birds_dist_fr, SR_dist = mm02d_birds_dist_sr)
dist_spill <- "CloseToPA"
dist1 <- "Far from PA"
dist2 <- "Close to PA"

# brodie
0.90*max(as.data.frame(bird_models$pd_efficacy_brodie$fitted)$fixed) # 0.94, 0.92
saveBoxPlot(bird_models$pd_dist_spillover_brodie, dist_spill, dist1, dist2,
            4.2, 4, c(1.5, 4.5), "PD",
            "results/figures/boxplot/brodie_dist_birds_pd.png")
saveBoxPlot(bird_models$fr_dist_spillover_brodie, dist_spill, dist1, dist2,
            350, 328, 0, "FR",
            "results/figures/boxplot/brodie_dist_birds_fr.png")
saveBoxPlot(bird_models$sr_dist_spillover_brodie, dist_spill, dist1, dist2,
            225, 210, 0, "SR",
            "results/figures/boxplot/brodie_dist_birds_sr.png")

# connec
0.90*max(as.data.frame(bird_models$sr_efficacy_connec$fitted)$fixed) # 0.94, 0.92
saveBoxPlot(bird_models$pd_dist_spillover_connec, dist_spill, dist1, dist2,
            4.2, 4, c(1.5, 4.5), "PD",
            "results/figures/boxplot/connec_dist_birds_pd.png")
saveBoxPlot(bird_models$fr_dist_spillover_connec, dist_spill, dist1, dist2,
            340, 310, 0, "FR",
            "results/figures/boxplot/connec_dist_birds_fr.png")
saveBoxPlot(bird_models$sr_dist_spillover_connec, dist_spill, dist1, dist2,
            210, 195, 0, "SR",
            "results/figures/boxplot/connec_dist_birds_sr.png")

# mammal dist models
#mammals_dist_models <- list(PD_dist = mm02d_mammals_dist_pd, FR_dist = mm02d_mammals_dist_fr, SR_dist = mm02d_mammals_dist_sr)

# brodie
0.90*max(as.data.frame(mammal_models$pd_dist_spillover_brodie$fitted)$fixed) # 0.94, 0.92
saveBoxPlot(mammal_models$pd_dist_spillover_brodie, dist_spill, dist1, dist2,
            3.8, 3.6, c(1, 4), "PD",
            "results/figures/boxplot/brodie_dist_mammals_pd.png")
saveBoxPlot(mammal_models$fr_dist_spillover_brodie, dist_spill, dist1, dist2,
            20.3, 19.1, c(3, 22), "FR",
            "results/figures/boxplot/brodie_dist_mammals_fr.png")
saveBoxPlot(mammal_models$sr_dist_spillover_brodie, dist_spill, dist1, dist2,
            18.7, 17.5, c(2, 20), "SR",
            "results/figures/boxplot/brodie_dist_mammals_sr.png")

# connec
0.90*max(as.data.frame(mammal_models$sr_dist_spillover_connec$fitted)$fixed) # 0.94, 0.92
saveBoxPlot(mammal_models$pd_dist_spillover_connec, dist_spill, dist1, dist2,
            2.85, 2.7, c(1, 3), "PD",
            "results/figures/boxplot/connec_dist_mammals_pd.png")
saveBoxPlot(mammal_models$fr_dist_spillover_connec, dist_spill, dist1, dist2,
            23, 21.5, c(3,25), "FR",
            "results/figures/boxplot/connec_dist_mammals_fr.png")
saveBoxPlot(mammal_models$sr_dist_spillover_connec, dist_spill, dist1, dist2,
            15.6, 14.6, c(4, 16.5), "SR",
            "results/figures/boxplot/connec_dist_mammals_sr.png")
