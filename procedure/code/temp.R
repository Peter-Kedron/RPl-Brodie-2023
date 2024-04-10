library(ggplot2)
library(ggthemes)
library(patchwork)
library(tufte)

loadfonts(device=c("all"))
import_econ_sans()

intercept <- bird_pd_models$connec_100$coefficients$fixed["(Intercept)"]
B_conn <- bird_pd_models$connec_100$coefficients$fixed["awf_ptg.z"]
B_connPA <- B_conn + bird_pd_models$connec_100$coefficients$fixed["PA:awf_ptg.z"]

# Making a Economist inspired theme for the plot title
title_theme <- theme(
    plot.title = element_text(
        family = "Econ Sans Cnd", 
        face = "bold",
        size = 12
    ),
    plot.subtitle = element_text(
        family = "Econ Sans Cnd",
        size = 10,
    )
)

ggplot(data = dat_PD_efficacy) +
    geom_point(aes(x = awf_ptg.z, y = asymptPD, color = as.factor(PA)), size = 0.4) +
    ggsci::scale_color_npg(name = 'PA Status', labels = c("Outside PA", "Inside PA")) +
    geom_abline(aes(slope = B_conn, intercept = intercept, linetype = "B_conn")) +
    geom_abline(aes(slope = B_connPA, intercept = intercept, linetype = "B_connPA")) +
    scale_linetype_manual(name = "Coef.", values = c("B_conn" = 1, "B_connPA" = 2)) +
    geom_rangeframe() +
    labs( x = "Connectivity",
          y = "Phylogenetic Diversity",
          title = "Connectivity Moderates PA Efficacy",
          subtitle = "On phylogenetic diversity, at a dispersal distance of 100km") +
    theme_tufte() + title_theme

ggsave("/Users/pinot/downloads/PD_conn.png", width = 5, height = 4, bg = "white")
