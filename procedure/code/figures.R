library(ggplot2)
library(ggthemes)
library(patchwork)
library(tufte)

intercept <- coef(bird_pd_models$connec_100$coefficients$fixed["(Intercept)"])
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

ggplot(data = dat_PD_efficacy, aes(x = awf_ptg.z, y = asymptPD, color = PA)) +
  geom_point() + 
  xlim(-0.2, 5) +
  geom_rangeframe() +
  labs( x = "Connectivity",
        y = "Biodiversity (PD)"
       ) +
  plot_annotation(title = "Moderating Effect", 
                  subtitle = "Connectivity Moderates PA Efficacy", 
                  theme = title_theme) +
  theme_tufte() +
  guides(col = "none")