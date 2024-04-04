library(ggplot2)
library(ggthemes)
library(tufte)

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
  labs( x = "Connectivity",
        y = "Biodiversity (PD)"
  )