library(tidyverse)
library(ggplot2)
library(ggridges)
library(dplyr)
library(forcats)
library(viridis)
library(tidygraph)
library(ggnetwork)
library(ggraph)
library(reshape2)
library(ggpubr)
library(colorspace) 

# Custom coordinate function
coord_radar <- function(theta = "x", start = 0, direction = 1) {
  theta <- match.arg(theta, c("x", "y"))
  r <- if (theta == "x") "y" else "x"
  ggproto(
    "CordRadar", CoordPolar,
    theta = theta, r = r,
    start = start, direction = sign(direction),
    is_linear = function(coord) TRUE
  )
}

# Normalization
normalit <- function(m) {
  (m - min(m)) / (max(m) - min(m))
}

# Load and normalize your data
df2 <- all_data_rh_byROI %>%
  summarise(
    SulcalDepth = normalit(SulcalDepth),
    SulcalSpan  = normalit(SurfaceArea),
    Curvature   = normalit(Curvature),
    Thickness   = normalit(Thickness),
    R1          = normalit(R1),
    AgeCat,
    ROI
  )

df_melt <- melt(df2, id = c("ROI", "AgeCat"))

# ROI colors
roi_base_cols <- c(
  "#b3eafc", "#9abdf9", "#548af9", "#2e5cd9", "#b1df85",
  "#66c3a2", "#43969d", "#265c67", "#fae04b", "#f1a039",
  "#c58038", "#8f5f22", "#f3b1e7", "#ea3394", "#901f5a"
)

# Map colors to ROIs
roi_names <- unique(df_melt$ROI)
base_map <- setNames(rep(roi_base_cols, length.out = length(roi_names)), roi_names)

# Make shades for age groups
make_shades <- function(base_col, n = 4) {
  if (n == 4) {
    # 0 represents black
    # 1 represents base color 
    mix_ratio <- c(0, 0.5, 0.8, 1)  
  } else {
    mix_ratio <- seq(0, 1, length.out = n)
  }
  
  shades <- sapply(mix_ratio, function(r) {
    mix <- colorspace::mixcolor(r,
                                colorspace::hex2RGB("#000000"),  
                                colorspace::hex2RGB(base_col))   
    colorspace::hex(mix)
  })
  # Reverse the order so black goes to 4m and base color goes to 1m
  return(rev(shades))
}


plot_list <- list()

for (roi in roi_names) {
  cat("Processing ROI:", roi, "\n")
  
  singleroi <- df_melt[df_melt$ROI == roi, ]
  singleroi$AgeCat <- factor(singleroi$AgeCat, levels = unique(singleroi$AgeCat))
  
  finaldata <- singleroi %>%
    group_by(variable, AgeCat) %>%
    summarise(
      n = n(),
      mean = mean(value),
      sd = sd(value),
      se = sd / sqrt(n),
      .groups = "drop"
    )
  
  n_groups <- length(levels(finaldata$AgeCat))
  roi_shades <- make_shades(base_map[[roi]], n = n_groups)
  names(roi_shades) <- levels(finaldata$AgeCat)

  
  current_plot <- finaldata %>%
    ggplot() +
    geom_ribbon(
      aes(
        x = variable,
        ymin = pmax(0, mean - se),
        ymax = mean + se,
        group = AgeCat,
        fill = AgeCat
      ),
      alpha = 0.4
    ) +
    geom_polygon(
      aes(x = variable, y = mean, group = AgeCat, color = AgeCat),
      size = 6,
      fill = NA
    ) +
    geom_point(
      aes(x = variable, y = mean, color = AgeCat),
      size = 1,
      shape = 16
    ) +
    coord_radar() +
    scale_y_continuous(
      labels = paste(seq(0, 90, by = 10), ""),
      breaks = seq(0, 0.9, by = 0.1),
      limits = c(0, 0.9)
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(
        family = "Arial", color = "gray", size = 16,
        hjust = 0.5, face = "bold"
      ),
      axis.title = element_blank(),
      axis.text.y = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks  = element_blank(),
      legend.title = element_text(family = "Arial", color = "black", size = 14),
      legend.text  = element_text(family = "Arial", color = "black", size = 12),
      legend.position = "none",
      panel.grid.minor = element_line(size = 0.5, color = "#D6D6D6"),
      panel.grid.major = element_line(size = 0.5, color = "#D6D6D6")
    ) +
    scale_color_manual(values = roi_shades) +
    scale_fill_manual(values  = roi_shades)
  
  plot_list[[roi]] <- current_plot

}
