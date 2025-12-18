# Figure 4B and Supplementary Figure 4A

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

# Custom function to create polar coordinate system for ggplot
# Acknowledgements to Ethan Willbrand for their help on this script

coord_radar <- function(theta = "x", start = 0, direction = 1) { 
  theta <- match.arg(theta, c("x", "y"))
  r <- if (theta == "x") "y" else "x"
  ggproto("CordRadar", CoordPolar, theta = theta, r = r, 
          start = start, direction = sign(direction), 
          is_linear = function(coord) TRUE)
}

# Function to normalize values to a [0,1] range
normalit <- function(m) { 
  (m - min(m)) / (max(m) - min(m))
}

# Change to 'all_data_rh_byROI' to create Supplementary Figure 4A (right hemisphere)
df2 <- all_data_lh_byROI %>% 
  summarise(
    SulcalDepth = normalit(SulcalDepth),   
    SulcalSpan = normalit(SulcalWidth),   
    Curvature = normalit(Curvature),       
    R1_gray = normalit(R1_gray),           
    Thickness = normalit(Thickness),       
    AgeCat,                                
    ROI                                    
  ) 

df_melt <- melt(df2, id = c("ROI", "AgeCat"))

# Change ROI name to analyze different brain regions
# Available ROIs: calcarine, pos, insula, central, cos, sts, sfs, ips, los, ifs, ots, its
singleroi <- df_melt[df_melt$ROI == "calcarine", ]

# Calculate summary statistics 
finaldata <- singleroi %>% 
  group_by(variable, AgeCat) %>% 
  summarise(
    n = n(),                              
    mean = mean(value),                   
    sd = sd(value),                       
    se = sd/sqrt(n)                       
  )

# Age group colors
group_col <- c("#B5D5F0", "#7EACD1", "#2C6A9F", "#000000")

# Create the polar plot with error bands
finaldata.plot <- finaldata %>% 
  ggplot() + 
  geom_ribbon(
    aes(x = variable, 
        ymin = pmax(0, mean - se),  
        ymax = mean + se,         
        group = AgeCat, 
        fill = AgeCat), 
    alpha = 0.4                    
  ) +
  geom_polygon(
    aes(x = variable, y = mean, group = AgeCat, color = AgeCat), 
    size = 0.6, 
    fill = NA                       
  ) + 
  geom_point(
    aes(x = variable, y = mean, color = AgeCat), 
    size = 1, 
    shape = 16
  ) + 
  # Convert to radar/polar coordinates
  coord_radar() + 
  scale_y_continuous(
    labels = paste(seq(0, 90, by = 10), ""),  
    breaks = seq(0, 0.9, by = 0.1),  
    limits = c(0, 0.90)
  ) +
  # Add labels
  labs(
    x = "Metric",
    y = "Units", 
    color = "Age Group",
    fill = "Age Group"
  ) + 
  theme_minimal() + 
  theme(
    plot.title = element_text(family = "Arial", color = "black", size = 12, hjust = 0.5), 
    axis.title = element_blank(),  
    axis.text.y = element_blank(),  
    axis.text.x = element_blank(),  
    axis.ticks = element_blank(), 
    legend.title = element_text(family = "Arial", color = "black", size = 14), 
    legend.text = element_text(family = "Arial", color = "black", size = 12), 
    strip.text = element_text(family = "Arial", color = "black", face = "italic", size = 12, hjust = 0.5), 
    panel.grid.minor = element_line(size = .1, color = "black"), 
    panel.grid.major = element_line(size = .1, color = "black")
  ) 

finaldata.plot + 
  scale_color_manual(values = group_col) +
  scale_fill_manual(values = group_col)  
