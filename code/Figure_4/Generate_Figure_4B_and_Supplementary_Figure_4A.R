# Figure 4B and Supplementary Figure 4A

library(tidyverse)
library(ggplot2)
library(ggridges)
library(ggplot2)
library(dplyr)
library(forcats)
library(viridis)
library(tidygraph)
library(ggnetwork)
library(ggraph)
library(reshape2)
library(ggpubr)
library(hrbrthemes)


setwd("~/Dropbox/TungST_Morphology_Manuscript/Polar_plot_files/")

# Custom function to create polar coordinate system for ggplot
# This transforms a standard x-y plot into a circular/polar representation
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

# Change 'all_data_lh_byROI' to 'all_data_rh_byROI' to create Supplementary Figure 4A (right hemisphere)
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

group_col<-c("#B5D5F0","#7EACD1","#2C6A9F","#000000")

# Create the polar plot
finaldata.plot <- finaldata %>% 
  ggplot() + 
  geom_polygon(
    aes(x = variable, y = mean, group = AgeCat, color = AgeCat, fill = AgeCat), 
    size = 4, 
    alpha = .01                          
  ) + 
  scale_color_brewer(palette = "Set3") +
  geom_point(
    aes(x = variable, y = mean, fill = NA), 
    size = 2, 
    shape = 2                            
  ) + 
  coord_radar() + 
  scale_y_continuous(
    labels = paste(seq(0, 90, by = 10), ""),  
    breaks = seq(0, 0.9, by = 0.1),  
    limits = c(0, 0.90)
  ) +
  labs(
    x = "Metric",
    y = "Units", 
    color = "Label",
    fill = "Label"
  ) + 
  theme_minimal() + 
  theme(
    plot.title = element_text(family = "Arial", color = "black", size = 12, hjust = 0.5), 
    axis.title = element_text(family = "Arial", color = "black", face = "italic", size = 0),  
    axis.text.y = element_text(family = "Arial", color = "black", size = 7),
    axis.text.x = element_text(family = "Arial", color = "black", size = 0),  # Hide variable names
    legend.title = element_text(family = "Arial", color = "black", size = 14), 
    legend.text = element_text(family = "Arial", color = "black", size = 12), 
    strip.text = element_text(family = "Arial", color = "black", face = "italic", size = 12, hjust = 0.5), 
    panel.grid.minor = element_line(size = .7, color = "black"), 
    panel.grid.major = element_line(size = .7, color = "black")
  ) 

# Apply the custom color palette for age groups and display the final plot
finaldata.plot + scale_color_manual(values = group_col)