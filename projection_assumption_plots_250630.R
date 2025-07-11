#
# Generates plots illustrating the projection assumptions for an OM
#
################################################################################

# Clear workspace
rm(list=ls())

# Load necessary libraries
library(MSEtool); library(SAMtool); library(tidyr); library(dplyr)
library(reshape2); library(ggplot2); library(patchwork)

# Choose the OM and set other important quantities
OM_file_name <- "M_decrease_in_M_F0.6.rds"
font_size <- 16 # Size of font for figures
fig.width <- 8 # Width of individual component figures
fig.height <- 6 # Height of individual component figures
fig.width.combined <- 6 # Width of the combined figure
fig.height.combined <- 10 # Height of the combined figure

# Grab the requested OM
OM_name <- sub("\\.(rds)$", "", OM_file_name)
OM <- readRDS(paste0("OM/",OM_file_name))

# Wt age plot
ymax <- max(OM@cpars$Wt_age[1,2:16,])*1.1
waa_plot <- melt(OM@cpars$Wt_age[1,2:16,], id.var = "id") %>% 
  ggplot(aes(x = Var2, y = value, col = as.factor(Var1))) + 
  geom_line(size = 1.5) + theme_classic() + xlab("Year") + 
  ylab("Weight (kg)") + theme(legend.position = "none", text = element_text(size = font_size)) +
  scale_color_viridis_d() + geom_vline(xintercept = 50.5, col = "red", lty = 2, size = 1.5) + 
  coord_cartesian(ylim = c(0, ymax), expand = F, xlim = c(0,105))
ggsave(paste0(getwd(),"/Figs/MSE/",OM_name,"/waa.jpeg"), plot=last_plot(), 
       width=fig.width, height=fig.height, units="in")

# M plot
ymax <- max(OM@cpars$M_ageArray[1,2:16,])*1.1
M_plot <- melt(OM@cpars$M_ageArray[1,2:16,], id.var = "id") %>%
  arrange(Var2) %>%
  group_by(Var1) %>%
  mutate(xend = lead(Var2), yend = value) %>%
  filter(!is.na(xend)) %>%
  ggplot() + 
  geom_segment(aes(x = Var2, xend = xend, y = value, yend = yend), 
               linetype = "solid", size = 1.5) +
  geom_vline(xintercept = 50.5, col = "red", lty = 2, size = 1.5) +
  coord_cartesian(ylim = c(0, ymax), expand = FALSE, xlim = c(0, 105)) +
  theme_classic() + xlab("Year") + ylab("M") +
  theme(legend.position = "none", text = element_text(size = font_size))
ggsave(paste0(getwd(),"/Figs/MSE/",OM_name,"/M.jpeg"), plot=last_plot(), 
       width=fig.width, height=fig.height, units="in")

# Recruit dev multiplier plot
ymax <- max(rep(OM@Misc$rec_timestep_factor, each = 25))*1.1
recdev_plot <- data.frame(year = 1:100, 
                          recdev = rep(OM@Misc$rec_timestep_factor, each = 25)) %>%
  arrange(year) %>%
  mutate(xend = lead(year), yend = recdev) %>%
  filter(!is.na(xend)) %>%
  ggplot() + 
  geom_segment(aes(x = year, xend = xend, y = recdev, yend = yend), 
               linetype = "solid", size = 1.5) +
  geom_vline(xintercept = 50.5, col = "red", lty = 2, size = 1.5) +
  coord_cartesian(ylim = c(0, ymax), expand = FALSE, xlim = c(0, 105)) +
  theme_classic() + xlab("Year") + ylab("Rec. dev. multiplier") +
  theme(legend.position = "none", text = element_text(size = font_size))
ggsave(paste0(getwd(),"/Figs/MSE/",OM_name,"/recdev.jpeg"), plot=last_plot(), 
       width=fig.width, height=fig.height, units="in")

# Combined plot of projection assumptions
M_plot + waa_plot + recdev_plot + plot_layout(ncol = 1, axes = "collect")
ggsave(paste0(getwd(),"/Figs/MSE/",OM_name,"/projections_assum.jpeg"), plot=last_plot(), 
       width=fig.width.combined, height=fig.height.combined, units="in")