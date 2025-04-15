#
# Code to combine the MSE plots into panel figures
#
################################################################################

# Load necessary libraries
library(patchwork); library(ggplot2)

# Set the font size and figure sizes
font_size <- 20
fig_width_c <- 14; fig_height_c <- 16 # For the SSB/SSB_MSY figures
fig_width_l <- 11.2; fig_height_l <- 19.2 # For the LRP/USR figures
fig_width_r <- 16; fig_height_r <- 10 # For the radar figures
fig_width_f <- 16; fig_height_f <- 10 # For the full time series figures
fig_width <- 16; fig_height <- 16 # For all other figures

# Load in plots
growthdown_lowF_plots <- readRDS(paste0(getwd(), "/Figs/MSE/M_decreasing_waa_rec_F0.3/M_decreasing_waa_rec_F0.3plots.rds"))
growthdown_highF_plots <- readRDS(paste0(getwd(), "/Figs/MSE/M_decreasing_waa_rec_F0.4/M_decreasing_waa_rec_F0.4plots.rds"))
growthup_lowF_plots <- readRDS(paste0(getwd(), "/Figs/MSE/M_increasing_waa_rec_F0.4/M_increasing_waa_rec_F0.4plots.rds"))
growthup_highF_plots <- readRDS(paste0(getwd(), "/Figs/MSE/M_increasing_waa_rec_F0.55/M_increasing_waa_rec_F0.55plots.rds"))
decreaseM_lowF_plots <- readRDS(paste0(getwd(), "/Figs/MSE/M_decrease_in_M_F0.55/M_decrease_in_M_F0.55plots.rds"))
decreaseM_highF_plots <- readRDS(paste0(getwd(), "/Figs/MSE/M_decrease_in_M_F0.6/M_decrease_in_M_F0.6plots.rds"))
increaseM_lowF_plots <- readRDS(paste0(getwd(), "/Figs/MSE/M_increase_in_M_F0.2/M_increase_in_M_F0.2plots.rds"))
increaseM_highF_plots <- readRDS(paste0(getwd(), "/Figs/MSE/M_increase_in_M_F0.3/M_increase_in_M_F0.3plots.rds"))

# Set up themes to call on for individual panels
theme_no_xy <- theme(axis.title = element_blank(), axis.ticks.x = element_blank(),
                    axis.text.x = element_blank(), axis.ticks.y = element_blank(),
                    axis.text.y = element_blank())
theme_no_x <- theme(axis.title.x = element_blank(), axis.ticks.x = element_blank(),
                    axis.text.x = element_blank())
theme_no_y <- theme(axis.title.y = element_blank(), 
                    axis.ticks.y = element_blank(), axis.text.y = element_blank(), 
                    axis.text.x = element_text(angle = 90, vjust = 0.5))

# Helper function to create a centered text plot
title_plot <- function(label_text, text_size = 7, font_face = "bold") {
  ggplot() +
    annotate("text", x = 0.5, y = 0.5, 
             label = label_text, 
             size = text_size, fontface = font_face, hjust = 0.5, vjust = 0.5) +
    theme_void()
}

############Projection-focused time series plots###########

# Set axis limits for SSB/SSB_MSY plots (increasing and decreasing scenarios)
increase_coord_row1 <- coord_cartesian(xlim = c(1,51), ylim = c(0,3.5))
increase_coord_row2 <- coord_cartesian(xlim = c(1,51), ylim = c(0,3.5))
increase_coord_row3 <- coord_cartesian(xlim = c(1,51), ylim = c(0,5.7))
increase_coord_row4 <- coord_cartesian(xlim = c(1,51), ylim = c(0,5.7))
decrease_coord_row1 <- coord_cartesian(xlim = c(1,51), ylim = c(0,1.2))
decrease_coord_row2 <- coord_cartesian(xlim = c(1,51), ylim = c(0,1.2))
decrease_coord_row3 <- coord_cartesian(xlim = c(1,51), ylim = c(0,2.9))
decrease_coord_row4 <- coord_cartesian(xlim = c(1,51), ylim = c(0,2.9))

# Set axis limits for LRP/USR plots (increasing and decreasing scenarios)
increase_coord_lrp <- coord_cartesian(xlim = c(51,100), ylim = c(0,60)) 
increase_coord_lrp2 <- coord_cartesian(xlim = c(51,100), ylim = c(0,50)) 
decrease_coord_lrp <- coord_cartesian(xlim = c(51,100), ylim = c(0,17))
decrease_coord_lrp2 <- coord_cartesian(xlim = c(51,100), ylim = c(0,30))

# Set axis limits for other time series plots (increasing and decreasing scenarios)
coord_LRP <- coord_cartesian(xlim = c(46,NA), ylim = c(0,1))
increase_coord_SSB <- coord_cartesian(xlim = c(46,NA), ylim = c(0,6))
increase_coord_AVY <- coord_cartesian(xlim = c(46,NA), ylim = c(0,4.5))
increase_coord_catch <- coord_cartesian(xlim = c(46,NA), ylim = c(0,15))
decrease_coord_SSB <- coord_cartesian(xlim = c(46,NA), ylim = c(0,2.25))
decrease_coord_AVY <- coord_cartesian(xlim = c(46,NA), ylim = c(0,4))
decrease_coord_catch <- coord_cartesian(xlim = c(46,NA), ylim = c(0,4.5))

# Arrange increasing scenarios' P(SSB>LRP), SSB/SSBMSYinitial, P(catch = 0), catch,
# and AVY into a grid and save
(growthup_lowF_plots[["Prob_SSB_LRP"]] + guides(col = "none") +
    theme_no_x + coord_LRP +
    ggtitle("Increasing productivity\n(G & R), non-rebuilding")) +
  (growthup_highF_plots[["Prob_SSB_LRP"]] + 
    theme_no_xy + guides(color = "none") + coord_LRP + 
    ggtitle("Increasing productivity\n(G & R), rebuilding")) + 
  (decreaseM_lowF_plots[["Prob_SSB_LRP"]] + 
    theme_no_xy + guides(color = "none") + coord_LRP + 
    ggtitle("Increasing productivity\n(M), non-rebuilding")) +
  (decreaseM_highF_plots[["Prob_SSB_LRP"]] + 
    theme_no_xy + guides(color = "none") + coord_LRP + 
    ggtitle("Increasing productivity\n(M), rebuilding")) +
  (growthup_lowF_plots[["Median_SSBMSY"]] + theme_no_x + 
     guides(col = guide_legend(nrow = 1)) + increase_coord_SSB) + 
  (growthup_highF_plots[["Median_SSBMSY"]] + theme_no_xy + guides(color = "none") +
    increase_coord_SSB) +
  (decreaseM_lowF_plots[["Median_SSBMSY"]] + theme_no_xy + guides(color = "none") +
    increase_coord_SSB) +
  (decreaseM_highF_plots[["Median_SSBMSY"]] + theme_no_xy + guides(color = "none") +
    increase_coord_SSB) +
  (growthup_lowF_plots[["Prob_Catch"]] + theme_no_x + guides(color = "none") + 
     coord_LRP) + 
  (growthup_highF_plots[["Prob_Catch"]] + theme_no_xy + guides(color = "none") + 
     coord_LRP) + 
  (decreaseM_lowF_plots[["Prob_Catch"]] + theme_no_xy + guides(color = "none") + 
     coord_LRP) +
  (decreaseM_highF_plots[["Prob_Catch"]] + theme_no_xy + guides(color = "none") + 
     coord_LRP) +
  (growthup_lowF_plots[["Median_Catch"]] + theme_no_x + guides(color = "none") + 
     increase_coord_catch) + 
  (growthup_highF_plots[["Median_Catch"]] + theme_no_xy + guides(color = "none") + 
     increase_coord_catch) + 
  (decreaseM_lowF_plots[["Median_Catch"]] + theme_no_xy + guides(color = "none") + 
     increase_coord_catch) +
  (decreaseM_highF_plots[["Median_Catch"]] + theme_no_xy + guides(color = "none") + 
     increase_coord_catch) +
  (growthup_lowF_plots[["AVY"]] + guides(color = "none") + 
     increase_coord_AVY + theme(axis.text.x = element_text(angle = 90, vjust = 0.5))) + 
  (growthup_highF_plots[["AVY"]] + theme_no_y + guides(color = "none") + 
     increase_coord_AVY) + 
  (decreaseM_lowF_plots[["AVY"]] + theme_no_y + guides(color = "none") + 
     increase_coord_AVY) +
  (decreaseM_highF_plots[["AVY"]] + theme_no_y + guides(color = "none") + 
     increase_coord_AVY) +
  plot_layout(guides = "collect", ncol = 4, nrow = 5) &
  theme(plot.title = element_text(size = font_size, face = "bold"),
        text = element_text(size = font_size),
        legend.text = element_text(size = font_size),
        legend.position = "bottom", legend.key.size = unit(1.5, "cm"))
ggsave(paste0(getwd(),"/Figs/MSE/combined_MSE_1A.jpeg"), plot=last_plot(), 
       width=fig_width, height=fig_height, units="in")

# Arrange decreasing scenarios' P(SSB>LRP), SSB/SSBMSYinitial, P(catch = 0), catch,
# and AVY into a grid and save
(growthdown_lowF_plots[["Prob_SSB_LRP"]] + guides(color = "none") +
    theme_no_x + coord_LRP +
    ggtitle("Decreasing productivity\n(G & R), non-rebuilding")) +
  (growthdown_highF_plots[["Prob_SSB_LRP"]] + 
     theme_no_xy + guides(color = "none") + coord_LRP + 
     ggtitle("Decreasing productivity\n(G & R), rebuilding")) + 
  (increaseM_lowF_plots[["Prob_SSB_LRP"]] + 
     theme_no_xy + guides(color = "none") + coord_LRP + 
     ggtitle("Decreasing productivity\n(M), non-rebuilding")) +
  (increaseM_highF_plots[["Prob_SSB_LRP"]] + 
     theme_no_xy + guides(color = "none") + coord_LRP + 
     ggtitle("Decreasing productivity\n(M), rebuilding")) +
  (growthdown_lowF_plots[["Median_SSBMSY"]] + theme_no_x + 
     guides(col = guide_legend(nrow = 1)) + decrease_coord_SSB) + 
  (growthdown_highF_plots[["Median_SSBMSY"]] + theme_no_xy + guides(color = "none") +
     decrease_coord_SSB) +
  (increaseM_lowF_plots[["Median_SSBMSY"]] + theme_no_xy + guides(color = "none") +
     decrease_coord_SSB) +
  (increaseM_highF_plots[["Median_SSBMSY"]] + theme_no_xy + guides(color = "none") +
     decrease_coord_SSB) +
  (growthdown_lowF_plots[["Prob_Catch"]] + theme_no_x + guides(color = "none") + 
     coord_LRP) + 
  (growthdown_highF_plots[["Prob_Catch"]] + theme_no_xy + guides(color = "none") + 
     coord_LRP) + 
  (increaseM_lowF_plots[["Prob_Catch"]] + theme_no_xy + guides(color = "none") + 
     coord_LRP) +
  (increaseM_highF_plots[["Prob_Catch"]] + theme_no_xy + guides(color = "none") + 
     coord_LRP) +
  (growthdown_lowF_plots[["Median_Catch"]] + theme_no_x + guides(color = "none") + 
     decrease_coord_catch) + 
  (growthdown_highF_plots[["Median_Catch"]] + theme_no_xy + guides(color = "none") + 
     decrease_coord_catch) + 
  (increaseM_lowF_plots[["Median_Catch"]] + theme_no_xy + guides(color = "none") + 
     decrease_coord_catch) +
  (increaseM_highF_plots[["Median_Catch"]] + theme_no_xy + guides(color = "none") + 
     decrease_coord_catch) +
  (growthdown_lowF_plots[["AVY"]] + guides(color = "none") + 
     decrease_coord_AVY + theme(axis.text.x = element_text(angle = 90, vjust = 0.5))) + 
  (growthdown_highF_plots[["AVY"]] + theme_no_y + guides(color = "none") + 
     decrease_coord_AVY) + 
  (increaseM_lowF_plots[["AVY"]] + theme_no_y + guides(color = "none") + 
     decrease_coord_AVY) +
  (increaseM_highF_plots[["AVY"]] + theme_no_y + guides(color = "none") + 
     decrease_coord_AVY) +
  plot_layout(guides = "collect", ncol = 4, nrow = 5) &
  theme(plot.title = element_text(size = font_size, face = "bold"),
        text = element_text(size = font_size),
        legend.text = element_text(size = font_size),
        legend.position = "bottom", legend.key.size = unit(1.5, "cm"))
ggsave(paste0(getwd(),"/Figs/MSE/combined_MSE_2A.jpeg"), plot=last_plot(), 
       width=fig_width, height=fig_height, units="in")

# Arrange increasing scenarios' SSB/SSB_MSY plots into a grid and save
(title_plot("Increasing\nproductivity\n(G & R),\nnon-rebuilding") + 
    growthup_lowF_plots[["SSB_SSBMSY_First5"]] + guides(col = guide_legend(nrow = 1)) +
    theme_no_x + theme(axis.title.y = element_blank()) + increase_coord_row1 +
    ggtitle(bquote(SSB / SSB[.("MSYinitial")]))) +
  (growthup_lowF_plots[["SSB_SSBMSY_EntireTS"]] + 
     theme_no_xy + guides(color = "none") + increase_coord_row1 + 
     ggtitle(bquote(SSB / SSB[.("MSYfull")]))) + 
  (growthup_lowF_plots[["SSB_SSBMSY_Last5"]] + 
     theme_no_xy + guides(color = "none") + increase_coord_row1 + 
     ggtitle(bquote(SSB / SSB[.("MSYrecent")]))) +
  (title_plot("Increasing\nproductivity\n(G & R),\nrebuilding")) +
     (growthup_highF_plots[["SSB_SSBMSY_First5"]] + 
     theme_no_x + guides(color = "none") + 
       theme(axis.title.y = element_blank()) + increase_coord_row2) +
  (growthup_highF_plots[["SSB_SSBMSY_EntireTS"]] + 
     theme_no_xy + guides(color = "none") + increase_coord_row2) + 
  (growthup_highF_plots[["SSB_SSBMSY_Last5"]] + 
     theme_no_xy + guides(color = "none") + increase_coord_row2) +
  (title_plot("Increasing\nproductivity (M),\nnon-rebuilding")) + 
     (decreaseM_lowF_plots[["SSB_SSBMSY_First5"]] + 
       theme_no_x + guides(color = "none") + 
       theme(axis.title.y = element_blank()) + increase_coord_row3) +
  (decreaseM_lowF_plots[["SSB_SSBMSY_EntireTS"]] + 
     theme_no_xy + guides(color = "none") + increase_coord_row3) + 
  (decreaseM_lowF_plots[["SSB_SSBMSY_Last5"]] + 
     theme_no_xy + guides(color = "none") + increase_coord_row3) +
  (title_plot("Increasing\nproductivity (M),\nrebuilding")) + 
     (decreaseM_highF_plots[["SSB_SSBMSY_First5"]] + 
       guides(color = "none") + theme(axis.title.y = element_blank(), 
                                      axis.text.x = element_text(angle = 90, 
                                                                 vjust = 0.5)) + 
       increase_coord_row4) +
  (decreaseM_highF_plots[["SSB_SSBMSY_EntireTS"]] + 
     theme_no_y + guides(color = "none") + increase_coord_row4) + 
  (decreaseM_highF_plots[["SSB_SSBMSY_Last5"]] + 
     theme_no_y + guides(color = "none") + increase_coord_row4) +
  plot_layout(guides = "collect", ncol = 4, nrow = 4,
              widths = c(0.85,1,1,1)) &
  theme(plot.title = element_text(size = font_size, face = "bold"),
        text = element_text(size = font_size),
        legend.text = element_text(size = font_size),
        legend.position = "bottom", legend.key.size = unit(1.5, "cm"))
ggsave(paste0(getwd(),"/Figs/MSE/combined_MSE_1C.jpeg"), plot=last_plot(), 
       width=fig_width_c, height=fig_height_c, units="in")

# Arrange decreasing scenarios' SSB/SSB_MSY plots into a grid and save
(title_plot("Decreasing\nproductivity\n(G & R),\nnon-rebuilding") + 
    growthdown_lowF_plots[["SSB_SSBMSY_First5"]] + guides(col = guide_legend(nrow = 1)) +
    theme_no_x + theme(axis.title.y = element_blank()) + decrease_coord_row1 +
    ggtitle(bquote(SSB / SSB[.("MSYinitial")]))) +
  (growthdown_lowF_plots[["SSB_SSBMSY_EntireTS"]] + 
     theme_no_xy + guides(color = "none") + decrease_coord_row1 + 
     ggtitle(bquote(SSB / SSB[.("MSYfull")]))) + 
  (growthdown_lowF_plots[["SSB_SSBMSY_Last5"]] + 
     theme_no_xy + guides(color = "none") + decrease_coord_row1 + 
     ggtitle(bquote(SSB / SSB[.("MSYrecent")]))) +
  (title_plot("Decreasing\nproductivity\n(G & R),\nrebuilding")) +
  (growthdown_highF_plots[["SSB_SSBMSY_First5"]] + 
     theme_no_x + guides(color = "none") + 
     theme(axis.title.y = element_blank()) + decrease_coord_row2) +
  (growthdown_highF_plots[["SSB_SSBMSY_EntireTS"]] + 
     theme_no_xy + guides(color = "none") + decrease_coord_row2) + 
  (growthdown_highF_plots[["SSB_SSBMSY_Last5"]] + 
     theme_no_xy + guides(color = "none") + decrease_coord_row2) +
  (title_plot("Decreasing\nproductivity (M),\nnon-rebuilding")) + 
  (increaseM_lowF_plots[["SSB_SSBMSY_First5"]] + 
     theme_no_x + guides(color = "none") + 
     theme(axis.title.y = element_blank()) + decrease_coord_row3) +
  (increaseM_lowF_plots[["SSB_SSBMSY_EntireTS"]] + 
     theme_no_xy + guides(color = "none") + decrease_coord_row3) + 
  (increaseM_lowF_plots[["SSB_SSBMSY_Last5"]] + 
     theme_no_xy + guides(color = "none") + decrease_coord_row3) +
  (title_plot("Decreasing\nproductivity (M),\nrebuilding")) + 
  (increaseM_highF_plots[["SSB_SSBMSY_First5"]] + 
     guides(color = "none") + theme(axis.title.y = element_blank(), 
                                    axis.text.x = element_text(angle = 90, 
                                                               vjust = 0.5)) + 
     decrease_coord_row4) +
  (increaseM_highF_plots[["SSB_SSBMSY_EntireTS"]] + 
     theme_no_y + guides(color = "none") + decrease_coord_row4) + 
  (increaseM_highF_plots[["SSB_SSBMSY_Last5"]] + 
     theme_no_y + guides(color = "none") + decrease_coord_row4) +
  plot_layout(guides = "collect", ncol = 4, nrow = 4,
              widths = c(0.85,1,1,1)) &
  theme(plot.title = element_text(size = font_size, face = "bold"),
        text = element_text(size = font_size),
        legend.text = element_text(size = font_size),
        legend.position = "bottom", legend.key.size = unit(1.5, "cm"))
ggsave(paste0(getwd(),"/Figs/MSE/combined_MSE_2C.jpeg"), plot=last_plot(), 
       width=fig_width_c, height=fig_height_c, units="in")

# Arrange increasing and decreasing scenario LRP/USR plots into a single grid and save
(growthup_lowF_plots[["Median_SSB_LRPs"]] + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
    ggtitle("Increasing productivity\n(G & R), non-rebuilding") + increase_coord_lrp) + 
  (growthup_highF_plots[["Median_SSB_LRPs"]] + theme_no_y + 
     theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
     ggtitle("Increasing productivity\n(G & R), rebuilding") + increase_coord_lrp) + 
  (decreaseM_lowF_plots[["Median_SSB_LRPs"]] + 
     theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
     ggtitle("Increasing productivity\n(M), non-rebuilding") + increase_coord_lrp2) +
  (decreaseM_highF_plots[["Median_SSB_LRPs"]] + theme_no_y + 
     theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
     ggtitle("Increasing productivity\n(M), rebuilding") + increase_coord_lrp2) +
  (growthdown_lowF_plots[["Median_SSB_LRPs"]] + 
     theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
     ggtitle("Decreasing productivity\n(G & R), non-rebuilding") + decrease_coord_lrp) + 
  (growthdown_highF_plots[["Median_SSB_LRPs"]] + theme_no_y +
     ggtitle("Decreasing productivity\n(G & R), rebuilding") + decrease_coord_lrp) + 
  (increaseM_lowF_plots[["Median_SSB_LRPs"]] + 
     theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
     ggtitle("Decreasing productivity\n(M), non-rebuilding") + decrease_coord_lrp2) +
  (increaseM_highF_plots[["Median_SSB_LRPs"]] + theme_no_y +
     ggtitle("Decreasing productivity\n(M), rebuilding") + decrease_coord_lrp2) +
  plot_layout(nrow = 4, guides = "collect") &
  theme(plot.title = element_text(size = font_size, face = "bold"),
        text = element_text(size = font_size),
        legend.text = element_text(size = font_size),
        legend.position = "right", legend.key.size = unit(1.5, "cm"))
ggsave("combined_MSE_LRPs_USRs.jpeg", plot=last_plot(), 
       width=fig_width_c*0.8, height=fig_height_c*1.2, units="in")

############Full time series plots###########

# Set axis limits for plots (increasing and decreasing scenarios)
increase_coord_SSB_zoomout <- coord_cartesian(ylim = c(0,65))
increase_coord_F_zoomout <- coord_cartesian(ylim = c(0,0.7))
increase_coord_catch_zoomout <- coord_cartesian(ylim = c(0,13))
decrease_coord_SSB_zoomout <- coord_cartesian(ylim = c(0,90))
decrease_coord_F_zoomout <- coord_cartesian(ylim = c(0,0.55))
decrease_coord_catch_zoomout <- coord_cartesian(ylim = c(0,11))

# Arrange increasing scenarios' SSB, F, and catch (historical and projected) 
# into a grid and save
(growthup_lowF_plots[["Median_SSB"]] + guides(col = guide_legend(nrow = 1)) +
    theme_no_x + increase_coord_SSB_zoomout +
    ggtitle("Increasing productivity\n(G & R), non-rebuilding")) +
  (growthup_highF_plots[["Median_SSB"]] + 
     theme_no_xy + guides(color = "none") + increase_coord_SSB_zoomout + 
     ggtitle("Increasing productivity\n(G & R), rebuilding")) + 
  (decreaseM_lowF_plots[["Median_SSB"]] + 
     theme_no_xy + guides(color = "none") + increase_coord_SSB_zoomout + 
     ggtitle("Increasing productivity\n(M), non-rebuilding")) +
  (decreaseM_highF_plots[["Median_SSB"]] + 
     theme_no_xy + guides(color = "none") + increase_coord_SSB_zoomout + 
     ggtitle("Increasing productivity\n(M), rebuilding")) +
  (growthup_lowF_plots[["Median_F"]] + theme_no_x + guides(color = "none") + 
     increase_coord_F_zoomout) + 
  (growthup_highF_plots[["Median_F"]] + theme_no_xy + guides(color = "none") +
     increase_coord_F_zoomout) +
  (decreaseM_lowF_plots[["Median_F"]] + theme_no_xy + guides(color = "none") +
     increase_coord_F_zoomout) +
  (decreaseM_highF_plots[["Median_F"]] + theme_no_xy + guides(color = "none") +
     increase_coord_F_zoomout) +
  (growthup_lowF_plots[["Median_Catch"]] + 
     theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) + 
     guides(color = "none") + increase_coord_catch_zoomout) + 
  (growthup_highF_plots[["Median_Catch"]] + theme_no_y + guides(color = "none") + 
     increase_coord_catch_zoomout) + 
  (decreaseM_lowF_plots[["Median_Catch"]] + theme_no_y + guides(color = "none") + 
     increase_coord_catch_zoomout) +
  (decreaseM_highF_plots[["Median_Catch"]] + theme_no_y + guides(color = "none") + 
     increase_coord_catch_zoomout) +
  plot_layout(guides = "collect", ncol = 4, nrow = 3) &
  theme(plot.title = element_text(size = font_size, face = "bold"),
        text = element_text(size = font_size),
        legend.text = element_text(size = font_size),
        legend.position = "bottom", legend.key.size = unit(1.5, "cm"))
ggsave(paste0(getwd(),"/Figs/MSE/combined_MSE_1B.jpeg"), plot=last_plot(), 
       width=fig_width_f, height=fig_height_f, units="in")

# Arrange decreasing scenarios' SSB, F, and catch (historical and projected) 
# into a grid and save
(growthdown_lowF_plots[["Median_SSB"]] + guides(col = guide_legend(nrow = 1)) +
    theme_no_x + decrease_coord_SSB_zoomout +
    ggtitle("Decreasing productivity\n(G & R), non-rebuilding")) +
  (growthdown_highF_plots[["Median_SSB"]] + 
     theme_no_xy + guides(color = "none") + decrease_coord_SSB_zoomout + 
     ggtitle("Decreasing productivity\n(G & R), rebuilding")) + 
  (increaseM_lowF_plots[["Median_SSB"]] + 
     theme_no_xy + guides(color = "none") + decrease_coord_SSB_zoomout + 
     ggtitle("Decreasing productivity\n(M), non-rebuilding")) +
  (increaseM_highF_plots[["Median_SSB"]] + 
     theme_no_xy + guides(color = "none") + decrease_coord_SSB_zoomout + 
     ggtitle("Decreasing productivity\n(M), rebuilding")) +
  (growthdown_lowF_plots[["Median_F"]] + theme_no_x + guides(color = "none") + 
     decrease_coord_F_zoomout) + 
  (growthdown_highF_plots[["Median_F"]] + theme_no_xy + guides(color = "none") +
     decrease_coord_F_zoomout) +
  (increaseM_lowF_plots[["Median_F"]] + theme_no_xy + guides(color = "none") +
     decrease_coord_F_zoomout) +
  (increaseM_highF_plots[["Median_F"]] + theme_no_xy + guides(color = "none") +
     decrease_coord_F_zoomout) +
  (growthdown_lowF_plots[["Median_Catch"]] + 
     theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) + 
     guides(color = "none") + decrease_coord_catch_zoomout) + 
  (growthdown_highF_plots[["Median_Catch"]] + theme_no_y + guides(color = "none") + 
     decrease_coord_catch_zoomout) + 
  (increaseM_lowF_plots[["Median_Catch"]] + theme_no_y + guides(color = "none") + 
     decrease_coord_catch_zoomout) +
  (increaseM_highF_plots[["Median_Catch"]] + theme_no_y + guides(color = "none") + 
     decrease_coord_catch_zoomout) +
  plot_layout(guides = "collect", ncol = 4, nrow = 3) &
  theme(plot.title = element_text(size = font_size, face = "bold"),
        text = element_text(size = font_size),
        legend.text = element_text(size = font_size),
        legend.position = "bottom", legend.key.size = unit(1.5, "cm"))
ggsave(paste0(getwd(),"/Figs/MSE/combined_MSE_2B.jpeg"), plot=last_plot(), 
       width=fig_width_f, height=fig_height_f, units="in")

############Radar plots###########

# Arrange increasing scenarios' radar plots into a grid and save
plot_spacer() + 
  title_plot("Increasing productivity\n(G & R), non-rebuilding") +
  title_plot("Increasing productivity\n(G & R), rebuilding") +
  title_plot("Increasing productivity\n(M), non-rebuilding") +
  title_plot("Increasing productivity\n(M), rebuilding") +
  title_plot("ST") + 
  growthup_lowF_plots[["Radar_Plot_ST"]] +
  growthup_highF_plots[["Radar_Plot_ST"]] +
  decreaseM_lowF_plots[["Radar_Plot_ST"]] +
  decreaseM_highF_plots[["Radar_Plot_ST"]] +
  plot_spacer() + 
  plot_spacer() + 
  plot_spacer() + 
  plot_spacer() + 
  plot_spacer() +
  title_plot("LT") +
  growthup_lowF_plots[["Radar_Plot_LT"]] +
  growthup_highF_plots[["Radar_Plot_LT"]] +
  decreaseM_lowF_plots[["Radar_Plot_LT"]] +
  decreaseM_highF_plots[["Radar_Plot_LT"]] +
  plot_layout(nrow = 4, ncol = 5, guides = "collect", 
              heights = c(0.3,1,0.012,1),
              widths = c(0.2,1,1,1,1)) &
  theme(plot.title = element_text(size = font_size, face = "bold"),
        text = element_text(size = font_size),
        legend.text = element_text(size = font_size),
        legend.position = "bottom", legend.key.size = unit(1.5, "cm"))
ggsave(paste0(getwd(),"/Figs/MSE/combined_MSE_radar1.jpeg"), plot=last_plot(), 
       width=fig_width_r, height=fig_height_r, units="in")

# Arrange decreasing scenarios' radar plots into a grid and save
plot_spacer() + 
  title_plot("Decreasing productivity\n(G & R), non-rebuilding") +
  title_plot("Decreasing productivity\n(G & R), rebuilding") +
  title_plot("Decreasing productivity\n(M), non-rebuilding") +
  title_plot("Decreasing productivity\n(M), rebuilding") +
  title_plot("ST") + 
  growthdown_lowF_plots[["Radar_Plot_ST"]] +
  growthdown_highF_plots[["Radar_Plot_ST"]] +
  increaseM_lowF_plots[["Radar_Plot_ST"]] +
  increaseM_highF_plots[["Radar_Plot_ST"]] +
  plot_spacer() + 
  plot_spacer() + 
  plot_spacer() + 
  plot_spacer() + 
  plot_spacer() +
  title_plot("LT") +
  growthdown_lowF_plots[["Radar_Plot_LT"]] +
  growthdown_highF_plots[["Radar_Plot_LT"]] +
  increaseM_lowF_plots[["Radar_Plot_LT"]] +
  increaseM_highF_plots[["Radar_Plot_LT"]] +
  plot_layout(nrow = 4, ncol = 5, guides = "collect", 
              heights = c(0.3,1,0.012,1),
              widths = c(0.2,1,1,1,1)) &
  theme(plot.title = element_text(size = font_size, face = "bold"),
        text = element_text(size = font_size),
        legend.text = element_text(size = font_size),
        legend.position = "bottom", legend.key.size = unit(1.5, "cm"))
ggsave(paste0(getwd(),"/Figs/MSE/combined_MSE_radar2.jpeg"), plot=last_plot(), 
       width=fig_width_r, height=fig_height_r, units="in")