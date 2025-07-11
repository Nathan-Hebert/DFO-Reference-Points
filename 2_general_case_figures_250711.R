#
# General case figures using 3 simulated stocks
#
################################################################################

# Load in necessary libraries
library(ggplot2); library(viridis); library(patchwork)
library(colorspace)

# Set up the 3 simulated stocks, and compute their RPs generated from a 
# grid of waa, mat, M, rec params, and sel... store their RPs in separate dataframes
stock_list <- c("S","M","L")
for (i in 1:length(stock_list))
{
  stock <- stock_list[i]
  source('1_general_case_setup_250127.R')
}

# Load in and combine the dataframes, and then calculate h
DF <- rbind(readRDS("simulated stock RPs/DF_stockS.rds"), 
            readRDS("simulated stock RPs/DF_stockM.rds"),
            readRDS("simulated stock RPs/DF_stockL.rds"))
DF$h <- (DF$alpha*DF$phi_0)/(DF$alpha*DF$phi_0 + 4)

# Rename factor levels for stock
levels(DF$stock) <- c("Short-lived", "Medium-lived", "Long-lived")

# Grab base case rows
BASE <- DF[which(DF$M_factor == 1.0 & DF$waa_factor == 1.0 & DF$sel_shift == 0
                 & DF$mat_shift == 0 & DF$recdev_factor == 1.0),]

# Values required for figures
font_size <- 16; point_size <- 2; line_width <- 1.2
fig_width <- 10; fig_height <- 6 # For everything except the SSB0 waa plot
theme_no_x <- theme(axis.title.x = element_blank(), axis.ticks.x = element_blank(),
                    axis.text.x = element_blank())
theme_no_y <- theme(axis.title.y = element_blank(), axis.ticks.y = element_blank(),
                    axis.text.y = element_blank())
theme_no_x_label <- theme(axis.title.x = element_blank())
RP_names <- c("phi_0", "SSB_0","F_MSY", "SSB_MSY", 
              "F_X.SPR", "R0", "MSY", "h")
ylabels <- c(expression("Relative \u03A6"["0"]), expression("Relative SSB"["0"]),
             expression("Relative F"["MSY"]), expression("Relative SSB"["MSY"]), 
             bquote("Relative F"[.(40)*"%SPR"]), expression("Relative R"["0"]),"Relative MSY",
             "Relative h")
color_palette <- darken(viridis(length(unique(DF$stock))), amount = 0.25) # For everything except 
                                                                          # the SSB0 waa plot

# Function to generate a plot illustrating the effect on an RP of changing
# one component of productivity while holding the others fixed... all 3
# simulated stocks are shown. DF_RP is the dataframe containing the data including the 
# RPs, BASE_DF is the subset of DF_RP containing the base case rows, col_name is 
# the name of the column containing the RP of interest, and vary is what productivity 
# component to change (i.e., either "M", "waa", "mat", "rec", or "sel")
relative_effects_plot <- function(DF_RP = DF, BASE_DF = BASE, 
                                  col_name, vary = "M", ylim, xlab, ylab,
                                  wid_lines = line_width, size_points = point_size,
                                  text_size = font_size)
{
  # Determine what to hold constant and what to vary
  switch(vary,
         "M" = {
           hold1 <- "waa_factor"
           hold2 <- "sel_shift"
           hold3 <- "mat_shift"
           hold4 <- "recdev_factor"
           DF_RP$sel_shift <- DF_RP$sel_shift + 1 # So sel_shift = 1 means no shift
           DF_RP$mat_shift <- DF_RP$mat_shift + 1 # So mat_shift = 1 means no shift
         },
         "waa" = {
           vary <- "waa_factor"
           hold1 <- "M_factor"
           hold2 <- "sel_shift"
           hold3 <- "mat_shift"
           hold4 <- "recdev_factor"
           DF_RP$sel_shift <- DF_RP$sel_shift + 1 # So sel_shift = 1 means no shift
           DF_RP$mat_shift <- DF_RP$mat_shift + 1 # So mat_shift = 1 means no shift
         },
         "sel" = {
           vary <- "sel_shift"
           hold1 <- "M_factor"
           hold2 <- "waa_factor"
           hold3 <- "mat_shift"
           hold4 <- "recdev_factor"
           DF_RP$mat_shift <- DF_RP$mat_shift + 1 # So mat_shift = 1 means no shift
         },
         "mat" = {
           vary <- "mat_shift"
           hold1 <- "M_factor"
           hold2 <- "waa_factor"
           hold3 <- "sel_shift"
           hold4 <- "recdev_factor"
           DF_RP$sel_shift <- DF_RP$sel_shift + 1 # So sel_shift = 1 means no shift
         },
         "rec" = {
           vary <- "recdev_factor"
           hold1 <- "M_factor"
           hold2 <- "waa_factor"
           hold3 <- "sel_shift"
           hold4 <- "mat_shift"
           DF_RP$sel_shift <- DF_RP$sel_shift + 1 # So sel_shift = 1 means no shift
           DF_RP$mat_shift <- DF_RP$mat_shift + 1 # So mat_shift = 1 means no shift
         },
         stop("'vary' must be either 'waa', 'sel', 'mat', 'rec', or 'M'.")
  )
  
  # Calculate relative values
  for (j in 1:nrow(BASE_DF))
  {
    dfx <- which(DF_RP$stock == BASE_DF$stock[j] & DF_RP[[hold1]] == 1.0 & 
                   DF_RP[[hold2]] == 1.0 & DF_RP[[hold3]] == 1.0 & DF_RP[[hold4]] == 1.0)
    DF_RP[dfx, "relative_val"] <- (DF_RP[dfx, col_name])/BASE_DF[j, col_name]
  }
  
  # Generate plot
  ggplot(DF_RP[which(DF_RP[[hold1]] == 1.0 & DF_RP[[hold2]] == 1.0 &
                       DF_RP[[hold3]] == 1.0 & DF_RP[[hold4]] == 1.0),], 
              aes(x = .data[[vary]],y = relative_val, 
                  col = as.factor(stock))) +
    geom_point(size = size_points) + geom_line(size = wid_lines) + theme_classic() + 
    theme(text = element_text(size = text_size)) + geom_hline(yintercept = 1) +
    scale_color_manual(values = color_palette, name = "Life-history") +
    labs(x = xlab, y = ylab) + scale_y_continuous(limits = ylim, expand = c(0,0))
}

##################Effect of Changing WAA on Reference Points####################

# Create individual plots of waa factor vs RPs with sel, rec, mat, and M fixed at base... 
# values are relative to base
plots <- list()
for (i in 1:length(ylabels)) {
  plots[[RP_names[i]]] <- relative_effects_plot(col_name = RP_names[i], 
                                                vary = "waa", ylim = c(0.2, 1.9), 
                                                xlab = "Weight-at-age multiplier c", 
                                                ylab = ylabels[i])
}
# Combine plots and save
(plots[["phi_0"]] + theme_no_x) + 
  (plots[["F_MSY"]] + theme_no_x) + 
  (plots[["F_X.SPR"]] + theme_no_x) +
  (plots[["SSB_0"]] + theme_no_x) + 
  (plots[["SSB_MSY"]] + theme_no_x) + 
  guide_area() + 
  (plots[["R0"]] + theme_no_x_label) +
  (plots[["MSY"]]) + 
  (plots[["h"]] + theme_no_x_label) +
  plot_layout(ncol = 3, nrow = 3, guides = "collect") & 
  theme(legend.background = element_rect(fill = "white", color = "black"))
ggsave(paste0(getwd(),"/figs/general case/waa_general_combined.jpeg"), 
       plot=last_plot(), width=fig_width, height=fig_height, units="in")

# Plot relationship between steepness h and new SSB_0/old SSB_0 when multiply phi_0
# by a constant c (such as above)... then save figure
c_values <- c(1.5,1.25,1.00,0.75,0.5)
color_palette2 <- darken(viridis(length(c_values)), amount = 0.25)
df <- data.frame(x = rep(seq(0.18, 1, length.out = 100), length(c_values)))
df$c <- rep(c_values, each = 100)
curve_function <- function(x, c) {
  ((4 * c * x - 1 + x) / (5 * x - 1))
}
ggplot(df, aes(x = x, y = curve_function(x, c), color = rev(as.factor(c)))) +
  theme_classic() +
  labs(x = "Original h", y = expression("New SSB"[0] / "Original SSB"[0])) +
  scale_x_continuous(breaks = seq(0.25, 1, 0.25), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(-3,5, 1), expand = c(0, 0)) +
  coord_cartesian(xlim = c(0.21, 1), ylim = c(0,4)) +
  geom_line(size = line_width) +
  scale_color_manual(values = rev(color_palette2), name = "c", labels = sprintf("%.2f", 
                                                                               c_values)) +  
  theme(text = element_text(size = font_size),
        panel.grid.major.y = element_line(color = "slategray", linetype = "dashed"),
        panel.grid.minor.y = element_line(color = "slategray", linetype = "dashed")) +
  theme(text = element_text(size = font_size), legend.position = c(0.87, 0.75),
        legend.background = element_rect(fill = "white", color = "black"))
ggsave(paste0(getwd(),"/figs/general case/waa_general_SSB0.jpeg"), 
       plot=last_plot(), width=8, height=6, units="in")

##################Effect of Shifting Maturity on Reference Points###############

# Create individual plots of mat shift vs RPs with M, waa, rec, and sel fixed at base... 
# values are relative to base
plots <- list()
for (i in 1:length(ylabels)) {
  plots[[RP_names[i]]] <- relative_effects_plot(col_name = RP_names[i], 
                                                vary = "mat", xlab = "\u394 in age", 
                                                ylab = ylabels[i],
                                                ylim = c(0.70, 1.5))
}
# Combine plots and save
(plots[["phi_0"]] + theme_no_x) + 
  (plots[["F_MSY"]] + theme_no_x) + 
  (plots[["F_X.SPR"]] + theme_no_x) +
  (plots[["SSB_0"]] + theme_no_x) + 
  (plots[["SSB_MSY"]] + theme_no_x) + 
  guide_area() + 
  (plots[["R0"]] + theme_no_x_label) +
  (plots[["MSY"]]) + 
  (plots[["h"]] + theme_no_x_label) +
  plot_layout(ncol = 3, nrow = 3, guides = "collect") & 
  theme(legend.background = element_rect(fill = "white", color = "black"))
ggsave(paste0(getwd(),"/figs/general case/mat_general_combined.jpeg"), 
       plot=last_plot(), width=fig_width, height=fig_height, units="in") 

###############Effect of Shifting Selectivity on Reference Points###############

# Create individual plots of sel shift vs RPs with M, rec, waa, and mat fixed at base... 
# values are relative to base
plots <- list()
for (i in 1:length(ylabels)) {
  plots[[RP_names[i]]] <- relative_effects_plot(col_name = RP_names[i], 
                                                vary = "sel", xlab = "\u394 in age", 
                                                ylab = ylabels[i],
                                                ylim = c(0.4, 1.8))
}
# Combine plots and save
(plots[["phi_0"]] + theme_no_x) + 
  (plots[["F_MSY"]] + theme_no_x) + 
  (plots[["F_X.SPR"]] + theme_no_x) +
  (plots[["SSB_0"]] + theme_no_x) + 
  (plots[["SSB_MSY"]] + theme_no_x) + 
  guide_area() + 
  (plots[["R0"]] + theme_no_x_label) +
  (plots[["MSY"]]) + 
  (plots[["h"]] + theme_no_x_label) +
  plot_layout(ncol = 3, nrow = 3, guides = "collect") & 
  theme(legend.background = element_rect(fill = "white", color = "black"))
ggsave(paste0(getwd(),"/figs/general case/sel_general_combined.jpeg"), 
       plot=last_plot(), width=fig_width, height=fig_height, units="in") 

##################Effect of Changing M on Reference Points######################

# Create individual plots of M vs RPs with mat, rec, sel, and waa fixed at base... 
# values are relative to base
plots <- list()
ylim <- list(ymin = c(0.01, 0.01, 0.52, 0.01, 0.52, 0.52, 0.01, 0.52),
             ymax = c(4, 4, 1.7, 4, 1.7, 1.7, 4, 1.7))
for (i in 1:length(ylabels)) {
  plots[[RP_names[i]]] <- relative_effects_plot(col_name = RP_names[i], 
                                                vary = "M", xlab = "M", 
                                                ylab = ylabels[i],
                                                ylim = c(ylim[["ymin"]][i], 
                                                         ylim[["ymax"]][i]))
}
# Combine plots and save
(plots[["phi_0"]] + theme_no_x) + 
  (plots[["F_MSY"]] + theme_no_x) + 
  (plots[["F_X.SPR"]] + theme_no_x) +
  (plots[["SSB_0"]] + theme_no_x) + 
  (plots[["SSB_MSY"]] + theme_no_x) + 
  guide_area() + 
  (plots[["R0"]] + theme_no_x_label) +
  (plots[["MSY"]]) + 
  (plots[["h"]] + theme_no_x_label) +
  plot_layout(ncol = 3, nrow = 3, guides = "collect") & 
  theme(legend.background = element_rect(fill = "white", color = "black"))
ggsave(paste0(getwd(),"/figs/general case/M_general_combined.jpeg"), 
       plot=last_plot(), width=fig_width, height=fig_height, units="in") 

#############Effect of Changing Recruitment on Reference Points#################

# Create individual plots of rec dev multiplier vs RPs with mat, M, sel, and waa fixed 
# at base... values are relative to base
plots <- list()
for (i in 1:length(ylabels)) {
  plots[[RP_names[i]]] <- relative_effects_plot(col_name = RP_names[i], 
                                                vary = "rec", 
                                                xlab = "Recruitment deviation multiplier c", 
                                                ylab = ylabels[i],
                                                ylim = c(0.12,1.92))
}
# Combine plots and save
(plots[["phi_0"]] + theme_no_x) + 
  (plots[["F_MSY"]] + theme_no_x) + 
  (plots[["F_X.SPR"]] + theme_no_x) +
  (plots[["SSB_0"]] + theme_no_x) + 
  (plots[["SSB_MSY"]] + theme_no_x) + 
  guide_area() + 
  (plots[["R0"]] + theme_no_x_label) +
  (plots[["MSY"]]) + 
  (plots[["h"]] + theme_no_x_label) +
  plot_layout(ncol = 3, nrow = 3, guides = "collect") & 
  theme(legend.background = element_rect(fill = "white", color = "black"))
ggsave(paste0(getwd(),"/figs/general case/rec_general_combined.jpeg"), 
       plot=last_plot(), width=fig_width, height=fig_height, units="in") 
