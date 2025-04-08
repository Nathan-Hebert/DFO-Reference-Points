#
# Generate figures illustrating an individual OM's MSE results
#
################################################################################

# Load necessary libraries
library(dplyr); library(ggplot2); library(DT); library(tidyr); library(forcats);
library(ggnewscale)

# Important parameters for figures
font_size <- 16
width <- 8 # Figure width in inches
height <- 6 # Figure height in inches
MP_labs <- c(expression(HCR[full]), expression(HCR[initial]), 
             expression(HCR[recent]), expression(HCR[recent-target]), 
             expression(HCR[recent-BPA]), expression(HCR["F=0"])) # Plot labels for MPs
maxyrp <- 50 # Number of years to plot projections for (up to OM@proyears)
LQ <- 0.05; UQ <- 0.95 # Quantiles for plotting alongside median in time series plots
q_alpha <- 0.2 # Alpha used for the quantiles in time series plots
lwd <- 1.1 # Line width used by plots
pt_size <- 2.5 # Point size used by plots
palette <- c("#cc660b", "#1f77b4", "#9467bd", "#2ca02c", # Color palette for MPs
             "#139BA6", "#7f7f7f") 

######Helper Functions######

# Generates a time series MSE figure (minus LCP or other reference lines)...
# proj_col is the name of the column of interest in the dataframe of projections proj_df, 
# hist_col is the analagous name from the historical dataframe hist_df...
# proj_col_lq and proj_col_uq are the names of the hist_df columns that contain the 
# lower and upper quantile values, respectively, to be plotted alongside proj_col
plot_ts <- function(proj_df = PlotDF, hist_df = hist, proj_col, hist_col, line_wid = lwd, 
                    q_alph = q_alpha, proj_col_lq, proj_col_uq, ylab, y_upper_lim = NA, 
                    text_size = font_size, show_historical = FALSE)
{
  # If show_historical == TRUE, create historical layer
  if(show_historical == FALSE)
  {
    historical_layer <- NULL
  } else
  {
    historical_layer <- geom_path(data = hist_df, aes(x = Year, y = .data[[hist_col]]), 
                                  size = line_wid)
  }
  
  # Generate the plot
  return(ggplot()  +
           geom_ribbon(data = proj_df, aes(x = Yr, ymin = .data[[proj_col_lq]], 
                                           ymax = .data[[proj_col_uq]], 
                                           fill = MP), alpha = q_alph) +
           geom_path(data = proj_df, aes(x = Yr, y = .data[[proj_col]], col = MP), 
                     size = line_wid) + historical_layer + 
           xlab("Year") + labs(col = "") + ylab(ylab) +
           theme_classic() + theme(text = element_text(size = text_size),
                                   legend.text = element_text(hjust = 0)) +
           scale_x_continuous(limits= c(1, 100), 
                              breaks = round(seq(10, 100, by = 10)),
                              labels = function(x) ifelse(x %% 20 == 0, x, "")) + 
           scale_y_continuous (expand = c (0, 0), limits = c (0, y_upper_lim)) +
           scale_color_manual(values = palette, labels = MP_labs) + 
           scale_fill_manual(values = palette, guide = "none"))
}

# Plots median projected SSB/SSBMSY, with projection variability... col = median
# projected SSB/SSBMSY, col_lq = lower quantile for SSB/SSBMSY, and col_uq = upper
# quantile for SSB/SSBMSY
plot_ssb_ssbmsy <- function(data = PlotDF, col, col_lq, col_uq, ssb_subscript, 
                            linewid = lwd, text_size = font_size, q_alph = q_alpha) {
  ssb_label <- bquote(SSB / SSB[.(ssb_subscript)])
  ggplot(data = data) +
    geom_ribbon(aes(x = Yr - 50, ymin = .data[[col_lq]], 
                    ymax = .data[[col_uq]], fill = MP), alpha = q_alph) + 
    geom_path(aes(x = Yr - 50, y = .data[[col]], col = MP), lwd = linewid) +
    geom_hline(yintercept = 0.4, linetype = "dashed", lwd = linewid) +
    geom_hline(yintercept = 0.8, linetype = "dashed", lwd = linewid) +
    ylab(ssb_label) + 
    xlab("Projection year") + 
    theme_classic() + labs(col = "", fill = "") +
    theme(text = element_text(size = text_size), legend.text = element_text(hjust = 0)) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) +
    scale_x_continuous(expand = c(0, 0), limits = c(1, NA)) +
    scale_color_manual(values = palette, labels = MP_labs) + 
    scale_fill_manual(values = palette, guide = "none")
}

# Calculates the recommended F at a given SSB, where the underlying HCR has a 
# prespecified LCP, UCP, and FCP (FMSY)
harvest_control_rule <- function(ssb, LCP, UCP, FCP) {
  slope <- FCP / (UCP - LCP)
  if (ssb < LCP) {
    return(0)
  } else if (ssb <= UCP) {
    return(slope * ssb - LCP * slope)
  } else {
    return(FCP)
  }
}

#######Plot Preparation/Setup######

# Create a list to store the plots
plot_list <- list()

# Truncate DF by capping number of years for projections
DFRP <- DF[DF$Yr<=maxyrp,] %>% mutate(I_Prob_Catch = if_else(C < 0.01, 1,0))

# Calculate short-term and long-term average annual variability in yield (AAVY) 
# for each sim-HCR-scenario combo... append annual catch difference to DF
AAVY_df <- data.frame(MP = character(), Sim = integer(), 
                      ST_AAVY = numeric(), LT_AAVY = numeric())
for (mp in MPs) {
  for (i in 1:OM@nsim) {
    catch <- DFRP$C[DFRP$MP == mp & DFRP$Sim == i]
    DFRP[DFRP$MP == mp & DFRP$Sim == i, "CatchDif"] <- c(NA, abs(diff(catch)))
    AAVY_value_ST <- mean(abs(diff(catch[1:4])))
    AAVY_value_LT <- mean(abs(diff(catch[25:49])))
    AAVY_df <- rbind(AAVY_df, data.frame(MP = mp, Sim = i, 
                                         ST_AAVY = AAVY_value_ST,
                                         LT_AAVY = AAVY_value_LT))
  }
}
# ... and then take the medians across sims, filter out "no fishing", and normalize
# to a 0-1 range, etc.
AAVY_df_med <- AAVY_df %>% group_by(MP) %>% 
  summarize(AAVY_ST = median(ST_AAVY), AAVY_LT = median(LT_AAVY)) %>% 
  filter(MP != "NFref") %>% 
  mutate_at(vars(contains("AAVY")), ~ (1-./max(.))+(1-max((1-./max(.))))) %>%
  pivot_longer(cols = c("AAVY_ST","AAVY_LT")) %>%
  mutate(term = ifelse(name == "AAVY_ST", "short","long"), variable = "AAVY") %>%
  select(!"name")

# Summarize DF for plotting
PlotDF <- DFRP %>%
  group_by(MP, Yr) %>%
  summarize(MedianSSB = median(SSB, na.rm = TRUE), 
            LQSSB = quantile(SSB, prob = LQ, na.rm = TRUE),
            UQSSB = quantile(SSB, prob = UQ, na.rm = TRUE),
            MedianCatch = median(C, na.rm = TRUE),
            LQCatch = quantile(C, prob = LQ, na.rm = TRUE),
            UQCatch = quantile(C, prob = UQ, na.rm = TRUE),
            MedianF = median(FM, na.rm = TRUE),
            LQF = quantile(FM, prob = LQ, na.rm = TRUE),
            UQF = quantile(FM, prob = UQ, na.rm = TRUE),
            AboveFirst5LRP = sum(I_Above_First5LRP)/OM@nsim,
            AboveLast5LRP = sum(I_Above_Last5LRP)/OM@nsim,
            AboveEntireTSLRP = sum(I_Above_EntireTSLRP/OM@nsim),
            MedianSSB_SSBMSYFirst5 = median(0.4*SSB/LRP_First5, na.rm = TRUE),
            LQSSB_SSBMSYFirst5 = quantile(0.4*SSB/LRP_First5, prob = LQ, na.rm = TRUE),
            UQSSB_SSBMSYFirst5 = quantile(0.4*SSB/LRP_First5, prob = UQ, na.rm = TRUE),
            MedianSSB_SSBMSYLast5 = median(0.4*SSB/LRP_Last5, na.rm = TRUE),
            LQSSB_SSBMSYLast5 = quantile(0.4*SSB/LRP_Last5, prob = LQ, na.rm = TRUE),
            UQSSB_SSBMSYLast5 = quantile(0.4*SSB/LRP_Last5, prob = UQ, na.rm = TRUE),
            MedianSSB_SSBMSYEntireTS = median(0.4*SSB/LRP_EntireTS, na.rm = TRUE),
            LQSSB_SSBMSYEntireTS = quantile(0.4*SSB/LRP_EntireTS, prob = LQ, na.rm = TRUE),
            UQSSB_SSBMSYEntireTS = quantile(0.4*SSB/LRP_EntireTS, prob = UQ, na.rm = TRUE),
            Prob_Catch = sum(I_Prob_Catch)/OM@nsim,
            MedianCDif = median(CatchDif, na.rm = TRUE), 
            LQCDif = quantile(CatchDif, prob = LQ, na.rm = TRUE), 
            UQCDif = quantile(CatchDif, prob = UQ, na.rm = TRUE)) %>%
  mutate(Yr = Yr + term_yr, MP = factor(MP, levels = levels(as.factor(DFRP$MP))))

# Add LRP_First5 to the hist dataframe
for(i in 1:OM@nsim)
{
  hist$LRP_First5[hist$Sim == i] <- DF$LRP_First5[DF$Sim == i]
}

# Format expanded historical sim data for plotting
hist <- hist %>%
  group_by(Year, MP) %>%
  summarize(MedianSSB = median(SSB), MedianCatch = median(CB), MedianF = median(FM),
            MedianSSB_SSBMSYFirst5 = median(0.4*SSB/LRP_First5))

# Add rows so the historical connects to the projected
PlotDF <- rbind(data.frame(MP = MPs, Yr = term_yr, 
                           MedianSSB = NA, 
                           LQSSB = NA,
                           UQSSB = NA,
                           MedianCatch = rep(rev(hist$MedianCatch)[1],MP_num),
                           LQCatch = NA,
                           UQCatch = NA,
                           MedianF = rep(rev(hist$MedianF)[1],MP_num),
                           LQF = NA,
                           UQF = NA,
                           AboveFirst5LRP = NA,
                           AboveLast5LRP = NA,
                           AboveEntireTSLRP = NA,
                           MedianSSB_SSBMSYFirst5 = NA,
                           LQSSB_SSBMSYFirst5 = NA,
                           UQSSB_SSBMSYFirst5 = NA,
                           MedianSSB_SSBMSYLast5 = NA,
                           LQSSB_SSBMSYLast5 = NA,
                           UQSSB_SSBMSYLast5 = NA,
                           MedianSSB_SSBMSYEntireTS = NA,
                           LQSSB_SSBMSYEntireTS = NA,
                           UQSSB_SSBMSYEntireTS = NA,
                           Prob_Catch = NA,
                           MedianCDif = NA,
                           LQCDif = NA,
                           UQCDif = NA), PlotDF)
hist <- rbind(hist, data.frame(Year = rep(term_yr + 1, MP_num), 
                               MedianSSB = sapply(MPs, function(lab) 
                                 PlotDF$MedianSSB[which(PlotDF$MP == lab&PlotDF$Yr == term_yr+1)]),
                               MedianSSB_SSBMSYFirst5 = sapply(MPs, function(lab) 
                                 PlotDF$MedianSSB_SSBMSYFirst5[which(PlotDF$MP == lab&PlotDF$Yr == term_yr+1)]),
                               MP = MPs)) 

# Calculate median LRPs/USRs across sims by year
RP_PLOTS <- DF %>%
  group_by(Yr) %>%
  summarize(across(starts_with(c("LRP", "USR")), median)) %>%
  pivot_longer(-Yr) %>%
  mutate(LRP = substr(name, 1, 3), 
         timeperiod = substr(name, 4, 5), 
         Yr = Yr + OM@nyears)
# Expand to get first 5 years LRP for historical years
RP_PLOTS <- rbind(RP_PLOTS[RP_PLOTS$timeperiod == "_F",], RP_PLOTS)
RP_PLOTS$Yr[1:(2*OM@nyears)] <- rep(1:50, each = 2)

# For radar plots - summarize long-term and short-term catch, SSB/SSBMSYfirst5
TRADEOFF_DF1 <- DFRP %>%
  filter(Yr <= 5 | (Yr > 25 & Yr <= 50)) %>%
  mutate(term = case_when(Yr <= 5 ~ "short", 
                          Yr > 25 & Yr <= 50 ~ "long"),
         MP = factor(MP, levels = levels(as.factor(DFRP$MP)))) %>%
  group_by(term, Sim, MP) %>%
  summarize(MedianCatch = median(C, na.rm = TRUE),
            MedianSSB_SSBMSY = median(0.4*SSB/LRP_First5, na.rm = TRUE)) %>%
  group_by(term, MP) %>%
  summarize(catch = median(MedianCatch, na.rm = TRUE),
            SSB_SSBMSY = median(MedianSSB_SSBMSY, na.rm = TRUE))

# For radar plots - summarize long-term and short-term P(SSB>first 5 years LRP) and
# P(catch = 0)
TRADEOFF_DF2 <- DFRP %>%
  filter(Yr <= 5 | (Yr > 25 & Yr <= 50)) %>%
  mutate(term = case_when(Yr <= 5 ~ "short", 
                          Yr > 25 & Yr <= 50 ~ "long"),
         MP = factor(MP, levels = levels(as.factor(DFRP$MP)))) %>%
  group_by(term, Yr, MP) %>%
  summarize(MedianLRP = median(sum(I_Above_First5LRP)/OM@nsim, na.rm = TRUE),
            MedianCatch0 = median(sum(I_Prob_Catch)/OM@nsim, na.rm = TRUE)) %>%
  group_by(term, MP) %>%
  summarize(LRP = median(MedianLRP, na.rm = TRUE),
            Catch0 = median(MedianCatch0, na.rm = TRUE))

# For radar plots - combine above two dataframes  plus AAVY, and ready for 
# plotting (by normalizing values to a 0-1 range, filtering out "no fishing", etc.) 
TRADEOFF_DF <- TRADEOFF_DF1 %>% left_join(TRADEOFF_DF2) %>% 
  filter(MP != "NFref") %>%
  mutate_at(vars(matches("^Catch$"), contains("SSB_SSBMSY")), ~ ./max(.)) %>%
  pivot_longer(cols = c("catch","SSB_SSBMSY","LRP", "Catch0"), 
               names_to = "variable", 
               values_to = "value") %>% rbind(AAVY_df_med) %>% 
  mutate(variable = factor(variable, 
                           levels = c("Catch0","SSB_SSBMSY","LRP","AAVY","catch")),
         value = ifelse(variable%in%c("Catch0")&is.na(value)==FALSE, 
                        1-value, value)) %>% 
  mutate(value = ifelse(is.na(value), 1, value)) %>%
  arrange(variable) %>% arrange(MP)

# Grab a subset of individual simulation runs to visualize later
set.seed(300)
sample <- sample(1:OM@nsim, size = 5)
DF_sample <- DF %>% filter(Sim%in%sample) %>%
  mutate(MP = factor(MP, levels = levels(as.factor(DF$MP))))

# Convert any Inf values to NA (i.e., ref points are 0)
hist[sapply(hist, is.infinite)] <- NA
PlotDF[sapply(PlotDF, is.infinite)] <- NA

######Harvest Control Rules Plot######

# Calculate the median control points across simulations for each year-HCR combo
fmsy <- do.call(c, lapply(CP_list[c("F_MSY_EntireTS","F_MSY_First5","F_MSY_Last5",
                                    "F_MSY_Last5Mixed", "F_MSY_Last5_MixedBPA")], 
                          function(x) apply(x, 1, median, na.rm = TRUE)))
lcp <- do.call(c, lapply(CP_list[c("SSB_MSY_EntireTS","SSB_MSY_First5","SSB_MSY_Last5",
                                   "SSB_MSY_First5Mixed", "SSB_MSY_First5_MixedBPA")], 
                         function(x) apply(0.4*x, 1, median, na.rm = TRUE)))
# A little more work is required for ucp due to BPA HCR
ucp <- do.call(c, lapply(CP_list[c("SSB_MSY_EntireTS","SSB_MSY_First5","SSB_MSY_Last5",
                                   "SSB_MSY_Last5Mixed")], 
                         function(x) apply(0.8*x, 1, median, na.rm = TRUE)))
mixed_BPA_ucp_matrix <- matrix(NA, nrow = nrow(CP_list[[1]]), ncol = ncol(CP_list[[1]]))
for(i in 1:nrow(CP_list[[1]]))
{
  for(j in 1:ncol(CP_list[[1]]))
  {
    mixed_BPA_ucp_matrix[i,j] <- 
      0.4*CP_list$SSB_MSY_First5_MixedBPA[i,j]*exp(1.645*CP_list$SSB_Sigma_MixedBPA[i,j])
  }
}
ucp <- c(ucp, apply(mixed_BPA_ucp_matrix, 1, median, na.rm = TRUE))

# Evaluate the resulting median HCRs across a range of SSB
HCR_DF <- data.frame(fmsy, lcp, ucp, year = 1:OM@proyears, 
                     MP = rep(as.character(MP_labs)[1:(length(MPs) - 1)], each = OM@proyears))
HCR_DF$MP <- factor(HCR_DF$MP, levels = as.character(MP_labs)[1:(length(MPs) - 1)])
SSB_values <- seq(0, round(max(ucp * 1.15)), by = 0.1)
HCR_DF <- HCR_DF[rep(seq_len(nrow(HCR_DF)), each = length(SSB_values)), ]
HCR_DF$SSB <- rep(SSB_values, times = nrow(HCR_DF) / length(SSB_values))
HCR_DF$FM <- mapply(harvest_control_rule, HCR_DF$SSB, HCR_DF$lcp, 
                    HCR_DF$ucp, HCR_DF$fmsy)

# Generate a plot of the median HCRs for each projection year
plot_list[["harvest_control_rules"]] <- 
  ggplot(data = HCR_DF[which(HCR_DF$year%%5==1),], 
       aes(x = SSB, y = FM, col = as.factor(year))) + geom_path(lwd = lwd) + 
  theme_classic() + scale_color_viridis_d() + facet_wrap(~MP, labeller = label_parsed) + 
  labs(y = "F", col = "Projection\nyear", x = "SSB (kt)") + 
  theme(text = element_text(size = font_size)) +
  scale_x_continuous(expand = c(0,0)) + scale_y_continuous(expand = c(0,0))
ggsave(paste0(getwd(),"/Figs/MSE/",OM_name,"/harvest_control_rules.jpeg"), plot=last_plot(), 
       width=width, height=height, units="in")

#######Radar Plots#######

# Generate radar plot illustrating LT trade-offs (exclude "no fishing")
LONG_TRADEOFF_DF <- TRADEOFF_DF %>% filter(term == "long")
hexagon_labels <- c("1-P(closure)","SSB/SSB[MSY]","P(SSB>LRP)","f(AAVY)","Catch")
hexagon_angles <- seq(0, 2 * pi, length.out = 6)[-6]+pi/10.5 # Exclude point to close polygon
hexagon_df <- data.frame(x = cos(hexagon_angles), y = sin(hexagon_angles),
                         label = hexagon_labels)
hexagon_df$label <- gsub("-", "*'-'*", hexagon_df$label)
LONG_TRADEOFF_DF$x <- rep(hexagon_df$x, MP_num-1)*LONG_TRADEOFF_DF$value
LONG_TRADEOFF_DF$y <- rep(hexagon_df$y, MP_num-1)*LONG_TRADEOFF_DF$value
plot_list[["Radar_Plot_LT"]] <- 
  ggplot(data = LONG_TRADEOFF_DF, aes(x = x, y = y, col = MP)) + 
  geom_polygon(data = hexagon_df, aes(x = x, y = y), 
               inherit.aes = F, fill = NA, col = "grey") + 
  geom_polygon(data = hexagon_df, aes(x = 0.75*x, y = 0.75*y), 
               inherit.aes = F, fill = NA, col = "grey") +
  geom_polygon(data = hexagon_df, aes(x = 0.5*x, y = 0.5*y), 
               inherit.aes = F, fill = NA, col = "grey") +
  geom_polygon(fill = NA, linewidth = lwd) +
  geom_text(data = hexagon_df[hexagon_df$x < (-0.9),], parse = TRUE,
            aes(x = 1.1*x, y = 1.1*y, label = label), 
            inherit.aes = F, size = 6, angle = 90, vjust = 0.5, hjust = 0.5) +
  geom_text(data = hexagon_df[hexagon_df$x > 0.9,], parse = TRUE,
            aes(x = 1.1*x, y = 1.1*y, label = label), 
            inherit.aes = F, size = 6, angle = -90, vjust = 0.5, hjust = 0.5) +
  geom_text(data = hexagon_df[!(hexagon_df$x > 0.9|hexagon_df$x < (-0.9)),], 
            parse = TRUE, aes(x = 1.1*x, y = 1.1*y, label = label), 
            inherit.aes = F, size = 6) +
  scale_color_manual(values = palette[-MP_num], 
                     labels = MP_labs[-(length(MP_labs))]) + 
  theme_classic() +
  theme(text = element_text(size = font_size),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        legend.position = "bottom") + 
  coord_fixed(clip = "off") + labs(col = "")
ggsave(paste0(getwd(),"/Figs/MSE/",OM_name,"/Radar_Plot_LT.jpeg"), plot=last_plot(), 
       width=width, height=height, units="in")

# Generate radar plot illustrating ST trade-offs (exclude "no fishing")
SHORT_TRADEOFF_DF <- TRADEOFF_DF %>% filter(term == "short")
hexagon_angles <- seq(0, 2 * pi, length.out = 6)[-6]+pi/10.5 # Exclude point to close polygon
hexagon_df <- data.frame(x = cos(hexagon_angles), y = sin(hexagon_angles),
                         label = hexagon_labels)
hexagon_df$label <- gsub("-", "*'-'*", hexagon_df$label)
SHORT_TRADEOFF_DF$x <- rep(hexagon_df$x, MP_num-1)*SHORT_TRADEOFF_DF$value
SHORT_TRADEOFF_DF$y <- rep(hexagon_df$y, MP_num-1)*SHORT_TRADEOFF_DF$value
plot_list[["Radar_Plot_ST"]] <- 
  ggplot(data = SHORT_TRADEOFF_DF, aes(x = x, y = y, col = MP)) + 
  geom_polygon(data = hexagon_df, aes(x = x, y = y), 
               inherit.aes = F, fill = NA, col = "grey") + 
  geom_polygon(data = hexagon_df, aes(x = 0.75*x, y = 0.75*y), 
               inherit.aes = F, fill = NA, col = "grey") +
  geom_polygon(data = hexagon_df, aes(x = 0.5*x, y = 0.5*y), 
               inherit.aes = F, fill = NA, col = "grey") +
  geom_polygon(fill = NA, linewidth = lwd) +
  geom_text(data = hexagon_df[hexagon_df$x < (-0.9),], parse = TRUE,
            aes(x = 1.1*x, y = 1.1*y, label = label), 
            inherit.aes = F, size = 6, angle = 90, vjust = 0.5, hjust = 0.5) +
  geom_text(data = hexagon_df[hexagon_df$x > 0.9,], parse = TRUE,
            aes(x = 1.1*x, y = 1.1*y, label = label), 
            inherit.aes = F, size = 6, angle = -90, vjust = 0.5, hjust = 0.5) +
  geom_text(data = hexagon_df[!(hexagon_df$x > 0.9|hexagon_df$x < (-0.9)),], 
            parse = TRUE, aes(x = 1.1*x, y = 1.1*y, label = label), 
            inherit.aes = F, size = 6) +
  scale_color_manual(values = palette[-MP_num], 
                     labels = MP_labs[-(length(MP_labs))]) + 
  theme_classic() +
  theme(text = element_text(size = font_size),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        legend.position = "bottom") + 
  coord_fixed(clip = "off") + labs(col = "")
ggsave(paste0(getwd(),"/Figs/MSE/",OM_name,"/Radar_Plot_ST.jpeg"), plot=last_plot(), 
       width=width, height=height, units="in")

#######Time Series Plots#######

# Median historical and projected SSB/(first 5 years SSBMSY)... include projection 
# variability
plot_list[["Median_SSBMSY"]] <- 
  plot_ts(proj_col = "MedianSSB_SSBMSYFirst5", hist_col = "MedianSSB_SSBMSYFirst5", 
        proj_col_lq = "LQSSB_SSBMSYFirst5", proj_col_uq = "UQSSB_SSBMSYFirst5", 
        ylab = bquote(SSB / SSB[.("MSY")]), show_historical = TRUE) + 
  geom_hline(yintercept = 0.4, linetype = "dashed", lwd = lwd) +
  geom_hline(yintercept = 0.8, linetype = "dashed", lwd = lwd)
ggsave(paste0(getwd(),"/Figs/MSE/",OM_name,"/Median_SSB_SSBMSYinitial.jpeg"), plot=last_plot(), 
       width=width, height=height, units="in")

# Median historical and projected SSB for first 3 HCRs... include projection 
# variability and 3 LRPs/USRs
PlotDFFilter <- PlotDF %>% filter(MP %in% c("MP_last5", "MP_first5", "MP_entirets"))
plot_list[["Median_SSB_LRPs"]] <- 
  plot_ts(proj_df = PlotDFFilter,
          proj_col = "MedianSSB", hist_col = "MedianSSB", proj_col_lq = "LQSSB", 
          proj_col_uq = "UQSSB", ylab = "SSB (kt)", show_historical = TRUE) +
  new_scale_color() + scale_color_manual(values = c("red","black"), guide = FALSE) +
  geom_line(data = RP_PLOTS[RP_PLOTS$Yr>50,], 
            aes(x = Yr, y = value, col = LRP, group = name), 
            linetype = 2, size = lwd*0.8, alpha = 0.6)
ggsave(paste0(getwd(),"/Figs/MSE/",OM_name,"/Median_SSB_LRPs.jpeg"), plot=last_plot(), 
       width=width, height=height, units="in")

# Median historical and projected SSB across the simulations... include projection 
# variability + median first 5 years LRP
plot_list[["Median_SSB"]] <- 
  plot_ts(proj_col = "MedianSSB", hist_col = "MedianSSB", proj_col_lq = "LQSSB", 
          proj_col_uq = "UQSSB", ylab = "SSB (kt)", show_historical = TRUE) +
  geom_hline(data = RP_PLOTS[RP_PLOTS$timeperiod == "_F"&RP_PLOTS$LRP == "LRP",], 
             aes(yintercept = value), linetype = 2, size = lwd)
ggsave(paste0(getwd(),"/Figs/MSE/",OM_name,"/Median_SSB.jpeg"), plot=last_plot(), 
       width=width, height=height, units="in")

# Hide "no fishing" from remaining plots
PlotDF <- PlotDF[which(PlotDF$MP != "NFref"),]
hist <- hist[which(hist$MP != "NFref"),]
palette <- palette[-MP_num]

# Median historical and projected F across the simulations + projection variability
plot_list[["Median_F"]] <- 
  plot_ts(proj_col = "MedianF", hist_col = "MedianF", show_historical = TRUE,
          proj_col_lq = "LQF", proj_col_uq = "UQF", ylab = "F")
ggsave(paste0(getwd(),"/Figs/MSE/",OM_name,"/Median_F.jpeg"), plot=last_plot(), 
       width=width, height=height, units="in")

# P(SSB > first 5 years LRP) for projection years ... dummies for hist_col, 
# proj_col_uq and proj_col_lq as they don't exist here
plot_list[["Prob_SSB_LRP"]] <- 
  plot_ts(proj_col = "AboveFirst5LRP", hist_col = "MedianSSB",
          proj_col_lq = "AboveFirst5LRP", proj_col_uq = "AboveFirst5LRP", 
          ylab = expression("P(SSB > LRP)"), y_upper_lim = 1)
ggsave(paste0(getwd(),"/Figs/MSE/",OM_name,"/Prob_SSB_LRP.jpeg"), plot=last_plot(), 
       width=width, height=height, units="in")

# Median historical and projected catch across the simulations + projection variability
plot_list[["Median_Catch"]] <- 
  plot_ts(proj_col = "MedianCatch", hist_col = "MedianCatch", show_historical = TRUE, 
          proj_col_lq = "LQCatch", proj_col_uq = "UQCatch", ylab = "Catch (kt)")
ggsave(paste0(getwd(),"/Figs/MSE/",OM_name,"/Median_Catch.jpeg"), plot=last_plot(), 
       width=width, height=height, units="in")

# P(Catch = 0 kt) for projection years ... dummies for hist_col, proj_col_uq and 
# proj_col_lq as they don't exist here
plot_list[["Prob_Catch"]] <- 
  plot_ts(proj_col = "Prob_Catch", hist_col = "MedianSSB",
          proj_col_lq = "Prob_Catch", proj_col_uq = "Prob_Catch", 
          ylab = expression("P(closure)"), y_upper_lim = 1)
ggsave(paste0(getwd(),"/Figs/MSE/",OM_name,"/Prob_Catch.jpeg"), plot=last_plot(), 
       width=width, height=height, units="in")

# AVY for projection years ... dummies for hist_col, proj_col_uq and proj_col_lq 
# as they don't exist here
plot_list[["AVY"]] <- 
  plot_ts(proj_col = "MedianCDif", hist_col = "Prob_Catch",
          proj_col_lq = "LQCDif", proj_col_uq = "UQCDif", 
          ylab = expression("AVY (kt)"))
ggsave(paste0(getwd(),"/Figs/MSE/",OM_name,"/AVY.jpeg"), plot=last_plot(), 
       width=width, height=height, units="in")

# Plot median projected SSB/(first 5 years SSBMSY), with projection variability
plot_list[["SSB_SSBMSY_First5"]] <- plot_ssb_ssbmsy(col = "MedianSSB_SSBMSYFirst5", 
                                                    col_lq = "LQSSB_SSBMSYFirst5",
                                                    col_uq = "UQSSB_SSBMSYFirst5", 
                                                    ssb_subscript = "MSYinitial")
ggsave(paste0(getwd(),"/Figs/MSE/",OM_name,"/SSB_SSBMSY_First5.jpeg"), plot=last_plot(), 
       width=width, height=height, units="in")

# Plot median projected SSB/(last 5 years SSBMSY), with projection variability
plot_list[["SSB_SSBMSY_Last5"]] <- plot_ssb_ssbmsy(col = "MedianSSB_SSBMSYLast5", 
                                                   col_lq = "LQSSB_SSBMSYLast5",
                                                   col_uq = "UQSSB_SSBMSYLast5", 
                                                   ssb_subscript = "MSYrecent")
ggsave(paste0(getwd(),"/Figs/MSE/",OM_name,"/SSB_SSBMSY_Last5.jpeg"), plot=last_plot(), 
       width=width, height=height, units="in")

# Plot median projected SSB/(entire time series SSBMSY), with projection variability
plot_list[["SSB_SSBMSY_EntireTS"]] <- plot_ssb_ssbmsy(col = "MedianSSB_SSBMSYEntireTS", 
                                                      col_lq = "LQSSB_SSBMSYEntireTS",
                                                      col_uq = "UQSSB_SSBMSYEntireTS", 
                                                      ssb_subscript = "MSYfull")
ggsave(paste0(getwd(),"/Figs/MSE/",OM_name,"/SSB_SSBMSY_EntireTS.jpeg"), plot=last_plot(), 
       width=width, height=height, units="in")

#######Simulation Runs Plots#######

# Plot the individual SSB trajectories from a sample of sim runs
plot_list[["Individual_Runs_SSB"]] <- ggplot(data = DF_sample) + 
  geom_line(aes(x = Yr, y = SSB, col = as.factor(Sim)), lwd = lwd) + 
  facet_wrap(~MP, ncol = 2) + theme_classic() + 
  theme(text = element_text(size = font_size)) + 
  scale_y_continuous(limits = c(0, NA), expand = c(0,0)) + 
  labs(col = "Simulation #", y = "SSB (kt)")
ggsave(paste0(getwd(),"/Figs/MSE/",OM_name,"/Individual_Runs_SSB.jpeg"), 
       plot=last_plot(), width=width, height=height, units="in")

# Plot the individual F trajectories from a sample of sim runs
plot_list[["Individual_Runs_F"]] <- ggplot(data = DF_sample) + 
  geom_line(aes(x = Yr, y = FM, col = as.factor(Sim)), lwd = lwd) + 
  facet_wrap(~MP, ncol = 2) + theme_classic() + 
  theme(text = element_text(size = font_size)) + 
  scale_y_continuous(limits = c(0, NA), expand = c(0,0)) + 
  labs(col = "Simulation #", y = "F")
ggsave(paste0(getwd(),"/Figs/MSE/",OM_name,"/Individual_Runs_F.jpeg"), 
       plot=last_plot(), width=width, height=height, units="in")

# Save list of plots so can combine with plots from other OM_names
saveRDS(plot_list, paste0(getwd(),"/Figs/MSE/",OM_name,"/",OM_name,"plots.rds"))