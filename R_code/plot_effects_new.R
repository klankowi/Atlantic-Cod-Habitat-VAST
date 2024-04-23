rm(list=ls())

# Load libraries
library(VAST)
library(effects)
library(tidyverse)
library(here)
library(splines)

# Set GGplot auto theme
theme_set(theme(panel.grid.major = element_line(color='lightgray'),
                panel.grid.minor = element_blank(),
                panel.background = element_blank(),
                panel.border = element_rect(color='black', linewidth=1, fill=NA),
                legend.position = "bottom",
                axis.text.x=element_text(size=12),
                axis.text.y=element_text(size=12),
                axis.title.x=element_text(size=14),
                axis.title.y=element_text(size=14, angle=90, vjust=2),
                plot.title=element_text(size=14, hjust = 0, vjust = 1.2),
                plot.caption=element_text(hjust=0, face='italic', size=12)))

# load data
load(here("VAST_runs/medium/Overall_BC/Overall_BC_medcod_allstrat_natsplin_fsON.Rdata"))

# Set WD for plotting
out_dir = here("VAST_runs/medium/Overall_BC/")

# Load functions
source(here("R_code/utilities/vast_functions.R"))
# Negate function
'%notin%' <- function(x,y)!('%in%'(x,y))

# Call parameter names
params <- colnames(scaled.covars)
params <- params[params %notin% c('Lon', 'Lat', 'Year')]

# Call category names
catnames <- fit$category_names
ncat = length(catnames)

  vast_covariate_effects<- get_vast_covariate_effects(vast_fit = fit, 
                                                 params_plot = c(params), 
                                                 params_plot_levels = 100, 
                                                 effects_pad_values = c(), 
                                                 nice_category_names = 'Medium Cod',
                                                 out_dir = out_dir,
                                                 category_to_use = 1,
                                                 ncat = ncat)
  
  # Remove intermediates
  rm(list=setdiff(ls(), c('vast_covariate_effects', 'fit', 'strata_use')))
  
  names_stay <- c("fit", "se", "lower", "upper", "Lin_pred")

  betternames=data.frame(
    Covariate = c('cobble_P', 'gravel_P', 'mud_P', 'sand_P', 'rugos', 'BATHY.DEPTH',
                  'h_bt', 'nao', 'amo'),
    Better = c('Cobble', 'Gravel', 'Mud', 'Sand', 'Rugosity', 'Depth (m)', 
               'Bottom Temp (C)', 'Monthly NAO', 'Annual AMO')
  )
  
  vast_cov_eff_l <- vast_covariate_effects %>%
    drop_na(Value)
  
  vast_cov_eff_l <- merge(vast_cov_eff_l, betternames)
  
  vast_cov_eff_l <- vast_cov_eff_l %>% 
    dplyr::select(-Covariate) %>% 
    rename(Covariate = Better)
  
  ylim_dat <- vast_cov_eff_l %>%
    group_by(., Lin_pred, Covariate) %>%
    summarize(., Min = min(lower, na.rm = TRUE), Max = max(upper, na.rm = TRUE))
  
  names_keep <- unique(vast_cov_eff_l$Covariate)
  
  colnames(fit$covariate_data)[4:12] <- c('Cobble', 'Gravel', 'Mud', 'Sand',
                                          'Rugosity', 'Depth (m)', 'Bottom Temp (C)',
                                          'Monthly NAO', 'Annual AMO')
  
  samp_dat <- fit$covariate_data %>% dplyr::select({
    {
      names_keep
    }
  }) %>% gather(., "Covariate", "Value")
  
  samp_quants <- samp_dat %>% 
    group_by(Covariate) %>% 
    mutate(
           q0.33=quantile(Value, 0.33),
           q0.66=quantile(Value, 0.66)) %>% 
    dplyr::select(-Value) %>% 
    unique() %>% 
    as.data.frame() %>% 
    pivot_longer(cols=c('q0.33', 'q0.66'))
  
  samp_quants <- samp_quants %>% 
    filter(Covariate != 'Cobble')
  
  samp_quants$name[samp_quants$Covariate == 'Sand'] <- 'q0.5'
  samp_quants$value[samp_quants$Covariate == 'Sand'] <- 
    quantile(samp_dat$Value[samp_dat$Covariate == 'Sand'], 0.5)
  
  samp_quants <- unique(samp_quants)
  
  vast_cov_eff_l$Covariate <- factor(vast_cov_eff_l$Covariate,
                                     levels = c('Cobble', 'Gravel', 'Sand', 'Mud',
                                                'Depth (m)', 'Bottom Temp (C)',
                                                'Rugosity', 'Monthly NAO', 
                                                'Annual AMO'))
  
  samp_quants$Covariate <- factor(samp_quants$Covariate,
                                  levels = c('Cobble', 'Gravel', 'Sand', 'Mud',
                                             'Depth (m)', 'Bottom Temp (C)',
                                             'Rugosity', 'Monthly NAO', 
                                             'Annual AMO'))
  
  samp_dat$Covariate <- factor(samp_dat$Covariate,
                               levels = c('Cobble', 'Gravel', 'Sand', 'Mud',
                                          'Depth (m)', 'Bottom Temp (C)',
                                          'Rugosity', 'Monthly NAO', 
                                          'Annual AMO'))
  
  plot_out <- ggplot() +
    geom_ribbon(data = vast_cov_eff_l[vast_cov_eff_l$Lin_pred == 'X1',], 
                aes(x = Value, ymin = lower, ymax = upper), 
                fill = "#bdbdbd", alpha=0.7) +
    geom_line(data = vast_cov_eff_l[vast_cov_eff_l$Lin_pred == 'X1',], aes(x = Value, y = fit)) +
    geom_vline(data=samp_quants, aes(xintercept=value),
               lty=2, col='darkgray', alpha=0.8) +
    xlab("Covariate value") +
    ylab("Linear predictor 1 fitted value") +
    facet_wrap(vars(Covariate), scales = "free_x", nrow=3) +
    theme(strip.text.x = element_text(size=12),
          strip.background = element_rect(fill='transparent',
                                          linewidth = 0),
          axis.text.x = element_text(size=11, angle=20))
  
  plot_out2 <- plot_out +
    geom_rug(data = samp_dat, aes(x = Value))

  plot_out2

  ggsave(here('VAST_runs/medium/Overall_BC/Medium_Cod_Effects_X1.png'),
         plot_out2,
         height=8.5, width = 10)  
  rm(plot_out, plot_out2)
  
  plot_out <- ggplot() +
    geom_ribbon(data = vast_cov_eff_l[vast_cov_eff_l$Lin_pred == 'X2',], 
                aes(x = Value, ymin = lower, ymax = upper), 
                fill = "#bdbdbd", alpha=0.7) +
    geom_line(data = vast_cov_eff_l[vast_cov_eff_l$Lin_pred == 'X2',], aes(x = Value, y = fit)) +
    geom_vline(data=samp_quants, aes(xintercept=value),
               lty=2, col='darkgray', alpha=0.8) +
    xlab("Covariate value") +
    ylab("Linear predictor 2 fitted value") +
    facet_wrap(vars(Covariate), scales = "free_x", nrow=3) +
    theme(strip.text.x = element_text(size=12),
          strip.background = element_rect(fill='transparent',
                                          linewidth = 0),
          axis.text.x = element_text(size=11, angle=20))
  
  plot_out2 <- plot_out +
    geom_rug(data = samp_dat, aes(x = Value))
  
  plot_out2
  
  ggsave(here('VAST_runs/medium/Overall_BC/Medium_Cod_Effects_X2.png'),
         plot_out2,
         height=8.5, width = 10)
  