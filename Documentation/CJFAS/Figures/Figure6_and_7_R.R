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
load(here("VAST_runs/medium/Overall_BC/ALL_Catchability2/Overall_BC_mediumcod_allstrat_natsplin_fsON_ALL_SVC.Rdata"))
medium <- fit
load(here("VAST_runs/small/Overall_BC/ALL_Catchability2/Overall_BC_smallcod_allstrat_natsplin_fsON_ALL_Catchability.Rdata"))
small <- fit
load(here("VAST_runs/large/Overall_BC/ALL_Catchability2/Overall_BC_largecod_allstrat_natsplin_fsON_ALL_Catchability_SVC.Rdata"))
large <- fit

# Clean workspace
rm(list=setdiff(ls(), c('small', 'medium', 'large')))

# Set WD for plotting
out_dir = here("VAST_runs/")

# Load functions
source(here("R_code/utilities/vast_functions.R"))
# Negate function
'%notin%' <- function(x,y)!('%in%'(x,y))

# Name size classes
sizes <- c('small', 'medium', 'large')

# Blank dataframe to fill
keep <- data.frame(
  fit=NA, se=NA, lower=NA, upper=NA, Lin_pred=NA, Value=NA,
  Covariate=NA, Size=NA
)

# Extract data, looping through sizes
for(i in 1:length(sizes)){
  size <- sizes[i]
  if(size == 'small'){fit <- small}
  if(size == 'medium'){fit <- medium}
  if(size == 'large'){fit <- large}
  
  # Call parameter names
  params <- colnames(fit$covariate_data)
  params <- params[params %notin% c('Lon', 'Lat', 'Year')]
  
  # Parameters not used in certain models
  if(size == 'small'){params <- params[params %notin% c('amo')]}
  if(size == 'large'){params <- params[params %notin% c('cobble_P',
                                                        'sand_P', 
                                                        'mud_P')]}
  # Call category names
  catnames <- fit$category_names
  ncat = length(catnames)
  
  # Extract effects
  vast_covariate_effects<- get_vast_covariate_effects(vast_fit = fit, 
                                                      params_plot = c(params), 
                                                      params_plot_levels = 100, 
                                                      effects_pad_values = c(), 
                                                      nice_category_names = paste0(size, ' Cod'),
                                                      out_dir = out_dir,
                                                      category_to_use = 1,
                                                      ncat = ncat)
  
  # Parameters we want to use
  names_stay <- c("fit", "se", "lower", "upper", "Lin_pred")
  
  # Merge non-informative names to informative name descriptions
  betternames=data.frame(
    Covariate = c('cobble_P', 'gravel_P', 'mud_P', 'sand_P', 'rugos', 'BATHY.DEPTH',
                  'h_bt', 'nao', 'amo'),
    Better = c('Cobble', 'Gravel', 'Mud', 'Sand', 'Rugosity', 'Depth (m)', 
               'BT (C)', 'NAO', 'AMO')
  )
  
  # Remove NA values
  vast_cov_eff_l <- vast_covariate_effects %>%
    drop_na(Value)
  
  # Merge to better names
  vast_cov_eff_l <- merge(vast_cov_eff_l, betternames)
  
  # Clean dataframe
  vast_cov_eff_l <- vast_cov_eff_l %>% 
    dplyr::select(-Covariate) %>% 
    rename(Covariate = Better)
  
  # Name size model
  vast_cov_eff_l$Size <- size
  
  # Bind to blank dataframe
  keep <- rbind(keep, vast_cov_eff_l)
  
  # Clean workspace
  rm(size, params, catnames, ncat, vast_covariate_effects, names_stay,
     betternames, vast_cov_eff_l)
  
}
# A raw copy (takes forever to reach this point, don't want to redo)
ejectbutton <- keep

# Remove blank initiation dataframe
keep <- keep[!is.na(keep$Size),]

# CLean
keep <- keep %>% 
  mutate(Size = str_to_sentence(Size)) %>% 
  mutate(Size = factor(Size, levels = c('Small', 'Medium', 'Large')),
         Covariate = factor(Covariate, levels=c('BT (C)',
                                                'Depth (m)',
                                                'Gravel',
                                                'Sand',
                                                'Mud',
                                                'Cobble',
                                                'Rugosity',
                                                'NAO',
                                                'AMO')))

# Function to tag facets without removing strips (based on egg::tag_facet)
tag_facet2 <- function(p, open = "(", close = ")", tag_pool = letters, x = -Inf, y = Inf, 
                       hjust = -0.5, vjust = 1.5, fontface = 2, family = "", ...) {
  
  gb <- ggplot_build(p)
  lay <- gb$layout$layout
  tags <- cbind(lay, label = paste0(open, tag_pool[lay$PANEL], close), x = x, y = y)
  p + geom_text(data = tags, aes_string(x = "x", y = "y", label = "label"), ..., hjust = hjust, 
                vjust = vjust, fontface = fontface, family = family, inherit.aes = FALSE)
}

# First linear predictor
fig6a <- ggplot(data=keep[keep$Lin_pred == 'X1',]) +
  geom_line(aes(x=Value, y=fit)) +
  geom_ribbon(aes(x=Value, ymin=lower, ymax=upper),
              alpha=0.3) +
  facet_grid2(c("Size", "Covariate"),
              scales='free_x', independent = 'none') +
  
  labs(x='Covariate value', y='Conditional effect on LP1') +
  
  scale_y_continuous(limits=c(-15, 6), expand=c(0,0)) +
  
  theme(strip.background = element_rect(fill='lightgray'),
        axis.text.x = element_text(angle= 90, size= 8),
        axis.text.y = element_text(size = 8))

ann_text <- data.frame(Value = c(-4.5, -4.5, -4.5, -4.5),
                       fit = c(0.0365, 0.626, 0.5, 0.5),
                       lab = rep('Cov.\nnot\nsig.',4),
                       Covariate = factor(c('AMO', 'Sand', 'Mud', 'Cobble'),
                                          levels = c('BT (C)',
                                                     'Depth (m)',
                                                     'Gravel',
                                                     'Sand',
                                                     'Mud',
                                                     'Cobble',
                                                     'Rugosity',
                                                     'NAO',
                                                     'AMO')),
                       Size = factor(c('Small', 'Large', 'Large', 'Large'),
                                     levels=c('Small', 'Medium', 'Large')))

fig6a <- fig6a + geom_text(
  data    = ann_text,
  mapping = aes(label = lab, x=fit, y=Value)
)

# Label facets
fig6a <- tag_facet2(fig6a, open='', close='', col='gray30', tag_pool = 
                      c(letters[1:8], '', letters[9:20], '', '', '', letters[21:23]))
fig6a

# Save local
ggsave(here('Documentation/CJFAS/Figures/Fig 6.pdf'),
       fig6a,
       height=23.7, width=18.2, units='cm',
       dpi = 600)

# Save to collaborators
ggsave("C:/Users/klankowicz/Box/Katie Lankowicz/Data_Analysis/Cod/Writings/CJFAS/Figures/Fig 6.pdf",
       fig6a,
       height=23.7, width=18.2, units='cm',
       dpi = 600)

# Second linear predictor
fig7a <- ggplot(data=keep[keep$Lin_pred == 'X2',]) +
  geom_line(aes(x=Value, y=fit)) +
  geom_ribbon(aes(x=Value, ymin=lower, ymax=upper),
              alpha=0.3) +
  facet_grid2(c("Size", "Covariate"),
              scales='free_x', independent = 'none') +
  
  labs(x='Covariate value', y='Conditional effect on LP2') +
  
  theme(strip.background = element_rect(fill='lightgray'),
        axis.text.x = element_text(angle= 90, size= 8),
        axis.text.y = element_text(size = 8))

ann_text <- data.frame(Value = c(-1.25, -1.25, -1.25, -1.25),
                       fit = c(0.0365, 0.626, 0.5, 0.5),
                       lab = rep('Cov.\nnot\nsig.',4),
                       Covariate = factor(c('AMO', 'Sand', 'Mud', 'Cobble'),
                                          levels = c('BT (C)',
                                                     'Depth (m)',
                                                     'Gravel',
                                                     'Sand',
                                                     'Mud',
                                                     'Cobble',
                                                     'Rugosity',
                                                     'NAO',
                                                     'AMO')),
                       Size = factor(c('Small', 'Large', 'Large', 'Large'),
                                     levels=c('Small', 'Medium', 'Large')))

fig7a <- fig7a + geom_text(
  data    = ann_text,
  mapping = aes(label = lab, x=fit, y=Value)
)

# Label facets
fig7a <- tag_facet2(fig7a, open='', close='', col='gray30', tag_pool = 
                      c(letters[1:8], '', letters[9:20], '', '', '', letters[21:23]))
fig7a

# Save local
ggsave(here('Documentation/CJFAS/Figures/Fig 7.pdf'),
       fig7a,
       height=23.7, width=18.2, units='cm',
       dpi = 600)

# Save to collaborators
ggsave("C:/Users/klankowicz/Box/Katie Lankowicz/Data_Analysis/Cod/Writings/CJFAS/Figures/Fig 7.pdf",
       fig7a,
       height=23.7, width=18.2, units='cm',
       dpi = 600)
