rm(list=ls())

library(here)
library(VAST)
library(tidyverse)
library(ggforce)
library(ggh4x)

# Negate function
'%notin%' <- function(x,y)!('%in%'(x,y))

# Set GGplot auto theme
theme_set(theme(panel.grid.major = element_line(color='lightgray'),
                panel.grid.minor = element_blank(),
                panel.background = element_blank(),
                panel.border = element_rect(color='black', linewidth=1, fill=NA),
                legend.position = "bottom",
                legend.background = element_rect(fill='transparent', colour = 'transparent'),
                axis.text.x=element_text(size=10),
                axis.text.y=element_text(size=10),
                axis.title.x=element_text(size=11),
                axis.title.y=element_text(size=11, angle=90, vjust=2),
                plot.title=element_text(size=12, hjust = 0, vjust = 1.2),
                plot.caption=element_text(hjust=0, face='italic', size=12)))

# Load data
load(here("VAST_runs/medium/Overall_BC/ALL/Overall_BC_mediumcod_allstrat_natsplin_fsON_ALL.Rdata"))
medium <- fit
load(here("VAST_runs/small/Overall_BC/ALL/Overall_BC_smallcod_allstrat_natsplin_fsON_ALL.Rdata"))
small <- fit
load(here("VAST_runs/large/Overall_BC/ALL/Overall_BC_largecod_allstrat_natsplin_fsON_ALL.Rdata"))
large <- fit

# Clean workspace
rm(list=setdiff(ls(), c('small', 'medium', 'large')))

# Function to tag facets without removing strips (based on egg::tag_facet)
tag_facet2 <- function(p, open = "(", close = ")", tag_pool = letters, x = -Inf, y = Inf, 
                       hjust = -0.5, vjust = 1.5, fontface = 2, family = "", ...) {
  
  gb <- ggplot_build(p)
  lay <- gb$layout$layout
  tags <- cbind(lay, label = paste0(open, tag_pool[lay$PANEL], close), x = x, y = y)
  p + geom_text(data = tags, aes_string(x = "x", y = "y", label = "label"), ..., hjust = hjust, 
                vjust = vjust, fontface = fontface, family = family, inherit.aes = FALSE)
}

# Name size classes
sizes <- c('small', 'medium', 'large')

# Blank dataframe to fill
keep <- data.frame(
  major_radius = NA,
  minor_radius = NA,
  angle = NA,
  pred = NA,
  Range1 = NA,
  Range2 = NA,
  Size = NA
)

for(i in 1:length(sizes)){
  
  if(sizes[i]=='small'){fit <- small}
  if(sizes[i]=='medium'){fit <- medium}
  if(sizes[i]=='large'){fit <- large}
  
  Obj = fit$tmb_list$Obj
  FileName = 'aniso' 
  ControlList = list(Width = 4, Height = 5, Res = 200, Units = "in")
  type = "ellipse"
  Report = fit$Report
  TmbData = Obj$env$data
  
  if (all(c("Options", "Options_vec") %in% names(TmbData))) {
    Options_vec = TmbData$Options_vec
    Options = TmbData$Options
  }
  if ("Options_list" %in% names(TmbData)) {
    Options_vec = TmbData$Options_list$Options_vec
    Options = TmbData$Options_list$Options
  }
  Map = Obj$env$map
  Params = Obj$env$last.par.best
  
  Eigen = eigen(Report$H)
  
  rss = function(V) sqrt(sum(V[1]^2 + V[2]^2))
  Major_1 = Minor_1 = Major_2 = Minor_2 = NA
  
  if (!("logkappa1" %in% names(Map)) || !is.na(Map$logkappa1)) {
    Major_1 = Eigen$vectors[, 1] * Eigen$values[1] * 
      Report$Range_raw1
    Minor_1 = Eigen$vectors[, 2] * Eigen$values[2] * 
      Report$Range_raw1
  }
  if (!("logkappa2" %in% names(Map)) || !is.na(Map$logkappa2)) {
    Major_2 = Eigen$vectors[, 1] * Eigen$values[1] * 
      Report$Range_raw2
    Minor_2 = Eigen$vectors[, 2] * Eigen$values[2] * 
      Report$Range_raw2
  }
  
  Range = 1.1 * c(-1, 1) * max(abs(cbind(Major_1, Minor_1, 
                                         Major_2, Minor_2)), na.rm = TRUE)
  
  el <- data.frame(
    major_radius = c(rss(Major_1), rss(Major_2)),
    minor_radius = c(rss(Minor_1), rss(Minor_2)),
    angle = c(90-(-1 * (atan(Major_1[1]/Major_1[2])/(2 * pi) * 360 - 90)),
              90-(-1 * (atan(Major_2[1]/Major_2[2])/(2 * pi) * 360 - 90))),
    pred = factor(c(1, 2)),
    Range1 = Range[1],
    Range2 = Range[2],
    Size = str_to_sentence(sizes[i])
  )
  
  keep <- rbind(keep, el)
  
  rm(el, Range, Map, Obj, Report, TmbData, ControlList, Eigen,
     Major_1, Minor_1, Major_2, Minor_2, Options, Options_vec,
     Params, type, FileName)
}

keep <- keep[!is.na(keep$major_radius),]

keep <- keep %>% 
  mutate(Size = factor(Size, levels=c('Small', 'Medium', 'Large')))

eli <- vector('list', length = nrow(keep))

region <- st_read(here('Data/GIS/cod_region_utm.shp'), quiet=T)
region <- st_make_valid(st_transform(region, "EPSG:32619"))

center <- sfheaders::sf_to_df(st_centroid(region))

for(i in 1:nrow(keep)){
  test <- ellipse(x = center$x, y = center$y,
                  sx = keep$major_radius[i]*1000, 
                  sy = keep$minor_radius[i]*1000,
                  rotation = keep$angle[i], 
                  n=100)
  
  test <- st_as_sf(as.data.frame(test), coords=c('x', 'y'),
                   crs="EPSG:32619")
  
  test <- test %>% 
    mutate(pred = keep$pred[i]) %>% 
    group_by(pred) %>% 
    summarise(geometry = st_combine(geometry)) %>% 
    st_cast('POLYGON') %>% 
    mutate(Size = keep$Size[i])
  
  eli[[i]] <- test
  
  rm(test)
  
}
eli <- do.call(rbind, eli)

coast <- st_transform(ecodata::coast, st_crs(eli))

supfig4 <- ggplot(data=eli) + 
  geom_sf(data=coast) +
  geom_sf(data=region, fill=NA) +
  
  geom_sf(aes(fill=pred, col=pred), alpha=0.7) +
  scale_fill_manual(values = c("#FFC20A", "#0C7BDC")) +
  scale_color_manual(values = c("#FFC20A", "#0C7BDC")) +
  
  coord_sf(xlim=c(-76, -65),
           ylim=c(36.75, 45),
           crs="EPSG:4326") +
  
  facet_grid2(cols=vars(Size)) +
  labs(col='Linear\npredictor', fill='Linear\npredictor',
       x=' ', y=' ') +
  theme(axis.text.x = element_text(angle = 20,
                                   vjust=0.75,
                                   hjust=0.75),
        legend.position = 'right',
        legend.box.margin = margin(-10, -5, -10, -10),
        legend.box.spacing = margin(0,10,0,0))

supfig4 <- tag_facet2(supfig4, open='', close='')
supfig4

# Save local
ggsave(here('Documentation/Figures/Supplementary/Supp Fig 4.pdf'),
       supfig4,
       height=(18.2 * 0.4), width=18.2, units='cm',
       dpi = 600)

# Save to collaborators
ggsave("C:/Users/klankowicz/Box/Katie Lankowicz/Data_Analysis/Cod/Writings/CJFAS/Figures/Supplementary/Supp Fig 4.pdf",
       supfig4,
       height=(18.2 * 0.4), width=18.2, units='cm',
       dpi = 600)
