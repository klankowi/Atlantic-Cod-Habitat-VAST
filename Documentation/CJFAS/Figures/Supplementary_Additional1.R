# Compare catches of BTS and BLLS #
#### Workspace setup ####
# Clear workspace
rm(list=ls())

# Load packages
suppressPackageStartupMessages(library(here))
suppressPackageStartupMessages(library(tidyverse, quietly=T, verbose=F))
suppressPackageStartupMessages(library(sf, quietly=T, verbose=F))

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

# Negate function
'%notin%' <- function(x,y)!('%in%'(x,y))

#### Add sample information and covars ####
# Load data
surveys <- read.csv(here("Data/Survey_Data/Bio_Data_Agesep.csv"))
surveys$month <- month(as.POSIXct(surveys$DATE, format='%Y-%m-%d'))

# Only larges
#surveys <- surveys[surveys$AGEGROUP == 'Age5+',]

# Add seasonality
surveys$season[surveys$month %in% c(1,2,9,10,11,12)] <- 'BFall'
surveys$season[surveys$month %in% c(3,4,5,6,7,8)] <- 'ASpring'

surveys$YEAR[surveys$month %in% c(1,2)] <- 
  surveys$YEAR[surveys$month %in% c(1,2)]-1

surveys$SEASON <- paste0(surveys$YEAR, ' ', surveys$season)

# Fix stupid error
surveys$TIME <- as.numeric(as.factor(paste0(surveys$YEAR, ' ',
                                            surveys$SEASON)))

# Remove points not in model domain (on land)
region_shape<- st_read(here("Data/GIS/cod_region_UTM.shp"), quiet=T)
region_shape <- st_transform(region_shape, "EPSG:4326")
region_shape <- st_make_valid(region_shape)
region_shape <- dplyr::select(region_shape, OBJECTID, geometry)

surveys_sf <- st_as_sf(surveys, coords=c('LON', 'LAT'))
st_crs(surveys_sf) <- 'EPSG:4326'

surveys_sf <- st_intersection(surveys_sf, region_shape)

surveys <- surveys[surveys$HAUL_ID %in% surveys_sf$HAUL_ID,]

# Test for missing values
surveys <- surveys[!is.na(surveys$rugos) &
                     !is.na(surveys$nao) &
                     !is.na(surveys$h_bt),]

surveys <- unique(surveys)

# Trim
surveys <- surveys %>% 
  dplyr::select(SURVEY, COD_N, AGEGROUP, AGE_N, SEASON) %>% 
  filter(SURVEY %in% c('NEFSC BTS', 'NEFSC BLLS')) %>% 
  separate(SEASON, into=c('YEAR', 'SEASON'))

surveys$SEASON[surveys$SEASON == 'ASpring'] <- 'Spring'
surveys$SEASON[surveys$SEASON == 'BFall'] <- 'Fall'

surveys$AGEGROUP[surveys$AGEGROUP == 'Age0-2'] <- 'Small'
surveys$AGEGROUP[surveys$AGEGROUP == 'Age2-5'] <- 'Medium'
surveys$AGEGROUP[surveys$AGEGROUP == 'Age5+'] <- 'Large'

surveys <- surveys %>%
  filter(YEAR >=2014) %>% 
  mutate(AGEGROUP = factor(AGEGROUP, levels=c('Small', 'Medium', 'Large')),
         SEASON = factor(SEASON, levels=c('Spring', 'Fall'))) %>% 
  mutate(LOGAGE_N = log(AGE_N))

outliers <- surveys %>% 
  group_by(AGEGROUP, SURVEY, SEASON) %>% 
  mutate(outlier = ifelse(AGE_N > median(AGE_N) + IQR(AGE_N) * 1.5, TRUE , FALSE)) %>% 
  filter(outlier)

fig <- ggplot() + 
  geom_boxplot(data=surveys,
               aes(x=AGEGROUP, y=AGE_N, col=SURVEY),
               outlier.shape = NA) +
  geom_jitter(data=outliers[outliers$SURVEY == 'NEFSC BLLS',],
             aes(x=(as.numeric(AGEGROUP)-0.2), y=AGE_N, fill=SURVEY),
             alpha=0.5, pch=21, stroke=NA,
             width = 0.15, height=0) + 
  geom_jitter(data=outliers[outliers$SURVEY == 'NEFSC BTS',],
              aes(x=(as.numeric(AGEGROUP)+0.2), y=AGE_N, fill=SURVEY),
              alpha=0.5, pch=21, stroke=NA,
              width = 0.15, height=0) + 
  facet_wrap(vars(SEASON)) +
  labs(x='Size class', y='Catch (n)', col='Survey', fill='Survey')

# Save local
ggsave(here('Documentation/CJFAS/Figures/Fig SA1.pdf'),
       fig,
       height=(18.2 * 0.75), width=18.2, units='cm',
       dpi = 600)

# Save to collaborators
ggsave("C:/Users/klankowicz/Box/Katie Lankowicz/Data_Analysis/Cod/Writings/CJFAS/Figures/Fig SA1.pdf",
       fig,
       height=(18.2 * 0.75), width=18.2, units='cm',
       dpi = 600)
