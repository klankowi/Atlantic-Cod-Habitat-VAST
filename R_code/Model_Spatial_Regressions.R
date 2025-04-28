rm(list=ls())

library(here)
library(tidyverse)

sizes <- c('small', 'medium', 'large')

for(i in 1:length(sizes)){
  cog <- read.csv(paste0(here('VAST_runs'), '/',
                        sizes[i], 
                        '/Overall_BC/ALL_Catchability2/COG_New_Process.csv'))
  cog$Size <- paste0(str_to_sentence(sizes[i]))
  assign(paste0(sizes[i]), cog)
  rm(cog)
}

cog <- rbind(small, medium, large)

cog <- cog %>% 
  mutate(Stratum = Region) %>% 
  mutate(easting = easting / 1000,
         northing= northing / 1000) %>% 
  mutate(Stratum = factor(Stratum, levels = c('ALL', 'EGOM', 
                                              'GBK', 'SNE', 'WGOM'))) %>% 
  mutate(Season = factor(Season, levels = c('Spring', 'Fall'))) %>% 
  #mutate(Size = factor(Size, levels = c('Small', 'Medium', 'Large'))) %>% 
  dplyr::select(-X, -Region) %>% 
  pivot_longer(cols = c('northing', 'easting'),
               names_to = 'Direction',
               values_to = 'Estimate') %>% 
  #separate(Direction, into=c('Direction', 'Trash'), 
  #         sep = "[^[:alnum:]]+") %>% 
  #dplyr::select(-Trash) %>% 
  as.data.frame()

out <- cog %>% 
  group_by(Size, Season, Stratum, Direction) %>% 
  summarise(Change = 
              as.numeric(
                lm(Estimate ~ Year)$coefficients[2]
              )) %>% 
  as.data.frame()

out

  