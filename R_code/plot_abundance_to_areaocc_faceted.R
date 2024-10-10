rm(list=ls())

library(VAST)
library(tidyverse)
library(here)
library(rstatix)

sizes <- c('small', 'medium', 'large')
strata <- c('ALL', 'EGOM', 'GBK', 'SNE', 'WGOM')

#### Center of gravity ####
cog <- data.frame(
  easting=NA, northing=NA, e.sd=NA, n.sd=NA,
  Year=NA, Season=NA, units=NA, strata=NA,
  Size=NA
)

outercog <- data.frame(
  easting=NA, northing=NA, e.sd=NA, n.sd=NA,
  Year=NA, Season=NA, units=NA, strata=NA,
  Size=NA
)

for(i in 1:length(sizes)){
  
  for(j in 1:length(strata)){
    innercog <- read.csv(here(
      paste0('VAST_runs/', sizes[i],
             '/Overall_BC/',
             strata[j],
             '/COG_', strata[j], '.csv')
    ))
    innercog$Size <- paste0(sizes[i])
    
    outercog <- rbind(outercog, innercog)
    rm(innercog)
  }
  cog <- rbind(cog, outercog)
}

cog <- cog %>% 
  filter(!is.na(Size)) %>% 
  dplyr::select(-units) %>% 
  unique() %>% 
  as.data.frame()
rm(outercog)

#### Relative abundance ####
ind <- data.frame(
  Size=NA, Year=NA, Season=NA, strata=NA, 
  RelAbun=NA, Std.Err=NA
)

for(i in 1:length(sizes)){
  innerind <- read.csv(here(
    paste0('VAST_runs/', sizes[i],
           '/Overall_BC/ALL/Index.csv')
    ))
  
  innerind <- innerind %>% 
    separate(Category, into=c('Size', 'crap')) %>% 
    separate(Time, into=c('Year', 'Season')) %>% 
    dplyr::select(-Units, -Std..Error.for.ln.Estimate., -crap) %>% 
    rename(strata = Stratum,
           RelAbun = Estimate,
           Std.Err = Std..Error.for.Estimate) %>% 
    mutate(Size=tolower(Size))
  ind <- rbind(ind, innerind)
  rm(innerind)

}
ind <- ind[!is.na(ind$Size),]

#### Area occupied ####
ao <- data.frame(
  AreaOcc=NA, St.Err.AO=NA, Year=NA, Season=NA, strata=NA,
  Size=NA
)

for(i in 1:length(sizes)){
  innerao <- read.csv(here(
    paste0('VAST_runs/', sizes[i],
           '/Overall_BC/ALL/AreaOcc.csv')
  ))
  
  innerao <- innerao %>% 
    mutate(Size=paste0(sizes[i])) %>% 
    dplyr::select(-units) %>% 
    rename(AreaOcc=area.occ,
           St.Err.AO = std.err)
  ao <- rbind(ao, innerao)
  rm(innerao)
  
}
ao <- ao[!is.na(ao$Size),]

df <- merge(cog, ind, by=c('Size', 'strata', 'Season', 'Year'))
df <- merge(df, ao, by=c('Size','strata', 'Season','Year'))
df$Size <- str_to_sentence(df$Size)

size <- c('Small','Medium','Large')
seasons <- c('Spring', 'Fall')

hold <- df

for(i in 1:length(seasons)){
  df <- hold[hold$Season == paste0(seasons[i]),]
  
  for(j in 1:length(size)){
    df2 <- df[df$Size == paste0(size[j]),]
    print(
      ggplot(data=df2) +
             geom_point(aes(x=RelAbun, y=AreaOcc, col=Year)) +
             facet_wrap(vars(strata), scales='free') +
             ggtitle(paste0(size[j], ' ', seasons[i]))
    )
    
    # print(
    # ggplot(data=df2) +
    #   geom_point(aes(x=northing, y=AreaOcc, col=Year)) +
    #   facet_wrap(vars(strata), scales='free') +
    #   ggtitle(paste0(size[j], ' ', seasons[i])) +
    #   labs(x='Northing', y='Area Occupied')
    # )
    test <- df2 %>% 
      group_by(Season, strata, Size) %>% 
      cor_test(AreaOcc, northing, method='spearman') %>% 
      as.data.frame() %>% 
      dplyr::select(Season, strata, Size, cor, p)
    test$sig[test$p<=0.05] <- '*'
    test$sig[test$p>0.05] <- ' '
    #message(paste0(size[j], ' ', seasons[i]))
    #print(test)
    
    rm(df2, test)
  }
}

