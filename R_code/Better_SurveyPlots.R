# Clear workspace
rm(list=ls())

#### Load libraries ####
library(tidyverse, quietly=T,verbose=F)
library(readxl, quietly=T,verbose=F)
library(reshape2, quietly=T, verbose=F)
library(ggh4x, quietly=T, verbose=F)
library(zoo, quietly=T, verbose=F)
library(RColorBrewer, quietly=T, verbose=F)
library(gridExtra, quietly=T, verbose=F)
library(here)

# Set seed for reproducibility
set.seed(123)

# Add function
scaleFUN <- function(x) sprintf("%.1f", x)

#### Load data ####
# Read in index data
indd <- read.csv(here('Data/Survey_Data/Index_Data.csv'))
# Remove extra DFO series
indd <- subset(indd, AREA != "5Z1-2")
indd <- subset(indd, AREA != "5ZJ-M")

# Load in length-frequency data
lenfq <- read.csv(here('Data/Survey_Data/LenFreq_Data.csv'))
# Remove extra DFO series
lenfq <- subset(lenfq, AREA !='5Z1-2')
lenfq <- subset(lenfq, AREA !='5ZJ-M')

# Load in numbers-at-age data
numage <- read.csv(here('Data/Survey_Data/NumAge_Data.csv'))
# Remove extra DFO series
numage <- subset(numage, AREA !='5Z1-2')
numage <- subset(numage, AREA !='5ZJ-M')

# Load in survey data
survd <- read.csv(here('Data/Survey_Data/Survey_Data.csv'))

# Load in metadata
metad <- read.csv(here('Data/Survey_Data/Metadata.csv'))

# Set seasonal color palettes
Spring <- "#3A836E"
Summer <- "#643187"
Fall <- "#EF3B2C"
Winter <- "#0476D0"

#### Index time series ####
# Set GGplot auto theme
theme_set(theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.background = element_blank(),
                
                axis.line = element_line(colour = "black"),
                legend.position = "n",
                
                axis.text.x=element_text(size=20),
                axis.text.y=element_text(size=20),
                axis.title.x=element_text(size=22),
                axis.title.y=element_text(size=22, angle=90, vjust=2),
                legend.text=element_text(size=14),
                
                plot.title=element_text(size=20, hjust = 0, vjust = 1.2),
                
                strip.text.x = element_text(size = 20)))

# Remove instances of 1 sample
indd <- subset(indd, INDEX_NAME !="SMAST Video Trawl_STELLWAGEN_FALL_A1+")

# Indices
for(f in 1:length(unique(indd$SURVEY))){
  survo <- unique(indd$SURVEY)[f]
  survo.sub <- subset(indd, SURVEY == paste0(survo))
  survo.list <- split(survo.sub, f=survo.sub$STOCK)
  
  for(r in 1:length(survo.list)){
    test.list <- survo.list[[r]]
    ages.tab <- table(test.list$AGES)
    
    stocko <- unique(test.list$STOCK)
    
    if(length(ages.tab)>1){
      test.list$NAME <- paste0(test.list$AREA, "_",
                               test.list$SEASON, "_",
                               test.list$AGES)
    }
    if(length(ages.tab)==1){
      test.list$NAME <- paste0(test.list$AREA, "_",
                               test.list$SEASON)
    }
    
    test.list$NAME <- gsub("_", " ", test.list$NAME)
    table(test.list$NAME) 
    
    test.list <- test.list[with(test.list, order(NAME,YEAR)),]
    
    for(i in 1:nrow(test.list)){
      if(test.list$SEASON[i] == "SPRING"){
        test.list$COLOR[i] <- Spring
      }
      if(test.list$SEASON[i] == "SUMMER"){
        test.list$COLOR[i] <- Summer
      }
      if(test.list$SEASON[i] == "FALL"){
        test.list$COLOR[i] <- Fall
      }
      if(test.list$SEASON[i] == "WINTER"){
        test.list$COLOR[i] <- Winter
      }
    }
    
    if(unique(test.list$SURVEY) == 'ASMFC Shrimp Trawl'){
      keepers <- test.list
      keepers$COLOR[keepers$YEAR != 2021] <- NA
      
      p <- ggplot(test.list, aes(x=YEAR, y=INDEX_NO)) +
        geom_ribbon(aes(ymin=CI_LO_NO, ymax=CI_HI_NO), fill=test.list$COLOR, 
                    alpha=0.2, linetype=0) +
        geom_line(lwd=1, col=test.list$COLOR) +
        geom_point(data=keepers,
                   aes(x=YEAR, y=INDEX_NO, col=COLOR)) +
        scale_color_manual(na.value='transparent',
                           values = rev(unique(test.list$COLOR))) +
        labs(x="Year",
             y="Index (No/ tow)",
             title=paste0(survo)) +
        scale_y_continuous(labels=scaleFUN) +
        facet_wrap2(vars(test.list$NAME), scales="free_y", ncol=2, nrow=2)
      
      p <- p+theme(axis.title.x = element_blank())
      
      q <- ggplot(test.list, aes(x=YEAR, y=INDEX_KG)) + 
        geom_ribbon(aes(ymin=CI_LO_KG, ymax=CI_HI_KG), fill=test.list$COLOR, 
                    alpha=0.2, linetype=0) +
        geom_line(lwd=1, col=test.list$COLOR) +
        labs(x="Year",
             y="Index (Kg/ tow)") +
        geom_point(data=test.list[test.list$YEAR == 2021,],
                   aes(x=YEAR, y=INDEX_KG),
                   col=test.list$COLOR[test.list$YEAR == 2021]) +
        scale_y_continuous(labels=scaleFUN) +
        theme(strip.background = element_blank(), strip.text.x = element_blank())
      
      tmp <- arrangeGrob((p),(q),ncol=1)
      ggsave(paste0('C:/Users/klankowicz/Desktop/Index_New_Plots/',
                    survo, " ", stocko,
                    '.png'),
             tmp)
      
      rm(p,q,tmp)
    }
    
    if(unique(test.list$SURVEY) == 'DFO Trawl'){

      p <- ggplot(test.list, aes(x=YEAR, y=INDEX_NO)) +
        geom_ribbon(aes(ymin=CI_LO_NO, ymax=CI_HI_NO), fill=test.list$COLOR, 
                    alpha=0.2, linetype=0) +
        geom_line(lwd=1, col=test.list$COLOR) +
        scale_color_manual(na.value='transparent',
                           values = rev(unique(test.list$COLOR))) +
        labs(x="Year",
             y="Population size",
             title=paste0(survo)) 
      
      ggsave(paste0('C:/Users/klankowicz/Desktop/Index_New_Plots/',
                    survo, " ", stocko,
                    '.png'),
             p,
             width = 3012, height = 1020, units='px')
      
      rm(p)
    }
    
    if(unique(test.list$SURVEY) == 'MADMF Industry'){
      
      keepers <- test.list
      keepers$COLOR[keepers$YEAR != 2021] <- NA
      
      p <- ggplot(test.list, aes(x=YEAR, y=INDEX_NO)) +
        geom_ribbon(aes(ymin=CI_LO_NO, ymax=CI_HI_NO), fill=test.list$COLOR, 
                    alpha=0.2, linetype=0) +
        geom_line(lwd=1, col=test.list$COLOR) +
        geom_point(data=keepers,
                   aes(x=YEAR, y=INDEX_NO, col=COLOR)) +
        scale_color_manual(na.value='transparent',
                           values = rev(unique(test.list$COLOR))) +
        labs(x="Year",
             y="Index (No/ tow)",
             title=paste0(survo)) +
        scale_y_continuous(labels=scaleFUN) +
        facet_wrap2(vars(test.list$NAME), scales="free_y", ncol=2, nrow=2)
      
      ggsave(paste0('C:/Users/klankowicz/Desktop/Index_New_Plots/',
                    survo, " ", stocko,
                    '.png'),
             p,
             width = 3012, height = 1020, units='px')
      
      rm(p)
    }
    
    if(unique(test.list$SURVEY) == 'MADMF Inshore Trawl' &
       unique(test.list$STOCK) == 'SNE'){
      
      keepers <- test.list
      keepers$COLOR[keepers$YEAR != 2021] <- NA
      
      p <- ggplot(test.list, aes(x=YEAR, y=INDEX_NO)) +
        geom_ribbon(aes(ymin=CI_LO_NO, ymax=CI_HI_NO), fill=test.list$COLOR, 
                    alpha=0.2, linetype=0) +
        geom_line(lwd=1, col=test.list$COLOR) +
        geom_point(data=test.list[test.list$YEAR == 2021,],
                   aes(x=YEAR, y=INDEX_NO), col=Spring) +
        labs(x="Year",
             y="Index (No/ tow)",
             title=paste0(survo)) +
        scale_y_continuous(labels=scaleFUN)
      
      ggsave(paste0('C:/Users/klankowicz/Desktop/Index_New_Plots/',
                    survo, " ", stocko,
                    '.png'),
             p,
             width = 3012, height = 1020, units='px')
      
      rm(p)
    }
    
    if(unique(test.list$SURVEY) == 'MADMF Inshore Trawl' &
       unique(test.list$STOCK) == 'WGOM'){
      
      keepers <- test.list
      keepers$COLOR[keepers$YEAR != 2021] <- NA
      
      test.list$SEASONAGE <- paste0(test.list$SEASON, ' ', test.list$AGES)
      
      p <- ggplot(test.list, aes(x=YEAR, y=INDEX_NO)) +
        geom_ribbon(aes(ymin=CI_LO_NO, ymax=CI_HI_NO), fill=test.list$COLOR, 
                    alpha=0.2, linetype=0) +
        geom_line(lwd=1, col=test.list$COLOR) +
        geom_point(data=test.list[test.list$YEAR == 2021,],
                   aes(x=YEAR, y=INDEX_NO, col=rev(COLOR))) +
        labs(x="Year",
             y="Index (No/ tow)",
             title=paste0(survo)) +
        scale_y_continuous(labels=scaleFUN) +
        facet_wrap2(vars(SEASONAGE), scales='free_y')
      
      ggsave(paste0('C:/Users/klankowicz/Desktop/Index_New_Plots/',
                    survo, " ", stocko,
                    '.png'),
             p,
             width = 3012*2, height = 1020*2, units='px')
      
      rm(p)
    }
    
    if(unique(test.list$SURVEY) == 'ME-NH Inshore Trawl' &
       unique(test.list$STOCK) == 'EGOM'){
      
      keepers <- test.list
      keepers$COLOR[keepers$YEAR != 2021] <- NA
      
      test.list$SEASONAGE <- paste0(test.list$SEASON, ' ', test.list$AGES)
      
      p <- ggplot(test.list, aes(x=YEAR, y=INDEX_NO)) +
        geom_ribbon(aes(ymin=CI_LO_NO, ymax=CI_HI_NO), fill=test.list$COLOR, 
                    alpha=0.2, linetype=0) +
        geom_line(lwd=1, col=test.list$COLOR) +
        geom_point(data=test.list[test.list$YEAR == 2021,],
                   aes(x=YEAR, y=INDEX_NO, col=rev(COLOR))) +
        labs(x="Year",
             y="Index (No/ tow)",
             title=paste0(survo)) +
        scale_y_continuous(labels=scaleFUN) +
        facet_wrap2(vars(SEASONAGE), scales='free_y')
      
      q <- ggplot(test.list, aes(x=YEAR, y=INDEX_KG)) + 
        geom_ribbon(aes(ymin=CI_LO_KG, ymax=CI_HI_KG), fill=test.list$COLOR, 
                    alpha=0.2, linetype=0) +
        geom_line(lwd=1, col=test.list$COLOR) +
        labs(x="Year",
             y="Index (Kg/ tow)") +
        geom_point(data=test.list[test.list$YEAR == 2021,],
                   aes(x=YEAR, y=INDEX_KG),
                   col=test.list$COLOR[test.list$YEAR == 2021]) +
        facet_wrap2(vars(SEASONAGE), scales='free_y')
        scale_y_continuous(labels=scaleFUN) +
        theme(strip.background = element_blank(), strip.text.x = element_blank())
      
      tmp <- arrangeGrob((p),(q),ncol=1)
      ggsave(paste0('C:/Users/klankowicz/Desktop/Index_New_Plots/',
                    survo, " ", stocko,
                    '.png'),
             tmp,
             width = 3012*2, height = 1020*2, units='px'
             )
      
      rm(p,q,tmp)
    }
    
    if(unique(test.list$SURVEY) == 'ME-NH Inshore Trawl' &
       unique(test.list$STOCK) == 'WGOM'){
      
      keepers <- test.list
      keepers$COLOR[keepers$YEAR != 2021] <- NA
      
      test.list$SEASONAGE <- paste0(test.list$SEASON, ' ', test.list$AGES)
      
      p <- ggplot(test.list, aes(x=YEAR, y=INDEX_NO)) +
        geom_ribbon(aes(ymin=CI_LO_NO, ymax=CI_HI_NO), fill=test.list$COLOR, 
                    alpha=0.2, linetype=0) +
        geom_line(lwd=1, col=test.list$COLOR) +
        geom_point(data=test.list[test.list$YEAR == 2021,],
                   aes(x=YEAR, y=INDEX_NO, col=rev(COLOR))) +
        labs(x="Year",
             y="Index (No/ tow)",
             title=paste0(survo)) +
        scale_y_continuous(labels=scaleFUN) +
        facet_wrap2(vars(SEASONAGE), scales='free_y')
      
      q <- ggplot(test.list, aes(x=YEAR, y=INDEX_KG)) + 
        geom_ribbon(aes(ymin=CI_LO_KG, ymax=CI_HI_KG), fill=test.list$COLOR, 
                    alpha=0.2, linetype=0) +
        geom_line(lwd=1, col=test.list$COLOR) +
        labs(x="Year",
             y="Index (Kg/ tow)") +
        geom_point(data=test.list[test.list$YEAR == 2021,],
                   aes(x=YEAR, y=INDEX_KG),
                   col=test.list$COLOR[test.list$YEAR == 2021]) +
        facet_wrap2(vars(SEASONAGE), scales='free_y') +
        scale_y_continuous(labels=scaleFUN) +
        theme(strip.background = element_blank(), strip.text.x = element_blank())
      
      tmp <- arrangeGrob((p),(q),ncol=1)
      ggsave(paste0('C:/Users/klankowicz/Desktop/Index_New_Plots/',
                    survo, " ", stocko,
                    '.png'),
             tmp,
             width = 3012*2, height = 1020*2, units='px'
      )
      
      rm(p,q,tmp)
    }
    
    
    # Base plot
    p <- ggplot(test.list, aes(x=YEAR, y=INDEX_NO)) + 
      geom_ribbon(aes(ymin=CI_LO_NO, ymax=CI_HI_NO), fill=test.list$COLOR, 
                  alpha=0.2, linetype=0) +
      geom_line(lwd=1, col=test.list$COLOR) +
      labs(x="Year",
           y="Index (No/ tow)", 
           title=paste0(survo)) +
      #scale_color_manual(values = linecolors) +
      #scale_fill_manual(values=linecolors) +
      scale_y_continuous(labels=scaleFUN) +
      facet_wrap2(vars(test.list$NAME), scales="free_y", ncol=2, nrow=2)
    
    if(length(unique(test.list$NAME)) >1){
      p <- ggplot(test.list, aes(x=YEAR, y=INDEX_NO)) + 
        geom_ribbon(aes(ymin=CI_LO_NO, ymax=CI_HI_NO), fill=test.list$COLOR, 
                    alpha=0.2, linetype=0) +
        geom_line(lwd=1, col=test.list$COLOR) +
        labs(x="Year",
             y="Index (No/ tow)", 
             title=paste0(survo)) +
        #scale_color_manual(values = linecolors) +
        #scale_fill_manual(values=linecolors) +
        scale_y_continuous(labels=scaleFUN) +
        facet_wrap2(vars(test.list$NAME), scales="free_y", ncol=2, nrow=2)
    }
    
    if(f==3){
      ggsave(paste0('C:/Users/klankowicz/Desktop/Index_New_Plots/',
                    survo, " ", stocko,
                    '.png'),
             p)
      rm(p)
      
      next()
    }
    
    truthfulness <- max(test.list$YEAR) == 2021 &
      is.na(test.list$INDEX_NO[test.list$YEAR == 2020]) &
      !is.na(test.list$INDEX_NO[test.list$YEAR == 2021]) &
      f!= 3
    
    if(length(table(truthfulness[TRUE])) >= 1){
      
      keepers <- test.list
      keepers$COLOR[keepers$YEAR != 2021] <- NA

      p <- ggplot(test.list, aes(x=YEAR, y=INDEX_NO)) +
        geom_ribbon(aes(ymin=CI_LO_NO, ymax=CI_HI_NO), fill=test.list$COLOR, 
                    alpha=0.2, linetype=0) +
        geom_line(lwd=1, col=test.list$COLOR) +
        geom_point(data=keepers,
                   aes(x=YEAR, y=INDEX_NO, col=COLOR)) +
        scale_color_manual(na.value='transparent',
                           values = rev(unique(test.list$COLOR))) +
        labs(x="Year",
             y="Index (No/ tow)",
             title=paste0(survo)) +
        scale_y_continuous(labels=scaleFUN) +
        facet_wrap2(vars(test.list$NAME), scales="free_y", ncol=2, nrow=2)
  
    }
    
    if(nrow(test.list[is.na(test.list$INDEX_KG)==TRUE,]) > 0 &
       nrow(test.list[is.na(test.list$INDEX_KG)==FALSE,]) == 0){
      #plot(p)
      
      ggsave(paste0('C:/Users/klankowicz/Desktop/Index_New_Plots/',
                    survo, " ", stocko,
                    '.png'),
             p)
      rm(p)
      
    }
    
    if(nrow(test.list[is.na(test.list$INDEX_KG)==FALSE,]) > 0){
      p <- p+theme(axis.title.x = element_blank())
      
      if(length(unique(test.list$NAME)) >1){
        q <- ggplot(test.list, aes(x=YEAR, y=INDEX_KG)) + 
          geom_ribbon(aes(ymin=CI_LO_KG, ymax=CI_HI_KG), fill=test.list$COLOR, 
                      alpha=0.2, linetype=0) +
          geom_line(lwd=1, col=test.list$COLOR) +
          labs(x="Year",
               y="Index (Kg/ tow)") +
          #scale_color_manual(values = linecolors) +
          #scale_fill_manual(values=linecolors) +
          scale_y_continuous(labels=scaleFUN) +
          facet_wrap2(vars(test.list$NAME), scales="free_y", ncol=2)+ 
          theme(strip.background = element_blank(), strip.text.x = element_blank())
      }
      
      if(max(test.list$YEAR) == 2021 &
         is.na(test.list$INDEX_NO[test.list$YEAR == 2020]) &
         !is.na(test.list$INDEX_NO[test.list$YEAR == 2021])){
        
          q <- ggplot(test.list, aes(x=YEAR, y=INDEX_KG)) + 
            geom_ribbon(aes(ymin=CI_LO_KG, ymax=CI_HI_KG), fill=test.list$COLOR, 
                        alpha=0.2, linetype=0) +
            geom_line(lwd=1, col=test.list$COLOR) +
            labs(x="Year",
                 y="Index (Kg/ tow)") +
            geom_point(data=test.list[test.list$YEAR == 2021,],
                       aes(x=YEAR, y=INDEX_KG),
                       col=test.list$COLOR[test.list$YEAR == 2021]) +
            #scale_color_manual(values = linecolors) +
            #scale_fill_manual(values=linecolors) +
            scale_y_continuous(labels=scaleFUN) +
            #facet_wrap2(vars(test.list$NAME), scales="free_y", ncol=2)+ 
            theme(strip.background = element_blank(), strip.text.x = element_blank())
      }
      

      
      tmp <- arrangeGrob((p),(q),ncol=1)
      ggsave(paste0('C:/Users/klankowicz/Desktop/Index_New_Plots/',
                    survo, " ", stocko,
                    '.png'),
             tmp)
      
      rm(p,q,tmp)
    }
  }
 rm(survo.list, survo.sub, test.list, ages.tab, survo, stocko) 
}
rm(f, i, r)

#### Age comp plots ####

# Set GGplot auto theme
theme_set(theme(panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.background = element_blank(),
                
                axis.line = element_line(colour = "black"),
                legend.position = "n",
                

                
                plot.title=element_text(size=24, hjust = 0, vjust = 1.2),
                
                strip.text.x = element_text(size = 22)))


for(f in 1:length(unique(numage$SURVEY))){
  survo <- unique(numage$SURVEY)[f]
  survo.sub <- subset(numage, SURVEY == paste0(survo))
  survo.list <- split(survo.sub, f=survo.sub$AREA)

  for(i in 1:length(survo.list)){
    temp <- survo.list[[i]]
    reshp <- dplyr::select(temp, INDEX_NAME, AREA, YEAR, SEASON, Age.0)
    reshp$AGE <- 0
    colnames(reshp) <- c('INDEX_NAME', 'AREA', 'YEAR', 'SEASON', 'VAL', 'AGE')
    
    for(j in 1:8){
      temp2 <- dplyr::select(temp, INDEX_NAME, AREA, YEAR, SEASON, paste0("Age.", j))
      colnames(temp2) <- c('INDEX_NAME', 'AREA', 'YEAR', 'SEASON', 'VAL')
      temp2$AGE <- paste0(j)
      reshp <- rbind(reshp, temp2)
    }
    
    temp3 <- dplyr::select(temp, INDEX_NAME, AREA, YEAR, SEASON, Age.9.)
    colnames(temp3) <- c('INDEX_NAME', 'AREA', 'YEAR', 'SEASON', 'VAL')
    temp3$AGE <- "9+"
    reshp <- rbind(reshp, temp3)
    
    temp4 <- dplyr::select(temp, INDEX_NAME, AREA, YEAR, SEASON, Unk)
    colnames(temp4) <- c('INDEX_NAME', 'AREA', 'YEAR', 'SEASON', 'VAL')
    temp4$AGE <- "UNK."
    reshp <- rbind(reshp, temp4)
    
    reshp <- reshp[with(reshp, order(INDEX_NAME, AREA, SEASON, YEAR, AGE)),]
    row.names(reshp) <- NULL
    
    reshp$AREASEASON <- paste0(reshp$AREA, " ", reshp$SEASON)
    
    stocko <- unique(temp$STOCK)
    
    if(length(unique(reshp$AREASEASON))==2){
      p <- ggplot(reshp, aes(x=AGE, y=YEAR, size = VAL)) +
        geom_point(alpha=0.5) +
        scale_size(range = c(.1, 24)) +
        theme_bw() + 
        theme(axis.text.x=element_text(size=24),
              axis.text.y=element_text(size=24),
              axis.title.x=element_text(size=26),
              axis.title.y=element_text(size=26, angle=90, vjust=2),
              
              legend.position = "n",
              
              plot.title=element_text(size=26, hjust = 0, vjust = 1.2),
              
              strip.text.x = element_text(size = 24)) +
        labs(x="Age (Years)",
             y="Year", 
             title=paste(survo)) +
        facet_wrap2(vars(reshp$AREASEASON))
      
      ggsave(paste0('C:/Users/klankowicz/Desktop/AgeComp_New_Plots/',
                    survo, " ", stocko,
                    '.png'),
             p)
      rm(p)
    }
    
    if(length(unique(reshp$AREASEASON))==1){
      add.tail <- reshp
      add.tail$VAL <- NA
      add.tail$AREASEASON <- 'ZZ'
      
      reshp <- rbind(reshp, add.tail)
      
      reshp$AREASEASON <- factor(reshp$AREASEASON)
      
      p <- ggplot(reshp, aes(x=AGE, y=YEAR, size = VAL)) +
        geom_point(alpha=0.5) +
        scale_size(range = c(.1, 24)) +
        theme_bw() + 
        theme(axis.text.x=element_text(size=24),
              axis.text.y=element_text(size=24),
              axis.title.x=element_text(size=26),
              axis.title.y=element_text(size=26, angle=90, vjust=2),
              
              legend.position = "n",
              
              plot.title=element_text(size=26, hjust = 0, vjust = 1.2),
              
              strip.text.x = element_text(size = 24)) +
        labs(x="  ",
             y="Year", 
             title=paste(survo)) +
        facet_grid(. ~ AREASEASON) +
        annotate("text", x=5, 
                 y=round(min(reshp$YEAR) - (0.1029412 * 
                                              (max(reshp$YEAR) - 
                                                 min(reshp$YEAR))),1), 
                 hjust=0,label="Age (Years)", size=7)+
        coord_cartesian(ylim=c(min(reshp$YEAR), max(reshp$YEAR)),clip="off")
      
      g <- ggplotGrob(p)
      # get the grobs that must be removed
      rm_grobs <- g$layout$name %in% c("panel-1-2", "axis-t-2", "axis-b-2",
                                       "strip-t-2")
      # remove grobs
      g$grobs[rm_grobs] <- NULL
      g$layout <- g$layout[!rm_grobs, ]
      
      ggsave(paste0('C:/Users/klankowicz/Desktop/AgeComp_New_Plots/',
                    survo, " ", stocko,
                    '.png'),
             g)
      rm(p, g)
    }
    
    rm(p, reshp, temp, temp2, temp3, temp4)
  }
}

#### Length Freq plots ####
# Set GGplot auto theme
theme_set(theme(panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.background = element_blank(),
                
                axis.line = element_line(colour = "black"),
                legend.position = "n",
                
                axis.text.x=element_text(size=18, angle=75),
                axis.text.y=element_text(size=18),
                axis.title.x=element_text(size=20),
                axis.title.y=element_text(size=20, angle = 90),
                legend.text=element_text(size=16),
                
                plot.title=element_text(size=20, hjust = 0, vjust = 1.2),
                
                strip.text.x = element_text(size = 16)))

for(f in 1:length(unique(lenfq$SURVEY))){
  survo <- unique(lenfq$SURVEY)[f]
  survo.sub <- subset(lenfq, SURVEY == paste0(survo))
  survo.sub$AREASEASON <- paste0(survo.sub$AREA," ", survo.sub$SEASON)
  ind.split <- split(survo.sub, f=survo.sub$AREASEASON)

  for(i in 1:length(ind.split)){
    stocko <- unique(ind.split[[i]]$AREASEASON)
    nrows <- ceiling(length(table(ind.split[[i]]$YEAR))/4)
    ncols <- 4
    
    if(nrows<=5){
      maxfreq <- plyr::round_any(max(ind.split[[i]]$FREQUENCY), 0.1, ceiling)
      
      yrtab <- unique(ind.split[[i]]$YEAR)[1:20]
      yrtab <- yrtab[is.na(yrtab)==FALSE]
      
      if(length(yrtab)<=4){colnums <- length(yrtab)}
      
      if(length(yrtab)>4){colnums <- 4}
      
      rownums <- ceiling(length(yrtab)/4)
      
      coll1 <- ind.split[[i]][ind.split[[i]]$YEAR %in% yrtab,]
      
      lenfreq1 <- ggplot(data=coll1,
                         aes(x=LENGTH, y=FREQUENCY)) +
        geom_col(fill="black",color="black") +
        labs(x="Length (cm)",
             y="Frequency", 
             title=paste(survo, 
                         coll1$AREA[1],
                         coll1$SEASON[1])) +
        facet_wrap2(vars(coll1$YEAR), ncol = ncols,
                    nrow = nrows) +
        force_panelsizes(rows=unit(rep(1.5, rownums), "in"),
                         col=unit(rep(1.65, colnums), "in"), TRUE) +
        ylim(0,maxfreq)
      
      ggsave(paste0('C:/Users/klankowicz/Desktop/LenFreq_New_Plots/',
                    survo, " ", stocko,
                    '.png'),
             width = 8.5,
             height = 11,
             units = c('in'),
             lenfreq1)
      rm(lenfreq1, coll1)
    }
    
    if(nrows>=6 & nrows<=10){
      maxfreq <- plyr::round_any(max(ind.split[[i]]$FREQUENCY), 0.1, ceiling)
      
      yrtab <- unique(ind.split[[i]]$YEAR)[1:20]
      yrtab <- yrtab[is.na(yrtab)==FALSE]
      
      if(length(yrtab)<=4){colnums <- length(yrtab)}
      
      if(length(yrtab)>4){colnums <- 4}
      
      rownums <- ceiling(length(yrtab)/4)

      coll1 <- ind.split[[i]][ind.split[[i]]$YEAR %in% yrtab,]
      
      lenfreq1 <- ggplot(data=coll1,
                         aes(x=LENGTH, y=FREQUENCY)) +
        geom_col(fill="black",color="black") +
        labs(x="Length (cm)",
             y="Frequency", 
             title=paste(survo, 
                         coll1$AREA[1],
                         coll1$SEASON[1])) +
        facet_wrap2(vars(coll1$YEAR), ncol = ncols,
                    nrow = nrows) + 
        force_panelsizes(rows=unit(rep(1.5, rownums), "in"),
                         col=unit(rep(1.65, colnums), "in"), TRUE) +
        ylim(0,maxfreq)
      
      ggsave(paste0('C:/Users/klankowicz/Desktop/LenFreq_New_Plots/',
                    survo, " ", stocko, " Set1.png"),
             width = 8.5,
             height = 11,
             units = c('in'),
             lenfreq1)
      
      rm(lenfreq1, coll1, yrtab, rownums, colnums)
      
      yrtab2 <- unique(ind.split[[i]]$YEAR)[21:40]
      yrtab2 <- yrtab2[is.na(yrtab2)==FALSE]
      
      if(length(yrtab2)<=4){colnums2 <- length(yrtab2)}
      
      if(length(yrtab2)>4){colnums2 <- 4}
      
      rownums2 <- ceiling(length(yrtab2)/4)
      
      coll2 <- ind.split[[i]][ind.split[[i]]$YEAR %in% yrtab2,]
      
      lenfreq2 <- ggplot(data=coll2,
                         aes(x=LENGTH, y=FREQUENCY)) +
        geom_col(fill="black",color="black") +
        labs(x="Length (cm)",
             y="Frequency", 
             title=paste(survo, 
                         coll2$AREA[1],
                         coll2$SEASON[1])) +
        facet_wrap2(vars(coll2$YEAR), ncol = ncols,
                    nrow = nrows) +
        force_panelsizes(rows=unit(rep(1.5, rownums2), "in"),
                         col=unit(rep(1.65, colnums2), "in"), TRUE) +
        ylim(0,maxfreq)
      
      ggsave(paste0('C:/Users/klankowicz/Desktop/LenFreq_New_Plots/',
                    survo, " ", stocko, " Set2.png"),
             width = 8.5,
             height = 11,
             units = c('in'),
             lenfreq2)
      
      rm(lenfreq2, coll2, yrtab2)

    }
    
    if(nrows>=11 & nrows<=15){
      maxfreq <- plyr::round_any(max(ind.split[[i]]$FREQUENCY), 0.1, ceiling)
      
      yrtab <- unique(ind.split[[i]]$YEAR)[1:20]
      yrtab <- yrtab[is.na(yrtab)==FALSE]
      
      if(length(yrtab)<=4){colnums <- length(yrtab)}
      
      if(length(yrtab)>4){colnums <- 4}
      
      rownums <- ceiling(length(yrtab)/4)
      
      coll1 <- ind.split[[i]][ind.split[[i]]$YEAR %in% yrtab,]
      
      lenfreq1 <- ggplot(data=coll1,
                         aes(x=LENGTH, y=FREQUENCY)) +
        geom_col(fill="black",color="black") +
        labs(x="Length (cm)",
             y="Frequency", 
             title=paste(survo, 
                         coll1$AREA[1],
                         coll1$SEASON[1])) +
        facet_wrap2(vars(coll1$YEAR), ncol = ncols,
                    nrow = nrows) + 
        force_panelsizes(rows=unit(rep(1.5, rownums), "in"),
                         col=unit(rep(1.65, colnums), "in"), TRUE) +
        ylim(0,maxfreq)
      
      ggsave(paste0('C:/Users/klankowicz/Desktop/LenFreq_New_Plots/',
                    survo, " ", stocko, " Set1.png"),
             width = 8.5,
             height = 11,
             units = c('in'),
             lenfreq1)
      
      rm(lenfreq1, coll1, yrtab, rownums, colnums)
      
      yrtab2 <- unique(ind.split[[i]]$YEAR)[21:40]
      yrtab2 <- yrtab2[is.na(yrtab2)==FALSE]
      
      if(length(yrtab2)<=4){colnums2 <- length(yrtab2)}
      
      if(length(yrtab2)>4){colnums2 <- 4}
      
      rownums2 <- ceiling(length(yrtab2)/4)
      
      coll2 <- ind.split[[i]][ind.split[[i]]$YEAR %in% yrtab2,]
      
      lenfreq2 <- ggplot(data=coll2,
                         aes(x=LENGTH, y=FREQUENCY)) +
        geom_col(fill="black",color="black") +
        labs(x="Length (cm)",
             y="Frequency", 
             title=paste(survo, 
                         coll2$AREA[1],
                         coll2$SEASON[1])) +
        facet_wrap2(vars(coll2$YEAR), ncol = ncols,
                    nrow = nrows) +
        force_panelsizes(rows=unit(rep(1.5, rownums2), "in"),
                         col=unit(rep(1.65, colnums2), "in"), TRUE) +
        ylim(0,maxfreq)
      
      ggsave(paste0('C:/Users/klankowicz/Desktop/LenFreq_New_Plots/',
                    survo, " ", stocko, " Set2.png"),
             width = 8.5,
             height = 11,
             units = c('in'),
             lenfreq2)
      
      rm(lenfreq2, coll2, yrtab2)
      
      yrtab3 <- unique(ind.split[[i]]$YEAR)[41:60]
      yrtab3 <- yrtab3[is.na(yrtab3)==FALSE]
      
      if(length(yrtab3)<=4){colnums3 <- length(yrtab3)}
      
      if(length(yrtab3)>4){colnums3 <- 4}
      
      rownums3 <- ceiling(length(yrtab3)/4)
      
      coll3 <- ind.split[[i]][ind.split[[i]]$YEAR %in% yrtab3,]
      
      lenfreq3 <- ggplot(data=coll3,
                         aes(x=LENGTH, y=FREQUENCY)) +
        geom_col(fill="black",color="black") +
        labs(x="Length (cm)",
             y="Frequency", 
             title=paste(survo, 
                         coll3$AREA[1],
                         coll3$SEASON[1])) +
        facet_wrap2(vars(coll3$YEAR), ncol = ncols,
                    nrow = nrows) +
        force_panelsizes(rows=unit(rep(1.5, rownums3), "in"),
                         col=unit(rep(1.65, colnums3), "in"), TRUE) +
        ylim(0,maxfreq)
      
      ggsave(paste0('C:/Users/klankowicz/Desktop/LenFreq_New_Plots/',
                    survo, " ", stocko, " Set3.png"),
             width = 8.5,
             height = 11,
             units = c('in'),
             lenfreq3)
      
      rm(lenfreq3, coll3, yrtab3)
      
    }
  }
}
