rm(list=ls())

sma <- read.csv(here('VAST_runs/small/Overall_BC/Index.csv'))
load(here('VAST_runs/small/Overall_BC/Overall_BC_smallcod_allstrat_natsplin_fsON.RData'))
sma.strat <- fit$settings$strata.limits
rm(list=setdiff(ls(), c('sma', 'sma.strat')))

med <- read.csv(here('VAST_runs/medium/Overall_BC/WGOM/Index.csv'))
load(here('VAST_runs/medium/Overall_BC/ALL/Overall_BC_medcod_allstrat_natsplin_fsON_ALL.RData'))
med.strat <- fit$settings$strata.limits
rm(list=setdiff(ls(), c('sma', 'sma.strat', 'med', 'med.strat')))

lar <- read.csv(here('VAST_runs/large/Overall_noBC/Index.csv'))
load(here('VAST_runs/large/Overall_noBC/Overall_noBC_larcod_allstrat_natsplin_fsON.RData'))
lar.strat <- fit$settings$strata.limits
rm(list=setdiff(ls(), c('sma', 'sma.strat', 'med', 'med.strat', 'lar', 'lar.strat')))

sma.strat$Stratum <- unique(sma$Stratum)
sma <- left_join(sma, sma.strat, by=c('Stratum'))

med.strat$Stratum <- unique(med$Stratum)
med <- left_join(med, med.strat, by=c('Stratum'))

lar.strat$Stratum <- unique(lar$Stratum)
lar <- left_join(lar, lar.strat, by=c('Stratum'))

full_data <- rbind(
  cbind(sma, Size='Small'),
  cbind(med, Size='Medium'),
  cbind(lar, Size='Large')
)

full_data <- full_data %>% 
  separate(Time, into=c('Year', 'Season')) %>% 
  mutate(Year = as.numeric(Year),
         Season = factor(Season, levels=c('Spring', 'Fall'))) %>% 
  dplyr::select(-Category, -Stratum, -Units) %>% 
  mutate(Size = factor(Size, levels=c('Small', 'Medium', 'Large'))) %>% 
  rename(sterr.est = Std..Error.for.Estimate,
         sterr.ln = Std..Error.for.ln.Estimate.,
         Stratum=STRATA)

allind <- ggplot() + 
  geom_line(data=full_data[full_data$Stratum != 'ALL',],
            aes(x=Year, y=Estimate, col=Stratum)) +
  geom_ribbon(data=full_data[full_data$Stratum != 'ALL',],
            aes(x=Year, ymin=Estimate-sterr.est, 
                ymax=Estimate+sterr.est,
                fill=Stratum),
            alpha=0.3) +
  facet_grid(vars(Size),vars(Season),  scales="free")  

ggsave(plot=allind,
       filename=here("VAST_Runs/All_Indices.png"),
       width=7.5, height=6, units='in')
