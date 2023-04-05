# Function to assign season by GoM cod spawning 
# Fall: Sept - Feb
# Spring: Mar - Aug

spring <- c(3, 4, 5, 6, 7, 8)
fall <- c(9, 10, 11, 12, 1, 2)

spawning_seasons <- function(date){
  if(lubridate::month(date) %in% spring){
    return('SPRING')
  }
  if(lubridate::month(date) %in% fall){
    return('FALL')
  }
}
