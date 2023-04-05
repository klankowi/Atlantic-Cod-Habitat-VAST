# Function to assign season by GoM cod spawning 
# Fall: Sept - Feb
# Spring: Mar - Aug

winter <- c(1, 2)
spring <- c(3, 4, 5, 6)
summer <- c(7, 8)
fall <- c(9, 10, 11, 12)

true_seasons <- function(date){
  if(lubridate::month(date) %in% winter){
    return('WINTER')
  }
  if(lubridate::month(date) %in% spring){
    return('SPRING')
  }
  if(lubridate::month(date) %in% summer){
    return('SUMMER')
  }
  if(lubridate::month(date) %in% fall){
    return('FALL')
  }
}
