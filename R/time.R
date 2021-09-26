get.sundial <- function(object,twilight="civil")
{
  twilight <- match.arg(twilight,c('nautical','civil','none'))
  if(twilight=='civil')
  { keep <- c("dawn","solarNoon","dusk","nadir") }
  else if(twilight=='nautical')
  { keep <- c("nauticalDawn","solarNoon","nauticalDusk","nadir") }
  else if(twilight=='none')
  { keep <- c("sunrise","solarNoon","sunset","nadir") }

  TODAY <- data.frame(date=object$timestamp,lat=object$latitude,lon=object$longitude)
  # calculate local timezones
  tz <- round(object$lon/15)
  # specify date in local timezone
  TODAY$date <- as.Date(TODAY$date,tz=tz)

  YESTERDAY <- TODAY
  YESTERDAY$date <- YESTERDAY$date - 1

  TOMORROW <- TODAY
  TOMORROW$date <- TOMORROW$date + 1

  YESTERDAY <- suncalc::getSunlightTimes(data=YESTERDAY,keep=keep,tz="UTC")
  TODAY <- suncalc::getSunlightTimes(data=TODAY,keep=keep,tz="UTC")
  TOMORROW <- suncalc::getSunlightTimes(data=TOMORROW,keep=keep,tz="UTC")

  YESTERDAY <- as.numeric(YESTERDAY)
  TODAY <- as.numeric(TODAY)
  TOMORROW <- as.numeric(TOMORROW)



}
