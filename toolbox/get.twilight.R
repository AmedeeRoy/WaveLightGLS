#------------------------------------------------------------------------------------------------------
# GET.TWILIGHT
#------------------------------------------------------------------------------------------------------
# - adapt function from TWGEOS for twilight estimation

get.twilight <- function(light, THRESHOLD, DARK.MIN){
  # twgeos
  resolution = as.double(light$dt[2]-light$dt[1])
  tagdata = light[, c("dt","value")]
  colnames(tagdata)=c("Date", "Light")
  tw <- findTwilights(tagdata, threshold = THRESHOLD, include= tagdata$Date, 
                      dark.min = DARK.MIN)
  tw <- twilightAdjust(tw, resolution*60)
  return(tw)
}

