#------------------------------------------------------------------------------------------------------
# CREATE.DATA.TEMP
#------------------------------------------------------------------------------------------------------
# - compute mean, var, min, max, nb.dive of temperature between each couple of twilight

create.data.temp <- function(trn, temperature){
  
  data.temp = matrix(ncol = 5, nrow = length(trn$tFirst))
  
  for (i in 1:length(trn$tFirst)){
    
    data = temperature[which(temperature$dt >= trn$tFirst[i] & 
                               temperature$dt <= trn$tSecond[i]),]$value
    
    if (length(data != 0)){
      mean = mean(data)
      var = var(data)
      min = min(data)
      max = max(data)
      nb = length(data)
      data.temp[i,] = c(mean, var, min, max, nb)
    }
  }
  
  data.temp <- data.frame(data.temp)
  colnames(data.temp)<-c("mean", "var", "min", "max", "nb")
  return(data.temp)
}