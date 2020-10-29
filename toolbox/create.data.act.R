#------------------------------------------------------------------------------------------------------
# CREATE.DATA.ACT
#------------------------------------------------------------------------------------------------------
# - compute mean, var, nb.dry, nb.wet, nb.dive of activity between each couple of twilight

create.data.act <- function(trn, activity){
  
  data.activity = matrix(ncol = 5, nrow = length(trn$tFirst))
  
  for (i in 1:length(trn$tFirst)){
    
    data = activity[which(activity$dt >= trn$tFirst[i] & 
                            activity$dt <= trn$tSecond[i]),]$value
    
    
    mean = mean(data)
    var = var(data)
    nb.dry = sum(data == 0)
    nb.wet = sum(data == 200)
    nb.dive = sum(data > 0 & data < 200)
    
    data.activity[i,] = c(mean, var, nb.dry, nb.wet, nb.dive)
  }
  
  data.activity <- data.frame(data.activity)
  colnames(data.activity)<-c("mean", "var", "nb.dry", "nb.wet" , "nb.dive")
  return(data.activity)
}
