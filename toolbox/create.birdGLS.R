create.birdGLS <- function(lig, act, tem, metadata, name){
  
  # Read CSV into R
  light <- read.csv(lig, header=TRUE, sep=",")
  activity <- read.csv(act, header=TRUE, sep=",")
  temperature <- read.csv(tem, header=TRUE, sep=",")
  
  for(k in colnames(metadata)){
    assign(k, metadata[metadata$ID == name, k])
  }
  Time.Deployment <- as.POSIXct(Time.Deployment, format = "%d/%m/%Y %H:%M", tz="GMT")
  Time.Recovery <- as.POSIXct(Time.Recovery, format = "%d/%m/%Y %H:%M", tz="GMT")
  
  # Formatting Light
  light <- light[, c(1,2,4)]
  colnames(light) <- c("quality","dt", "value")
  light$dt <-as.POSIXct(strptime(as.character(light$dt),"%d/%m/%y %H:%M:%S"), "GMT")
  light$quality <- as.character(light$quality)
  light <- light[which(light$dt >= Time.Deployment & light$dt <= Time.Recovery),]
  
  # Formatting Temperature
  temperature[,4] <- temperature[, 4] + temperature[,5]/10^6
  temperature <- temperature[, c(1,2,4)]
  colnames(temperature) <- c("quality","dt", "value")
  temperature$dt <-as.POSIXct(strptime(as.character(temperature$dt),"%d/%m/%y %H:%M:%S"), "GMT")
  temperature$quality <- as.character(temperature$quality)
  temperature <- temperature[which(temperature$dt >= Time.Deployment & temperature$dt <= Time.Recovery),]
  
  # Formatting Activity
  if (Model=="MK3005"){
    # For MK3005 sensors : modification of activity data
    colnames(activity) <- c("Quality", "DateTime", "Day01011900", "Duration", "WD")
    activity$DateTime <- as.POSIXct(strptime(as.character(activity$DateTime),"%d/%m/%y %H:%M:%S"), "GMT") 
    activity <- activity.converter(activity)
    activity <- activity[which(activity$dt >= Time.Deployment & activity$dt <= Time.Recovery),]
    
  } else {
    # For MK3006 sensors : no modification
    colnames(activity) <- c("quality","dt","a", "value","temperature")
    activity$a <- NULL
    activity <- activity[!is.na(activity$value),]
    activity$dt <- as.POSIXct(strptime(as.character(activity$dt),"%d/%m/%y %H:%M:%S"), "GMT") 
  }
  activity$quality = as.character(activity$quality)
  
  birdGLS <- list(light = light, activity = activity, temperature = temperature, ID = as.character(ID),
                  Species = as.character(Species),
                  Model = as.character(Model),
                  Time.Deployment = Time.Deployment,
                  Time.Recovery = Time.Recovery, Pos.Deployment = c(Long.Deployment, Lat.Deployment),
                  Pos.Recovery = c(Long.Recovery, Lat.Recovery), Sex = as.character(Sex))
  
  class(birdGLS) <- "birdGLS"
  return(birdGLS)
}

print.birdGLS <- function(birdGLS){
  cat("\n-----------------------\n", "birdGLS data", "\n----------------------- \n")
  cat("Name : ", birdGLS$ID, "\n")
  cat("Lenght of Flight : ", difftime(birdGLS$Time.Recovery, birdGLS$Time.Deployment, units = "days"), " days \n")
  
  cat("\n---\nQuality Light Data : \n")
  print(table(birdGLS$light$quality))
  cat("\n---\nQuality Activity Data : \n")
  print(table(birdGLS$activity$quality))
  cat("\n---\nQuality Temperature Data : \n")
  print(table(birdGLS$temperature$quality))
}

select.birdGLS <- function(birdGLS, Time.Deployment, Time.Recovery, Pos.Deployment = NULL, Pos.Recovery = NULL){
  
  birdGLS$Time.Deployment = Time.Deployment
  birdGLS$Time.Recovery = Time.Recovery
  birdGLS$Pos.Deployment = Pos.Deployment
  birdGLS$Pos.Recovery = Pos.Recovery
  
  birdGLS$light <- birdGLS$light[which(birdGLS$light$dt >= Time.Deployment & birdGLS$light$dt <= Time.Recovery),] 
  birdGLS$activity <- birdGLS$activity[which(birdGLS$activity$dt >= Time.Deployment & birdGLS$activity$dt <= Time.Recovery),]
  birdGLS$temperature <- birdGLS$temperature[which(birdGLS$temperature$dt >= Time.Deployment & birdGLS$temperature$dt <= Time.Recovery),]
  
  return(birdGLS)
}

plot.birdGLS <- function(birdGLS){
  layout(matrix(c(1,2,3), nrow = 3))
  
  plot(birdGLS$light$dt, birdGLS$light$value, pch = 19, xlab = "time", ylab = "light", 
       col = ifelse(birdGLS$light$quality == "ok", 0, 2))
  lines(birdGLS$light$dt, birdGLS$light$value)
  
  plot(birdGLS$activity$dt, birdGLS$activity$value, pch = 19, xlab = "time", ylab = "activity", 
       col = ifelse(birdGLS$activity$quality == "ok", 0, 2))
  lines(birdGLS$activity$dt, birdGLS$activity$value)
  
  plot(birdGLS$temperature$dt, birdGLS$temperature$value, pch= 19, xlab = "time", ylab = "temp", 
       col = ifelse(birdGLS$temperature$quality == "ok", 0, 2))
  lines(birdGLS$temperature$dt, birdGLS$temperature$value)
}
