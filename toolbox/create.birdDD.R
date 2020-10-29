#' This function converts activity data from MK3005
#' 
#' @param x A number.
#' @param y A number.
#' @return The sum of \code{x} and \code{y}.
#' @examples
#' add(1, 1)
#' add(10, 1)
#' 

create.birdDD <- function(birdGLS, trn){
  
  DDactivity <- create.data.act(trn, birdGLS$activity)
  DDtemperature <- create.data.temp(trn, birdGLS$temperature)
  
  DD <- list(ID = birdGLS$ID, Sex = birdGLS$Sex, days = trn, activity = DDactivity, temperature = DDtemperature)
  class(DD) <- "birdDD"
  return(DD)
}

print.birdDD <- function(birdDD){
  cat("\n-----------------------\n", "birdDD data", "\n----------------------- \n")
  cat("Name : ", birdDD$ID, "\n")
  cat("Lenght of Flight : ", nrow(birdDD$days), " days and nights \n")
  
  cat("\n----- Few Statistics -----")
  cat("\nDays/Nights with high wet activity (>= 75%) : ", sum(birdDD$activity$mean >= 150))
  cat("\nDays/Nights with high dry activity (<= 25%) : ", sum(birdDD$activity$mean <= 50))
  cat("\n")
}

select.birdDD <- function(birdDD, idx){
  birdDD$days <- birdDD$days[idx,]
  birdDD$activity <- birdDD$activity[idx,]
  birdDD$temperature <- birdDD$temperature[idx,]
  return(birdDD)
}
