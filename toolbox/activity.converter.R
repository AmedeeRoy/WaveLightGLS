activity.converter <- function(a){
  # Modifying activity structure
  if (length(a$Duration)>1){
    # remove error
    for (kk in which(is.na(a$Duration))){
      a$Duration[kk-1] <- difftime(a$DateTime[kk+1], a$DateTime[kk-1], units= "secs")
    }
    a <- a[-which(is.na(a$Duration)),]
    n <- nrow(a)
    a$WDn <- as.numeric(a$WD)-2  # 0:dry ; 1:wet
    
    # if problem remove suspect data
    if (min(a$Duration) <=0){
      a<- a[a$Quality=="ok",]
    }
  } else {
    a$WDn <- as.numeric(a$WD)-1  # 0:dry
  }
  
  a$Qn <- ifelse(a$Quality == "ok", 1, 0)
  
  WD <- rep(a$WDn, times=a$Duration)  # wet/dry, 1 point par sec
  Q <- rep(a$Qn, times=a$Duration)    # Quality, 1 point par sec
  
  dft <- data.frame(DateTime=1:length(WD)-1, WD, Q) # la base temps (secs depuis le debut)
  lim <- seq(0, length(dft$WD)+60*10, by=60*10)     # base temps 10 min pour decoupage
  cutdft <- cut(dft$DateTime, breaks=lim, include.lowest = TRUE)  # decoupage, tourne en 15 secs
  sdft <- as.numeric(tapply(dft$WD, cutdft, sum))/3  # calcul de la somme (= nombre de wet par 10 min, /3 pour max=200), tourne en 6 secs
  qdft <- ifelse(as.numeric(tapply(dft$Q, cutdft, sum)) >= 590, "ok", "SUSPECT")
  a10 <- data.frame(quality = qdft, dt = lim[1:(length(lim)-1)] + a$DateTime[1], value = sdft)  # assemblage. DateTime = debut du segment de 10 min

  return(a10)
}