plot.kde.coord <- function(coord, alpha = 0.2, H=NULL, N=NULL, eez = EEZ_br, map = world, title = "", col){
  
  coord_noNA = coord[!is.na(rowSums(coord)),]
  
  fdn <- data.frame(longitude = c(-32.43, -33.80, -29.35), latitude = c(-3.85, -3.86, 0.92), 
                    z=c("Fernando de Noronha","Atol das Rocas","Sao Pedro Sao Paolo"))
  
  country.label <- data.frame(longitude = c(-50, -47, -37.25), latitude = c(-12, 9, -15), 
                              name = c("Brazil", "Atlantic Ocean", "EEZ"))
  
  CRDS <- data.frame(x = coord_noNA[,1], y = coord_noNA[,2])

  f1 <- with(CRDS, kde2d(CRDS[,1], CRDS[,2], n = N, h = H, lims = c(-60, 0, -30, 20)))
  LKL <- data.frame(expand.grid(x = f1$x, y = f1$y), z = as.vector(f1$z))
  LKL[,3] <- LKL[,3]/sum(LKL[,3])
  colnames(LKL) <- c("longitude", "latitude", "norm")
  
  
  ### COMPUTE PRESENCE IN ZEE
  p = Polygon(cbind(eez$long, eez$lat))
  ps = Polygons(list(p),1)
  sps = SpatialPolygons(list(ps))
  PRESENCE = sum(LKL[!is.na(over(SpatialPoints(LKL[,1:2]), sps)), 3])
    
  
  MAX <- max(LKL$norm, na.rm = TRUE)
  MIN <- min(LKL$norm, na.rm = TRUE)
  
  prob.pre <- function(x){
    return(abs(sum(LKL$norm[LKL$norm>= x], na.rm = TRUE) - quantile))
  }
  bk <- NULL
  for (quantile in c(0.1, 0.25, 0.5, 0.75, 0.9)){
    max <-  optimize(f = prob.pre, interval = c(1e-100,MAX))
    bk <- c(bk, max$minimum)
  }
  
  gp <- ggplot() + ggtitle(title, subtitle = paste0("Presence in EEZ : " , round(PRESENCE*100), "%")) +
    geom_polygon(data = eez[-(250:1550),], mapping = aes(long,lat, alpha = 1), fill = "gray",
                 colour = "black", linetype = 2) +
    theme_bw() + theme(legend.position = "none", plot.margin = unit(c(.2,.2,.2,.2), units = "lines"))+  xlab("") + ylab("") +
    geom_point(data = CRDS, mapping= aes(x, y),
               size = 1, col = col, alpha= alpha) +
    scale_color_gradientn(colours = c("white", col,
                                      adjustcolor(col,red.f=0.75, blue.f=0.75, green.f=0.75),
                                      adjustcolor(col,red.f=0.5, blue.f=0.5, green.f=0.5)),
                            aesthetics = c("colour", "fill"), na.value = "white") +
    geom_tile(aes(x = longitude, y = latitude, fill = norm, alpha= 1), data = LKL)+
    geom_contour(aes(x = longitude, y = latitude, z = norm), data = LKL, 
                 colour = adjustcolor(col,red.f=0.75, blue.f=0.75, green.f=0.75),
                 na.rm=TRUE, breaks=bk) +
    geom_map(data = map, map=map, aes(map_id=id),
           fill="darkgray", color="#7f7f7f", size=0.5) +
    coord_fixed(xlim = c(-55, -25), ylim = c(-20, 10)) +
    geom_text(data=country.label, mapping = aes(longitude, latitude, label = name), fontface="italic") +
    geom_point(data = fdn, mapping= aes(longitude, latitude, shape= factor(z)),
               size = 2.5)
  return(gp)
}

plot.map <- function(m, eez = EEZ_br, map = world, title = title){
  

  fdn <- data.frame(longitude = c(-32.39, -33.80, -29.35), latitude = c(-3.8, -3.86, 0.92), 
                    z=c("Fernando de Noronha","Atol das Rocas","Sao Pedro Sao Paolo"))
  
  country.label <- data.frame(longitude = c(-50, -47, -37.25), latitude = c(-12, 9, -15), 
                              name = c("Brazil", "Atlantic Ocean", "EEZ"))
  
  
  LKL <- data.frame(expand.grid(x = m$x, y = m$y), z = as.vector(m$z))
  LKL[,3] <- LKL[,3]/sum(LKL[,3])
  colnames(LKL) <- c("longitude", "latitude", "norm")
  
  ### COMPUTE PRESENCE IN ZEE
  p = Polygon(cbind(EEZ_br$long, EEZ_br$lat))
  ps = Polygons(list(p),1)
  sps = SpatialPolygons(list(ps))
  PRESENCE = sum(LKL[!is.na(over(SpatialPoints(LKL[,1:2]), sps)), 3])
  
  MAX <- max(LKL$norm, na.rm = TRUE)
  MIN <- min(LKL$norm, na.rm = TRUE)
  
  prob.pre <- function(x){
    return(abs(sum(LKL$norm[LKL$norm>= x], na.rm = TRUE) - quantile))
  }
  bk <- NULL
  for (quantile in c(0.25, 0.5, 0.75, 0.95)){
    max <-  optimize(f = prob.pre, interval = c(1e-100,MAX))
    bk <- c(bk, max$minimum)
  }
  
  gp <- ggplot() + ggtitle(title, subtitle = paste0("Presence in ZEE : " , round(PRESENCE, 2))) +
    geom_polygon(data = EEZ_br[-(250:1550),], mapping = aes(long,lat, alpha = 1), fill = "gray",
                 colour = "black", linetype = 2) +
    scale_fill_gradientn("Probability of presence", colours= c("white", "yellow", "red"), 
                         na.value = "white", values=c(0, 0.015, 1)) +
    geom_tile(aes(x = longitude, y = latitude, fill = norm, alpha= 1), data = LKL)+
    geom_contour(aes(x = longitude, y = latitude, z = norm), data = LKL, colour = "black", na.rm=TRUE,
                 breaks=bk) +
    geom_map(data=map, map=map, aes(map_id=id),
             fill="darkgray", color="#7f7f7f", size=0.5) + 
    coord_fixed(xlim = c(-55, -25), ylim = c(-20, 10)) + 
    # geom_text(data=country.label, mapping = aes(longitude, latitude, label = name), fontface="italic") +
    # geom_point(data = crds, mapping= aes(x, y),
    #            size = 1, col = "red", alpha= 0.05) +
    geom_point(data = fdn, mapping= aes(longitude, latitude, shape= factor(z)), 
               size = 3) + 
    theme_bw() + theme(legend.position = "none")+  xlab("") + ylab("") 
  return(gp)
}
