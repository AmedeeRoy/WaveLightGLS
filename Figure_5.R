# WaveLightGLS
#
# Figure 5
#---------------------------------

list.data = list.files("./results", pattern=".RData", all.files=FALSE,full.names=TRUE)
list.data <- list.data[grepl("SUDA", list.data)]

twilight.dev <- list()

## Calibration function
calib <- function(twl_c, lat, zenith.start = 96) {
  z0    <- seq(zenith.start-10, zenith.start+10, by = 0.2)
  crds1 <- lapply(cbind(z0), function(x) thresholdPath(twl_c$Twilight, twl_c$Rise, zenith = x)$x)
  dist1 <- unlist(lapply(crds1, function(x1) median(abs(x1[,2]-lat))))
  z0[which.min(dist1)]
}

###################################
## Calculation of Twiligth Error ##
###################################
lon.breed <- -32.4255
lat.breed <- -3.8496

tm <- seq(as.POSIXct("2017-05-04", tz = "GMT"), as.POSIXct("2018-04-23", tz = "GMT"), by = "day")
rise <- rep(c(TRUE, FALSE), length(tm))
c.dat <- data.frame(Twilight = twilight(rep(tm, each = 2), lon = lon.breed, lat = lat.breed, 
                                        rise = rise, zenith = 93), Rise = rise)

calib.tm  <- c(as.POSIXct("2017-05-10", tz = "GMT"), as.POSIXct("2017-06-15", tz = "GMT"))  

CALIBRATION <- NULL
DATA <- NULL
i = 1
for (data in list.data){
  load(data)
  
  ### CALIBRATION
  twl <- geolight.convert(birdDD$days$tFirst, birdDD$days$tSecond, birdDD$days$type)
  twl_calib <- subset(twl, Twilight>=calib.tm[1] & Twilight<=calib.tm[2])
  sun  <- solar(twl_calib[,1])
  z    <- refracted(zenith(sun, lon.breed, lat.breed))
  
  twl_t   <- twilight(twl_calib[,1], lon.breed, lat.breed, rise = twl_calib[,2], zenith = max(z)+0.1)
  twl_dev <- ifelse(twl_calib$Rise, as.numeric(difftime(twl_calib[,1], twl_t, units = "mins")),
                    as.numeric(difftime(twl_t, twl_calib[,1], units = "mins")))
  
  png(paste0('./calibration/', birdGLS$ID, '.png'))
  hist(twl_dev, main = birdGLS$ID, freq = F, breaks = 26)
  seq <- seq(0, 80, length = 100)
  fitml_ng <- fitdistr(twl_dev, "gamma")
  lines(seq, dgamma(seq, fitml_ng$estimate[1], fitml_ng$estimate[2]), col = "firebrick", lwd = 3, lty = 2)
  dev.off()
  
  
  out <- data.frame(bird = birdGLS$ID, zenith.median = median(z), zenith.max = max(z),
                    shape = fitml_ng$estimate[1], scale = fitml_ng$estimate[2],
                    model = birdGLS$Model, twl_dev = twl_dev)
  CALIBRATION <- rbind(CALIBRATION, out)
  
  
  ### ALL DEPLOYMENT
  twl_dev_all0 <- twilight(twl[,1], lon.breed, lat.breed, rise = twl[,2], zenith = max(z)+0.1)
  twl_dev_all  <- ifelse(twl$Rise, as.numeric(difftime(twl[,1], twl_dev_all0, units = "mins")),
                         as.numeric(difftime(twl_dev_all0, twl[,1], units = "mins")))
  
  zenith <- calib(twl_calib, lat.breed, 96)
  crds <- thresholdPath(twl$Twilight, twl$Rise, zenith = zenith)
  
  if(nrow(birdDD$activity)+1 == nrow(crds$x)){
    act <- c(NA,birdDD$activity$mean)
  } else{
    act <- rep(NA, nrow(crds$x))
  }
  
  out <- data.frame(bird = birdGLS$ID, time = crds$time, zenithT = median(z), zenith = max(z), 
                    tw_error = twl_dev_all, lon = crds$x[,1], lat =  crds$x[,2],
                    act = act)
  DATA <- rbind(DATA, out)
  
  cat(i, ' out of ', length(list.data), '\n')
  i = i+1
}

### TEMPERATURE ERROR AND DEVIATIONS
DATA$temp_fdn_sat <-  getSSTPoint(path = "./data/METOFFICE-GLO-SST-L4-REP-OBS-SST_1590497114049.nc",
                                  coord = matrix(rep(birdGLS$Pos.Deployment, nrow(DATA)), byrow = TRUE, ncol = 2),
                                  time = DATA$time) -273.15

DATA$temp_sat <-  getSSTPoint(path = "./data/METOFFICE-GLO-SST-L4-REP-OBS-SST_1590497114049.nc",
                              coord = DATA[,c("lon", "lat")],
                              time = DATA$time) -273.15

DATA$diff <- DATA$temp_sat - DATA$temp_fdn_sat


# Read Metadata and select relevant data
options(stringsAsFactors = FALSE)
metadata = "./data/Metadata_GLS.csv"
metadata <- read.csv(metadata, header=TRUE, sep=",")
# selection
DATA$sex <- sapply(DATA$bird, function(b){metadata$Sex[metadata$ID == b]})


START = "2017-05-05"
END = "2018-04-23"
days <- seq(as.Date(START), as.Date(END), by ="day")
TIME <- data.frame( date = rep(days, 2),
                    type = c(rep(1, length(days)), rep(2, length(days))))
TIME <- TIME[order(TIME$date, TIME$type),]

### Utilization Distribution MALE/FEMALE
# load graphic data
data(wrld_simpl)
wrld_simpl@data$id <- wrld_simpl@data$NAME
world <- fortify(wrld_simpl)

eez<-readOGR("./data/World_EEZ.shp", "World_EEZ")
EEZ <-  fortify(eez)
EEZ_br <- EEZ[which(EEZ$id==163 & !EEZ$hole),]

eez = EEZ_br
map = world
fdn <- data.frame(longitude = c(-32.43, -33.80, -29.35), latitude = c(-3.85, -3.86, 0.92), 
                  z=c("Fernando de Noronha","Atol das Rocas","Sao Pedro Sao Paolo"))

country.label <- data.frame(longitude = c(-50, -47, -37.25), latitude = c(-12, 9, -15), 
                            name = c("Brazil", "Atlantic Ocean", "EEZ"))

p = Polygon(cbind(eez$long, eez$lat))
ps = Polygons(list(p),1)
sps = SpatialPolygons(list(ps))

get_lkl <- function(coord){
  coord_noNA = coord[!is.na(rowSums(coord)),]
  CRDS <- data.frame(x = coord_noNA[,1], y = coord_noNA[,2])
  f1 <- with(CRDS, kde2d(CRDS[,1], CRDS[,2], n = 100, h = 2, lims = c(-60, 0, -30, 20)))
  LKL <- data.frame(expand.grid(x = f1$x, y = f1$y), z = as.vector(f1$z))
  LKL[,3] <- LKL[,3]/sum(LKL[,3])
  colnames(LKL) <- c("longitude", "latitude", "norm")
  
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
  
  center = data.frame(lon = mean(coord_noNA$lon), lat = mean(coord_noNA$lat))
  return(list(grid = LKL, bk = bk, center = center))
}


library(ggplotify)
library(gridExtra)
library(grid)


## PERIOD (a)
start = TIME$date[110]
end = TIME$date[200]
dd <- DATA[which(DATA$time>= start & DATA$time<= end),]
act.M <- dd$act[dd$sex == 'M']/2
act.F <- dd$act[dd$sex == 'F']/2
hm <- table(round(act.M, digits = -1))
hf <- table(round(act.F, digits = -1))
a <- data.frame(M = c(hm/sum(hm)), F = c(hf/sum(hf)))
p1_a <- as.ggplot(~barplot(as.matrix(t(a)), col = c('firebrick', '#9ECAE1'), beside = TRUE, las = 1,
                         xlab = '', ylab = '', ylim = c(0,0.6), 
                         cex.lab = 1) +
                  mtext(text = 'Proportion of time in water (%)', 1, 2.5) + 
                  mtext(text = 'Frequency', 2, 2.5),
                scale = 1)
dev.M <- dd$lon[which(dd$sex == 'M' & dd$act >150)] -lon.breed
dev.F <- dd$lon[which(dd$sex == 'F' & dd$act >150)] -lon.breed
d.M <- density(dev.M , n = 1024)
d.F <- density(dev.F , n = 1024)
test <- t.test(dev.M, dev.F)
xx.M <-  (d.M$x > -1) & (d.M$x < 4) 
xx.F <-  (d.F$x > -1) & (d.F$x < 4) 
p2_a <- as.ggplot(~plot(d.M$x[xx.M], d.M$y[xx.M], xlim = c(-1, 4), ylim = c(0, 0.7), type = 'l', col = 'firebrick',
                        las = 1, lwd = 3, xlab = "Eastward Deviation from FdN (°) \n when mean wet time >75%", ylab = "Density", cex.lab = 1, mgp=c(3,1,0))+
                    lines(d.F$x[xx.F], d.F$y[xx.F], lwd = 3, col = '#9ECAE1') + 
                    abline(v = 0, lty = 2) +
                    text(2.8, 0.67, paste('p-value :', round(test$p.value, 2)), cex = 1),
                  scale = 0.95)

select <- (DATA$time >= start) & (DATA$time <= end)  & 
  (abs(DATA$diff) <0.5)

data_m <- DATA[ select & (DATA$sex == 'M') ,c("lon", "lat", "act")]
data_f <- DATA[ select & (DATA$sex == 'F') ,c("lon", "lat", "act")]

lkl_m <- get_lkl(data_m)
lkl_f <- get_lkl(data_f)

map_a <- ggplot() +
  geom_polygon(data = eez[-(250:1550),], mapping = aes(long,lat, alpha = 1), fill = "gray",
               colour = "white", linetype = 1, show.legend = FALSE) +
  theme_bw() + theme(legend.position = c(0.78, 0.89), 
                     legend.text = element_text(size=10),
                     plot.margin = unit(c(.2,.2,.2,.2), units = "lines"))+  xlab("") + ylab("") +
  geom_contour(aes(x = longitude, y = latitude, z = norm), data = lkl_m$grid,
               size = 1.5,
               colour = adjustcolor('firebrick',red.f=0.75, blue.f=0.75, green.f=0.75),
               na.rm=TRUE, breaks=lkl_m$bk[5], show.legend = FALSE) +
  geom_contour(aes(x = longitude, y = latitude, z = norm), data = lkl_f$grid, 
               size = 1.5,
               colour = adjustcolor('#9ECAE1',red.f=0.75, blue.f=0.75, green.f=0.75),
               na.rm=TRUE, breaks=lkl_m$bk[5], show.legend = FALSE) +
  geom_map(data = map, map=map, aes(map_id=id),
           fill="darkgray", color="#7f7f7f", size=0.5, show.legend = FALSE) +
  coord_fixed(xlim = c(-55, -25), ylim = c(-20, 10)) + 
  geom_point(data = lkl_m$center, mapping= aes(lon, lat), size = 3, stroke = 2, shape = 3, col = 'firebrick', show.legend = FALSE)+
  geom_point(data = lkl_f$center, mapping= aes(lon, lat), size = 3, stroke = 2, shape = 3, col = '#9ECAE1', show.legend = FALSE)+
  geom_point(data = fdn, mapping= aes(longitude, latitude, shape= factor(z), size = 1.5), show.legend = FALSE) +
  geom_text(data=country.label, mapping = aes(longitude, latitude, label = name), fontface="italic")+
  ggtitle("(a) Early chick-rearing")


## PERIOD (b)
start = TIME$date[230]
end = TIME$date[450]
dd <- DATA[which(DATA$time>= start & DATA$time<= end),]
act.M <- dd$act[dd$sex == 'M']/2
act.F <- dd$act[dd$sex == 'F']/2
hm <- table(round(act.M, digits = -1))
hf <- table(round(act.F, digits = -1))
a <- data.frame(M = c(hm/sum(hm)), F = c(hf/sum(hf)))
p1_b <- as.ggplot(~barplot(as.matrix(t(a)), col = c('firebrick', '#9ECAE1'), beside = TRUE, las = 1,
                           xlab = '', ylab = '', ylim = c(0,0.6), 
                           cex.lab = 1) +
                    mtext(text = 'Proportion of time in water (%)', 1, 2.5) + 
                    mtext(text = 'Frequency', 2, 2.5),
                  scale = 1)
dev.M <- dd$lon[which(dd$sex == 'M' & dd$act >150)] -lon.breed
dev.F <- dd$lon[which(dd$sex == 'F' & dd$act >150)] -lon.breed
d.M <- density(dev.M , n = 1024)
d.F <- density(dev.F , n = 1024)
test <- t.test(dev.M, dev.F)
xx.M <-  (d.M$x > -1) & (d.M$x < 4) 
xx.F <-  (d.F$x > -1) & (d.F$x < 4) 
p2_b <- as.ggplot(~plot(d.M$x[xx.M], d.M$y[xx.M], xlim = c(-1, 4), ylim = c(0, 0.7), type = 'l', col = 'firebrick',
                        las = 1, lwd = 3, xlab = "Eastward Deviation from FdN (°) \n when mean wet time >75%", ylab = "Density", cex.lab = 1, mgp=c(3,1,0))+
                    lines(d.F$x[xx.F], d.F$y[xx.F], lwd = 3, col = '#9ECAE1') + 
                    abline(v = 0, lty = 2) +
                    text(2.8, 0.67, paste('p-value :', round(test$p.value, 2)), cex = 1),
                  scale = 0.95)

select <- (DATA$time >= start) & (DATA$time <= end)  & 
  (abs(DATA$diff) <0.5)

data_m <- DATA[ select & (DATA$sex == 'M') ,c("lon", "lat", "act")]
data_f <- DATA[ select & (DATA$sex == 'F') ,c("lon", "lat", "act")]

lkl_m <- get_lkl(data_m)
lkl_f <- get_lkl(data_f)

map_b <- ggplot() +
  geom_polygon(data = eez[-(250:1550),], mapping = aes(long,lat, alpha = 1), fill = "gray",
               colour = "white", linetype = 1, show.legend = FALSE) +
  theme_bw() + theme(legend.position = c(0.78, 0.89), 
                     legend.text = element_text(size=10),
                     plot.margin = unit(c(.2,.2,.2,.2), units = "lines"))+  xlab("") + ylab("") +
  geom_contour(aes(x = longitude, y = latitude, z = norm), data = lkl_m$grid,
               size = 1.5,
               colour = adjustcolor('firebrick',red.f=0.75, blue.f=0.75, green.f=0.75),
               na.rm=TRUE, breaks=lkl_m$bk[5], show.legend = FALSE) +
  geom_contour(aes(x = longitude, y = latitude, z = norm), data = lkl_f$grid, 
               size = 1.5,
               colour = adjustcolor('#9ECAE1',red.f=0.75, blue.f=0.75, green.f=0.75),
               na.rm=TRUE, breaks=lkl_m$bk[5], show.legend = FALSE) +
  geom_map(data = map, map=map, aes(map_id=id),
           fill="darkgray", color="#7f7f7f", size=0.5, show.legend = FALSE) +
  coord_fixed(xlim = c(-55, -25), ylim = c(-20, 10)) + 
  geom_point(data = lkl_m$center, mapping= aes(lon, lat), size = 3, stroke = 2, shape = 3, col = 'firebrick', show.legend = FALSE)+
  geom_point(data = lkl_f$center, mapping= aes(lon, lat), size = 3, stroke = 2, shape = 3, col = '#9ECAE1', show.legend = FALSE)+
  geom_point(data = fdn, mapping= aes(longitude, latitude, shape= factor(z), size = 1.5), show.legend = FALSE) +
  geom_text(data=country.label, mapping = aes(longitude, latitude, label = name), fontface="italic")+
  ggtitle("(b) Late chick-rearing to post-breeding")


## PERIOD (c)
start = TIME$date[480]
end = TIME$date[630]
dd <- DATA[which(DATA$time>= start & DATA$time<= end),]
act.M <- dd$act[dd$sex == 'M']/2
act.F <- dd$act[dd$sex == 'F']/2
hm <- table(round(act.M, digits = -1))
hf <- table(round(act.F, digits = -1))
a <- data.frame(M = c(hm/sum(hm)), F = c(hf/sum(hf)))
p1_c <- as.ggplot(~barplot(as.matrix(t(a)), col = c('firebrick', '#9ECAE1'), beside = TRUE, las = 1,
                           xlab = '', ylab = '', ylim = c(0,0.6), 
                           cex.lab = 1) +
                    mtext(text = 'Proportion of time in water (%)', 1, 2.5) + 
                    mtext(text = 'Frequency', 2, 2.5),
                  scale = 1)
dev.M <- dd$lon[which(dd$sex == 'M' & dd$act >150)] -lon.breed
dev.F <- dd$lon[which(dd$sex == 'F' & dd$act >150)] -lon.breed
d.M <- density(dev.M , n = 1024)
d.F <- density(dev.F , n = 1024)
test <- t.test(dev.M, dev.F)
xx.M <-  (d.M$x > -1) & (d.M$x < 4) 
xx.F <-  (d.F$x > -1) & (d.F$x < 4) 
p2_c <- as.ggplot(~plot(d.M$x[xx.M], d.M$y[xx.M], xlim = c(-1, 4), ylim = c(0, 0.7), type = 'l', col = 'firebrick',
                        las = 1, lwd = 3, xlab = "Eastward Deviation from FdN (°) \n when mean wet time >75%", ylab = "Density", cex.lab = 1, mgp=c(3,1,0))+
                    lines(d.F$x[xx.F], d.F$y[xx.F], lwd = 3, col = '#9ECAE1') + 
                    abline(v = 0, lty = 2) +
                    text(2.8, 0.67, paste('p-value :', round(test$p.value, 4)), cex = 1),
                  scale = 0.95)

select <- (DATA$time >= start) & (DATA$time <= end)  & 
  (abs(DATA$diff) <0.5)

data_m <- DATA[ select & (DATA$sex == 'M') ,c("lon", "lat", "act")]
data_f <- DATA[ select & (DATA$sex == 'F') ,c("lon", "lat", "act")]

lkl_m <- get_lkl(data_m)
lkl_f <- get_lkl(data_f)

map_c <- ggplot() +
  geom_polygon(data = eez[-(250:1550),], mapping = aes(long,lat, alpha = 1), fill = "gray",
               colour = "white", linetype = 1, show.legend = FALSE) +
  theme_bw() + theme(legend.position = 'right', 
                     legend.text = element_text(size=10),
                     plot.margin = unit(c(.2,.2,.2,.2), units = "lines"))+  xlab("") + ylab("") +
  geom_contour(aes(x = longitude, y = latitude, z = norm), data = lkl_m$grid,
               size = 1.5,
               colour = adjustcolor('firebrick',red.f=0.75, blue.f=0.75, green.f=0.75),
               na.rm=TRUE, breaks=lkl_m$bk[5], show.legend = FALSE) +
  geom_contour(aes(x = longitude, y = latitude, z = norm), data = lkl_f$grid, 
               size = 1.5,
               colour = adjustcolor('#9ECAE1',red.f=0.75, blue.f=0.75, green.f=0.75),
               na.rm=TRUE, breaks=lkl_m$bk[5], show.legend = FALSE) +
  geom_map(data = map, map=map, aes(map_id=id),
           fill="darkgray", color="#7f7f7f", size=0.5, show.legend = FALSE) +
  coord_fixed(xlim = c(-55, -25), ylim = c(-20, 10)) + 
  geom_point(data = lkl_m$center, mapping= aes(lon, lat, col = 'M'), size = 3, stroke = 2, shape = 3, show.legend = TRUE)+
  geom_point(data = lkl_f$center, mapping= aes(lon, lat, col = 'F'), size = 3, stroke = 2, shape = 3,  show.legend = TRUE)+
  geom_point(data = fdn, mapping= aes(longitude, latitude, shape= factor(z), size = 1.5), show.legend = TRUE) + 
  geom_text(data=country.label, mapping = aes(longitude, latitude, label = name), fontface="italic")+
  ggtitle("(c) Pre-breeding")+
  scale_shape_discrete(
    name = NULL,
    breaks=c("Fernando de Noronha","Atol das Rocas","Sao Pedro Sao Paolo"),
    labels = c("Fernando de Noronha","Atol das Rocas","São Pedro São Paolo")) +
  scale_color_manual(
    name = NULL,
    values = c("#9ECAE1", "firebrick"),
    labels = c("M", "F")) +

  guides(size = FALSE, alpha = FALSE, 
         col = guide_legend(override.aes = list(size = 3, shape = 3, col = c("firebrick", "#9ECAE1"))),
         shape = guide_legend(override.aes = list(size=3, shape=c(17,16,15))))


legend <- g_legend(map_c)

png("./figure/Figure_5.png", width = 1024, height = 768)
blank <- grid.rect(gp=gpar(col="white"))
grid.arrange(map_a, map_b, map_c+theme(legend.position = 'none'), legend, 
             p1_a, p1_b, p1_c, blank, 
             p2_a, p2_b, p2_c, blank, 
             nrow=3, ncol = 4,
             widths = c(2/7, 2/7, 2/7, 1/7), heights = c(4/10, 3/10, 3/10))
dev.off()









