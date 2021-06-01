# WaveLightGLS
#
# Figure 1 & 2
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
  
  if(nrow(birdDD$temperature)+1 == nrow(crds$x)){
    temp <- c(NA,birdDD$temperature$mean)
  } else{
    temp <- rep(NA, nrow(crds$x))
  }
  
  out <- data.frame(bird = birdGLS$ID, time = crds$time, zenithT = median(z), zenith = max(z), 
                    tw_error = twl_dev_all, lon = crds$x[,1], lat =  crds$x[,2],
                    act = act, temp = temp)
  DATA <- rbind(DATA, out)
  
  cat(i, ' out of ', length(list.data), '\n')
  i = i+1
}


### TEMPERATURE ERROR AND DEVIATIONS

DATA$temp_fdn_sat <- getSSTPoint(path = "./data/METOFFICE-GLO-SST-L4-REP-OBS-SST_1590497114049.nc",
                                coord = matrix(rep(c(lon.breed, lat.breed), nrow(DATA)), ncol = 2, byrow = TRUE),
                                time = DATA$time) -273.15

DATA$temp_sat <-  getSSTPoint(path = "./data/METOFFICE-GLO-SST-L4-REP-OBS-SST_1590497114049.nc",
                                  coord = DATA[,c("lon", "lat")],
                                  time = DATA$time) -273.15

diff <- DATA$temp - DATA$temp_fdn_sat
sel <- DATA$time >= calib.tm[1] & DATA$time <= calib.tm[2]
nights <- hour(DATA$time) < 12

# temp_dt <- getSSTPoint(path = "./data/METOFFICE-GLO-SST-L4-REP-OBS-SST_1590497114049.nc", 
#             coord = matrix(rep(c(lon.breed, lat.breed), 354), ncol = 2, byrow = TRUE), 
#             time = seq(min(DATA$time), max(DATA$time), by = 'days')) -273.15
# 
# plot(temp_dt)



### FIGURE 1

png('./figure/Figure_1.png', width = 760, height = 750)
par(mfrow = c(3,2), mar = c(5,4,1,2))
dev <- CALIBRATION$twl_dev
hist(dev, xlim = c(-15, 40), ylim = c(0, 0.1), breaks = seq(-500, 1500, by = 2.5), freq = F,
     main = "", col = "grey",
     xlab = "", ylab ="density")
mtext(side=1, line=2, at=5, adj=0, cex=0.8, "(minutes)")
mtext(side=3, line=-2, at=40, adj=1, cex=1, "(a) Twilight Deviation")
mtext(side=3, line=-3.5, at=40, adj=1, cex=0.9, "Calibration period")
gamma <- unique(CALIBRATION[,c("shape", "scale", "model")])
# for ( i in 1:nrow(gamma)){
#   xx <- seq(0, 40, by = 0.1)
#   yy <- dgamma(xx, gamma$shape[i], gamma$scale[i])
#   lines(xx, yy, col = "grey", lty = 2)
#   i = i+1
# }
seq <- seq(0, 40, length = 100)
fit_g <- fitdistr(dev, "gamma")
lines(seq, dgamma(seq, fit_g$estimate[1], fit_g$estimate[2]), col = "firebrick", lwd = 2.5, lty = 2)

seq_ = seq(-2,2,by = 0.01)
hist(diff[sel], freq = F, xlim = c(-3, 3), ylim = c(0,2), breaks = seq(-10, 10, by = 0.25),
     main = "", col = "grey",
     xlab = "", ylab ="density")
mtext(side=1, line=2, at=0, adj=0, cex=0.8, "(celsius)")
mtext(side=3, line=-2, at=3, adj=1, cex=1, "(b) Temperature Deviation")
mtext(side=3, line=-3.5, at=3, adj=1, cex=0.9, "Calibration period")
lines(seq_, 1/1.75*(seq_>-0.25&seq_<1.5), col = "firebrick", lwd = 2.5, lty = 2)


hist(DATA$tw_error, xlim = c(-15, 40), ylim = c(0, 0.1), breaks = seq(-500, 1500, by = 2.5), freq = F,
     main = "", col = "#9ECAE1",
     xlab = "", ylab="density")
mtext(side=1, line=2, at=5, adj=0, cex=0.8, "(minutes)")
mtext(side=3, line=-2, at=40, adj=1, cex=1, "(c) Twilight Deviation")
mtext(side=3, line=-3.5, at=40, adj=1, cex=0.9, "Year-round data")
lines(seq, dgamma(seq, fit_g$estimate[1], fit_g$estimate[2]), col = "firebrick", lwd = 2.5, lty = 2)

hist(diff, freq = F, xlim = c(-3, 3), ylim = c(0,2), breaks = seq(-10, 10, by = 0.25),
     main = "", col = "#9ECAE1",
     xlab = "", ylab="density")
mtext(side=3, line=-2, at=3, adj=1, cex=1, "(d) Temperature Deviation")
mtext(side=3, line=-3.5, at=3, adj=1, cex=0.9, "Year-round data")
mtext(side=1, line=2, at=0, adj=0, cex=0.8, "(celsius)")
lines(seq_, 1/1.75*(seq_>-0.25&seq_<1.5), col = "firebrick", lwd = 2.5, lty = 2)

### HISTOGRAMS WITH HIGH ACTIVITY
hist(DATA$tw_error[DATA$act>150], xlim = c(-15, 40), ylim = c(0, 0.1), breaks = seq(-500, 1500, by = 2.5), freq = F,
     main = "", col = "#FDAE6B",
     xlab = "", ylab="density")
mtext(side=1, line=2, at=5, adj=0, cex=0.8, "(minutes)")
mtext(side=3, line=-2, at=40, adj=1, cex=1, "(e) Twilight Deviation")
mtext(side=3, line=-3.5, at=40, adj=1, cex=0.9, "Time spent in water >75%")
lines(seq, dgamma(seq, fit_g$estimate[1], fit_g$estimate[2]), col = "firebrick", lwd = 2.5, lty = 2)

hist(diff[DATA$act>150], freq = F, xlim = c(-3, 3), ylim = c(0,2), breaks = seq(-10, 10, by = 0.25),
     main = "", col = "#FDAE6B",
     xlab = "", ylab="density")
mtext(side=3, line=-2, at=3, adj=1, cex=1, "(f) Temperature Deviation")
mtext(side=3, line=-3.5, at=3, adj=1, cex=0.9, "Time spent in water >75%")
mtext(side=1, line=2, at=0, adj=0, cex=0.8, "(celsius)")
lines(seq_, 1/1.75*(seq_>-0.25&seq_<1.5), col = "firebrick", lwd = 2.5, lty = 2)
dev.off()

### FIGURE 2
# load graphic data
data(wrld_simpl)
wrld_simpl@data$id <- wrld_simpl@data$NAME
world <- fortify(wrld_simpl)

eez<-readOGR("./data/World_EEZ.shp", "World_EEZ")
EEZ <-  fortify(eez)
EEZ_br <- EEZ[which(EEZ$id==163 & !EEZ$hole),]

### YEAR-ROUND ERROR RANGE AT FDN
days <- seq(min(as.Date(DATA$time)), max(as.Date(DATA$time)), by = "days")
days_rise <- twilight(days, lon.breed, lat.breed,
                      rise = TRUE, zenith = 96, iters = 3)
days_fall <- twilight(days, lon.breed, lat.breed, 
                      rise = FALSE, zenith = 96, iters = 3)

twilights <- data.frame(Twilight = c(days_rise, days_fall),
                        Rise = c(rep(TRUE, length(days_rise)), rep(FALSE, length(days_fall))))
twilights <- twilights[order(twilights$Twilight),]


COORD <- NULL
for (k in 1:100){
  tw <- twilights
  tw$Twilight[tw$Rise] <- tw$Twilight[tw$Rise] + seconds(round( 60*(rgamma(sum(tw$Rise), fit_g$estimate[1], fit_g$estimate[2]))))
  tw$Twilight[!tw$Rise] <- tw$Twilight[!tw$Rise] - seconds(round( 60*(rgamma(sum(!tw$Rise), fit_g$estimate[1], fit_g$estimate[2]))))
  
  zenith <- calib(tw, lat.breed, 96)
  crds <- thresholdPath(tw$Twilight, tw$Rise, zenith = zenith)
  out <- data.frame(lon = crds$x[,1], lat = crds$x[,2], time = crds$time)
  
  out$temp_sat <- getSSTPoint(path = "./data/METOFFICE-GLO-SST-L4-REP-OBS-SST_1590497114049.nc",
                              coord = crds$x,
                              time = crds$time) -273.15
  ### TEMPERATURE ERROR AND DEVIATIONS
  out$temp_fdn_sat <-  getSSTPoint(path = "./data/METOFFICE-GLO-SST-L4-REP-OBS-SST_1590497114049.nc",
                                    coord = matrix(rep(c(lon.breed, lat.breed), nrow(out)), byrow = TRUE, ncol = 2),
                                    time = crds$time) -273.15
  COORD <- rbind(COORD, out)
}


diff <- DATA$temp - DATA$temp_sat

map1_th <- plot.kde.coord(COORD[,c("lon", "lat")], H=2, N=100, alpha = 0, eez = EEZ_br, 
                          title = "(a) Error Range Estimation", col = "firebrick")
map2_th <- plot.kde.coord(COORD[abs(COORD$temp_sat-COORD$temp_fdn_sat)<=0.5,c("lon", "lat")],
                          H=2, N=100, alpha = 0, eez = EEZ_br,
                          title = "(b) Error Range Estimation", col = "firebrick")

map1 <- plot.kde.coord(DATA[,c("lon", "lat")], H=2, N=100, alpha = 0.1, eez = EEZ_br,
                       title = "(c) Positions Distribution", col = "#9ECAE1")
map2 <- plot.kde.coord(DATA[ diff<=1.5 & diff >= -0.25,c("lon", "lat")], H=2, N=100, alpha = 0.1, eez = EEZ_br,
                       title = "(d) Positions Distribution", col = "#9ECAE1")


map1_act <- plot.kde.coord(DATA[which(DATA$act > 150),c("lon", "lat")], H=2, N=100, alpha = 0.1,
                           eez = EEZ_br, title = "(e) Wet Positions Distribution", col = "#FDAE6B")
map2_act <- plot.kde.coord(DATA[which(DATA$act > 150 & diff<=1.5 & diff >= -0.25),c("lon", "lat")], H=2, N=100,
                           alpha = 0.1, eez = EEZ_br, 
                           title = "(f) Wet Positions Distribution", col = "#FDAE6B")

legend <- g_legend(map2_act)

png("./figure/Figure_2.png", width = 1270, height = 800)
grid.arrange(map1_th + theme(legend.position = 'none'),
             map1+ theme(legend.position = 'none'),
             map1_act+ theme(legend.position = 'none'),
             legend,
             map2_th+ theme(legend.position = 'none'),
             map2+ theme(legend.position = 'none'),
             map2_act+ theme(legend.position = 'none'), ncol=4, nrow = 2, 
             widths = c(2/7, 2/7, 2/7, 1/7))
dev.off()



## Bhattacharyya coefficient
N = 100
H = 2
CRDS <- COORD[,c("lon", "lat")]
# CRDS <- COORD[abs(COORD$temp_sat-COORD$temp_fdn_sat)<=0.5,c("lon", "lat")]
CRDS <- CRDS[!is.na(rowSums(CRDS)),]
f1 <- with(CRDS, kde2d(CRDS[,1], CRDS[,2], n = N, h = H, lims = c(-60, 0, -30, 20)))
CRDS <- DATA[,c("lon", "lat")]
# CRDS <- DATA[abs(diff)<=0.5,c("lon", "lat")]
CRDS <- CRDS[!is.na(rowSums(CRDS)),]
f2 <- with(CRDS, kde2d(CRDS[,1], CRDS[,2], n = N, h = H, lims = c(-60, 0, -30, 20)))

sum(sqrt(f1$z*f2$z/sum(f1$z)/sum(f2$z)))


###
R = 6378

mean(abs(pi * R * (COORD$lon - lon.breed) / 180))
sd(abs(pi * R * (COORD$lon - lon.breed) / 180))


mean(abs(pi * R * (COORD$lat - lat.breed) / 180))
sd(abs(pi * R * (COORD$lat - lat.breed) / 180))



mean(abs(pi * R * (COORD$lon[abs(COORD$temp_sat-COORD$temp_fdn_sat)<=0.5] - lon.breed) / 180), na.rm = TRUE)
sd(abs(pi * R * (COORD$lon[abs(COORD$temp_sat-COORD$temp_fdn_sat)<=0.5] - lon.breed) / 180), na.rm = TRUE)


mean(abs(pi * R * (COORD$lat[abs(COORD$temp_sat-COORD$temp_fdn_sat)<=0.5] - lat.breed) / 180), na.rm = TRUE)
sd(abs(pi * R * (COORD$lat[abs(COORD$temp_sat-COORD$temp_fdn_sat)<=0.5] - lat.breed) / 180), na.rm = TRUE)


