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

diff <- DATA$temp_sat - DATA$temp_fdn_sat
sel <- which(DATA$time >= calib.tm[1] & DATA$time <= calib.tm[2])

### FIGURE 1

png('./figure/Figure_1.png', width = 760, height = 750)
par(mfrow = c(3,2), mar = c(5,4,1,2))
dev <- CALIBRATION$twl_dev
hist(dev, xlim = c(-15, 40), ylim = c(0, 0.1), breaks = seq(-500, 1500, by = 2.5), freq = F,
     main = "", col = "grey",
     xlab = "", ylab ="density")
mtext(side=1, line=2, at=5, adj=0, cex=0.8, "(minutes)")
mtext(side=3, line=-2, at=40, adj=1, cex=1, "Twilight Deviation")
mtext(side=3, line=-3, at=40, adj=1, cex=0.8, "Calibration period")
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
hist(diff[sel], freq = F, xlim = c(-3, 3), ylim = c(0,1.5), breaks = seq(-5, 5, by = 0.25),
     main = "", col = "grey",
     xlab = "", ylab ="density")
mtext(side=1, line=2, at=0, adj=0, cex=0.8, "(celsius)")
mtext(side=3, line=-2, at=3, adj=1, cex=1, "Temperature Deviation")
mtext(side=3, line=-3, at=3, adj=1, cex=0.8, "Calibration period")
lines(seq_, 1*(abs(seq_)<0.5), col = "firebrick", lwd = 2.5, lty = 2)

hist(DATA$tw_error, xlim = c(-15, 40), ylim = c(0, 0.1), breaks = seq(-500, 1500, by = 2.5), freq = F,
     main = "", col = "#9ECAE1",
     xlab = "", ylab="density")
mtext(side=1, line=2, at=5, adj=0, cex=0.8, "(minutes)")
mtext(side=3, line=-2, at=40, adj=1, cex=1, "Twilight Deviation")
mtext(side=3, line=-3, at=40, adj=1, cex=0.8, "Year-round data")
lines(seq, dgamma(seq, fit_g$estimate[1], fit_g$estimate[2]), col = "firebrick", lwd = 2.5, lty = 2)

hist(diff, freq = F, xlim = c(-3, 3), ylim = c(0,1.5), breaks = seq(-10, 10, by = 0.25),
     main = "", col = "#9ECAE1",
     xlab = "", ylab="density")
mtext(side=3, line=-2, at=3, adj=1, cex=1, "Temperature Deviation")
mtext(side=3, line=-3, at=3, adj=1, cex=0.8, "Year-round data")
mtext(side=1, line=2, at=0, adj=0, cex=0.8, "(celsius)")
lines(seq_, 1*(abs(seq_)<0.5), col = "firebrick", lwd = 2.5, lty = 2)

### HISTOGRAMS WITH HIGH ACTIVITY
hist(DATA$tw_error[DATA$act>150], xlim = c(-15, 40), ylim = c(0, 0.1), breaks = seq(-500, 1500, by = 2.5), freq = F,
     main = "", col = "#FDAE6B",
     xlab = "", ylab="density")
mtext(side=1, line=2, at=5, adj=0, cex=0.8, "(minutes)")
mtext(side=3, line=-2, at=40, adj=1, cex=1, "Twilight Deviation")
mtext(side=3, line=-3, at=40, adj=1, cex=0.8, "Time spent in water >75%")
lines(seq, dgamma(seq, fit_g$estimate[1], fit_g$estimate[2]), col = "firebrick", lwd = 2.5, lty = 2)

hist(diff[DATA$act>150], freq = F, xlim = c(-3, 3), ylim = c(0,1.5), breaks = seq(-10, 10, by = 0.25),
     main = "", col = "#FDAE6B",
     xlab = "", ylab="density")
mtext(side=3, line=-2, at=3, adj=1, cex=1, "Temperature Deviation")
mtext(side=3, line=-3, at=3, adj=1, cex=0.8, "Time spent in water >75%")
mtext(side=1, line=2, at=0, adj=0, cex=0.8, "(celsius)")
lines(seq_, 1*(abs(seq_)<0.5), col = "firebrick", lwd = 2.5, lty = 2)
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


map1_th <- plot.kde.coord(COORD[,c("lon", "lat")], H=2, N=100, alpha = 0, eez = EEZ_br, 
                          title = "(a) Error Range Estimation", col = "firebrick")
map2_th <- plot.kde.coord(COORD[abs(COORD$temp_sat-COORD$temp_fdn_sat)<=0.5,c("lon", "lat")],
                          H=2, N=100, alpha = 0, eez = EEZ_br,
                          title = "(b) Error Range Estimation", col = "firebrick")

map1 <- plot.kde.coord(DATA[,c("lon", "lat")], H=2, N=100, alpha = 0.1, eez = EEZ_br,
                       title = "(a) Habitat Estimation", col = "#9ECAE1")
map2 <- plot.kde.coord(DATA[abs(diff)<=0.5,c("lon", "lat")], H=2, N=100, alpha = 0.1, eez = EEZ_br,
                       title = "(b) Habitat Estimation", col = "#9ECAE1")


map1_act <- plot.kde.coord(DATA[which(DATA$act > 150),c("lon", "lat")], H=2, N=100, alpha = 0.1,
                           eez = EEZ_br, title = "(a) Wet Habitat Estimation", col = "#FDAE6B")
map2_act <- plot.kde.coord(DATA[which(DATA$act > 150 & abs(diff)<=0.5),c("lon", "lat")], H=2, N=100,
                           alpha = 0.1, eez = EEZ_br, 
                           title = "(b) Wet Habitat Estimation", col = "#FDAE6B")

png("./figure/Figure_2.png", width = 1270, height = 800)
grid.arrange(map1_th, map1, map1_act, map2_th, map2, map2_act, ncol=3, nrow = 2)
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


## mean distance to the colony --- For Figure 4

# Read Metadata and select relevant data
options(stringsAsFactors = FALSE)
metadata = "./data/Metadata_GLS.csv"
metadata <- read.csv(metadata, header=TRUE, sep=",")
# selection
DATA$sex <- sapply(DATA$bird, function(b){metadata$Sex[metadata$ID == b]})

lon.dev.M <- (DATA[which(DATA$act > 150 & DATA$sex == "M" & 
                        DATA$time>= as.Date("2018-01-01") &
                        DATA$time<= as.Date("2018-03-15")),c("lon")])
lon.dev.F <- (DATA[which(DATA$act > 150 & DATA$sex == "F" & 
                           DATA$time>= as.Date("2018-01-01") &
                           DATA$time<= as.Date("2018-03-15")),c("lon")])
# orthodromic distance
distance_ortho_robuste<-function(lat1m,lat2m,lon1m,lon2m){
  R <- 6377726
  dist.m<-R*2*asin(
    ((sin((lat1m*pi/180-lat2m*pi/180)/2))^2+
       cos(lat1m*pi/180)*cos(lat2m*pi/180)*(sin((lon1m*pi/180-lon2m*pi/180)/2))^2)^.5)
  return(dist.m)
  end
}


d.M <- 
  sapply(lon.dev.M, function(l){
  distance_ortho_robuste(lat.breed, lat.breed, lon.breed, l)
})

d.F <- 
  sapply(lon.dev.F, function(l){
    distance_ortho_robuste(lat.breed, lat.breed, lon.breed, l)
  })


HM <- hist(d.M/1000, breaks = seq(0, 5000, by = 30), plot = FALSE)
HF <- hist(d.F/1000, breaks = seq(0, 5000, by = 30), plot = FALSE)

plot(HF, col = adjustcolor("pink", alpha.f = 0.2), xlim = c(0, 700))
plot(HM, col = adjustcolor("blue", alpha.f = 0.2), add = TRUE)

median(d.M)
median(d.F)

d.M <- density(lon.dev.M - lon.breed, n = 1024)
d.F <- density(lon.dev.F - lon.breed, n = 1024)

plot(d.M$x, d.M$y, xlim = c(-1.5, 5), ylim = c(0, 0.5), lwd = 2, type = 'l', lty = 2, 
     las = 1, xlab = "Eastward Deviation from FdN (degrees)", ylab = "Density",
     main = "Positions Estimates of days in high wet environments( >75%)")
lines(d.F$x, d.F$y, lwd = 2)
legend("topright", c("F", "M"), lwd = 2, lty = c(1,2), bty = "n")

t.test(lon.dev.M- lon.breed, lon.dev.F- lon.breed)

save(d.M, d.M, file = paste0("./results/sex_difference_distance.RData"))
