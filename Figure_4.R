# WaveLightGLS
#
# Figure 3
#----------------------------------

load("./results/sex_difference_distance.RData" )

# Read Metadata and select relevant data
metadata = "./data/Metadata_GLS.csv"
metadata <- read.csv(metadata, header=TRUE, sep=",")

list.data = list.files("./results", pattern=".RData", all.files=FALSE,full.names=TRUE)
list.data <- list.data[grepl("SUDA", list.data)]

# selection
idx <- unlist(sapply(metadata$ID[which(metadata$Sex == "F")], function(i){grep(i, list.data)}))
list.sex <- rep("M", length(list.data))
list.sex[idx] <- "F"

START = "2017-05-05"
END = "2018-04-23"

days <- seq(as.Date(START), as.Date(END), by ="day")
TIME <- data.frame( date = rep(days, 2),
                    type = c(rep(1, length(days)), rep(2, length(days)))
                  )
TIME <- TIME[order(TIME$date, TIME$type),]

load(list.data[2])
WC <- array(0, c(nrow(TIME), wc$nr, length(list.data)))
PV <- array(1, c(nrow(TIME), wc$nr, length(list.data)))
ANGLE <- array(0, c(nrow(TIME), wc$nr, length(list.data)))
NB <- array(0, c(nrow(TIME), length(list.data)))

for (data in list.data){
  load(data)
  i = which(data == list.data)
  datetime = birdDD$days$tFirst + seconds(difftime(birdDD$days$tSecond, birdDD$days$tFirst, units = "secs")/2)
  datetime <- datetime[!is.nan(birdDD$activity$mean)]
  
  start = min(datetime)
  start_day = as.Date(start)
  start_type = ifelse(hour(start)<12, 1, 2)
  start_idx = which(TIME$date == start_day & TIME$type == start_type)
  
  end = max(datetime)
  end_day = as.Date(end)
  end_type = ifelse(hour(end)<12, 1, 2)
  end_idx = which(TIME$date == end_day & TIME$type == end_type)
  
  try({
  WC[start_idx:end_idx,,i] <-t(apply(wc$Power.xy, 2, rev))
  ANGLE[start_idx:end_idx,,i] <- Arg(t(apply(wc$Wave.xy, 2, rev)))
  PV[start_idx:end_idx,,i] <- t(apply(wc$p.value, 2, rev))
  NB[start_idx:end_idx, i] <- 1
  })
}

plot(rowSums(NB[,list.sex == "M"]), type = 'l')
lines(rowSums(NB[,list.sex == "F"]), lty = 2)
abline(v = 5, lty = 2, col = "firebrick")
abline(v = 690, lty = 2, col = "firebrick")
legend("topright", c("M", "F"), lty = c(1,2))

time <- TIME[5:690,]
# male
wc_m <- apply(WC[5:690,,list.sex == "M"], MARGIN=c(1, 2), sum)
angle_m <- apply(ANGLE[5:690,,list.sex == "M"], MARGIN=c(1, 2), sum)
pv_m <- apply(PV[5:690,,list.sex == "M"] < 0.01, MARGIN=c(1, 2), sum)
nb_m <-  apply(NB[5:690,list.sex == "M"], MARGIN=1, sum)
wm = wc_m/nb_m
am = angle_m/nb_m

# female
wc_f <- apply(WC[5:690,,list.sex == "F"], MARGIN=c(1, 2), sum)
angle_f <- apply(ANGLE[5:690,,list.sex == "F"], MARGIN=c(1, 2), sum)
pv_f <- apply(PV[5:690,,list.sex == "F"] < 0.01, MARGIN=c(1, 2), sum)
nb_f <-  apply(NB[5:690,list.sex == "F"], MARGIN=1, sum)
wf = wc_f/nb_f
af = angle_f/nb_f


png('./figure/Figure_4.png', width = 600, height = 756)
### AVERAGED CROSS-POWER
layout(matrix(c(1,1,1,
                1,1,1,
                2,2,2,
                2,2,2,
                3,3,3,
                4,5,6), ncol = 3, byrow = TRUE))

## FEMALE
par(mar = c(0,5,1,10))
image( z = wf, x = 1:nrow(time), y = 1:length(wc$Period),
         col = viridis(20, begin = 0, end = 1),
         xlab = "", ylab = "", main = "",
         xaxt = 'n', yaxt = 'n',
       breaks = seq(0,0.5, length.out = 21))
axis(2, at = which(wc$Period%%2==0), labels =  c(16, 8, 4, 2, 1), las = 1,
     cex.axis = 1.2)
mtext("Periods (x 24h)", 2, line = 3, cex = 1)
mtext("Averaged \n Cross-Wavelet Power", 4, line = 3, cex = 1)
contour(z = pv_f>0, x = 1:nrow(time), y=1:length(wc$Period), add = TRUE,
        nlevels = 1, col = 'gray90', drawlabels = FALSE)
A.row = seq(1, nrow(af), round(nrow(af)/30))
A.col = seq(1, ncol(af), round(ncol(af)/30))
af[-A.row, ] = NA
af[, -A.col] = NA
for (i in 1:nrow(af)) {
  for (j in 1:ncol(af)) {
    if (pv_f[i,j]>0 & !is.na(af[i,j])){
      arrow(i, j, w = 0.03, alpha = af[i,j], col.arrow = "gray90")
    }
  }
}
### COI
coi = COI(1, 1, nrow(time), ncol(WC), wc$Period)
xx = coi$x
yy = wc$nr-coi$y/max(coi$y)*wc$nr
yy[1] = 0
yy[length(yy)] = 0
polygon(coi$x, yy, border = 'white', col = rgb(1,1,1, alpha = 0.6))
text(labels = "(F)", x = 25, y = 75, col = "grey", cex = 2, lwd = 2)


## MALE
par(mar = c(0,5,1,5))
image2D( z = wm, x = 1:nrow(time), y = 1:length(wc$Period),
         col = viridis(20, begin = 0, end = 1),
         xlab = "", ylab = "", main = "",
         xaxt = 'n', yaxt = 'n',
         breaks = seq(0,0.5, length.out = 21))
axis(2, at = which(wc$Period%%2==0), labels =  c(16,8,4,2,1), las = 1,
     cex.axis = 1.2)
mtext("Periods (x 24h)", 2, line = 3, cex = 1)
contour(z = pv_m>0, x = 1:nrow(time), y=1:length(wc$Period), add = TRUE,
        nlevels = 1, col = 'gray90', drawlabels = FALSE)
A.row = seq(1, nrow(am), round(nrow(am)/30))
A.col = seq(1, ncol(am), round(ncol(am)/30))
am[-A.row, ] = NA
am[, -A.col] = NA
for (i in 1:nrow(af)) {
  for (j in 1:ncol(af)) {
    if (pv_m[i,j]>0 & !is.na(am[i,j])){
      arrow(i, j, w = 0.03, alpha = am[i,j], col.arrow = "gray90")
    }
  }
}
### COI
coi = COI(1, 1, nrow(time), ncol(WC), wc$Period)
xx = coi$x
yy = wc$nr-coi$y/max(coi$y)*wc$nr
yy[1] = 0
yy[length(yy)] = 0
polygon(coi$x, yy, border = 'white', col = rgb(1,1,1, alpha = 0.6))
text(labels = "(M)", x = 25, y = 75, col = "grey", cex = 2, lwd = 2)


## phenology masked booby
laying <- as.Date(c("2017-03-15", "2017-05-10", "2018-03-15", "2018-05-10"))
incubating <- as.Date(c("2017-03-15", "2017-06-15", "2018-03-15", "2018-06-15"))
hatching <- as.Date(c("2017-05-01", "2017-06-15", "2018-05-01", "2018-06-15"))  # éclosion
rearing <- as.Date(c("2017-05-01", "2017-11-10", "2018-05-01", "2018-11-10"))   # élevage
fledging <- as.Date(c("2017-09-01", "2017-11-10", "2018-09-01", "2018-11-10"))  # envol


par(mar = c(7,5,1,10))
plot(days,rep(-1, length(days)), ylim = c(0,5),  
     bty = 'n',
     xaxs = 'i', yaxs = 'i', xaxt = 'n', yaxt = 'n', ylab = "", xlab = "")
axis(4, at = c(0.5, 1.5, 2.5, 3.5, 4.5), las = 1, cex.axis = 1, lwd = 0, lwd.ticks = 0,
     labels = c("laying", "incubating", "hatching", "rearing", "fledging"))
axis(1, at=range(days), labels=c("",""), lwd.ticks=0)
axis(1, at = days[which(day(days) == 1)],
     labels = month.abb[month(days[which(day(days) == 1)])],
     cex.axis = 1.5, lwd=0, lwd.ticks=1)
arrows(TIME$date[110], -3, TIME$date[200], -3, lwd = 2, col = "firebrick", code = 3, length = 0.1, xpd = TRUE)
arrows(TIME$date[230], -3, TIME$date[450], -3, lwd = 2, col = "grey", code = 3, length = 0.1, xpd = TRUE)
arrows(TIME$date[480], -3, TIME$date[630], -3, lwd = 2, col = "gray22", code = 3, length = 0.1, xpd = TRUE)
text(labels = c("(a)", "(b)", "(c)"), 
     x = TIME$date[c(310/2, 680/2, 1110/2)],
     y = -4,
     cex = 1.2,
     col = c("firebrick", "grey", "gray22"),
     xpd = TRUE)
#laying
rect(laying[1],0, laying[2], 1, density = NULL, border = NA, col =  adjustcolor("#AB1F17", alpha.f = 0.4))
rect(laying[3],0, laying[4], 1, density = NULL, border = NA, col =  adjustcolor("#AB1F17", alpha.f = 0.4))
#incubating
rect(incubating[1],1, incubating[2], 2, density = NULL, border = NA, col = adjustcolor("#AB1F17", alpha.f = 0.55))
rect(incubating[3],1, incubating[4], 2, density = NULL, border = NA, col = adjustcolor("#AB1F17", alpha.f = 0.55))
#hatching
rect(hatching[1],2, hatching[2], 3, density = NULL, border = NA, col =  adjustcolor("#AB1F17", alpha.f = 0.7))
rect(hatching[3],2, hatching[4], 3, density = NULL, border = NA, col =  adjustcolor("#AB1F17", alpha.f = 0.7))
#rearing
rect(rearing[1],3, rearing[2], 4, density = NULL, border = NA, col =  adjustcolor("#AB1F17", alpha.f = 0.85))
rect(rearing[3],3, rearing[4], 4, density = NULL, border = NA, col =  adjustcolor("#AB1F17", alpha.f = 0.85))
#fledging
rect(fledging[1],4, fledging[2], 5, density = NULL, border = NA, col =  adjustcolor("#AB1F17", alpha.f = 1))
rect(fledging[3],4, fledging[4], 5, density = NULL, border = NA, col =  adjustcolor("#AB1F17", alpha.f = 1))


par(mar = c(5,4.5,0,4))
plot(rev(wc$Period)/2, colMeans(wf[110:200,], na.rm = TRUE), xlab = "", ylab = "", 
     type = 'l', col = "firebrick", lwd = 2, ylim = c(0, 0.3))
lines(rev(wc$Period)/2, colMeans(wm[110:200,], na.rm = TRUE), col = "firebrick", lwd = 2, lty = 2)
mtext("Periods", 1, line =3)
mtext("Averaged Power", 2, line =3)
mtext("(a)", 3, line =1)
# abline(v = wc$Period[wc$nr-which.max(colMeans(wf[110:200,]))], lty = 2)

plot(rev(wc$Period)/2, colMeans(wf[230:450,], na.rm = TRUE), xlab = "", ylab = "", 
     type = 'l', col = "grey", lwd = 2, ylim = c(0, 0.3))
lines(rev(wc$Period)/2, colMeans(wm[230:450,], na.rm = TRUE), col = "grey", lwd = 2, lty = 2)
mtext("Periods", 1, line =3)
mtext("Averaged Power", 2, line =3)
mtext("(b)", 3, line =1)
# abline(v = wc$Period[wc$nr-which.max(colMeans(wf[220:430,]))], lty = 2)

plot(rev(wc$Period)/2, colMeans(wf[480:630,], na.rm = TRUE), xlab = "", ylab = "", 
     type = 'l', col = "gray22", lwd = 2, ylim = c(0, 0.3))
lines(rev(wc$Period)/2, colMeans(wm[480:630,], na.rm = TRUE), col = "gray22", lwd = 2, lty = 2)
mtext("Periods", 1, line =3)
mtext("Averaged Power", 2, line =3)
mtext("(c)", 3, line =1)
legend("topright", c("F", "M"), lwd = 2, lty = c(1,2), bty = "n")
# abline(v = wc$Period[wc$nr-which.max(colMeans(wf[460:620,]))], lty = 2)
  
# par(mar = c(5,2,0,4))
# plot(d.M$x, d.M$y, xlim = c(-1.5, 5), ylim = c(0, 0.5), lwd = 2, type = 'l', lty = 2, 
#      las = 1, xlab = "", ylab = "",
#      main = "")
# lines(d.F$x, d.F$y, lwd = 2)
# legend("topright", c("F", "M"), lwd = 2, lty = c(1,2), bty = "n")
# mtext("Eastward Deviation (degrees)", 1, line =3)
# mtext("Density", 2, line =3)
# mtext("(c')", 3, line =1)

dev.off()