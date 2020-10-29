# WaveLightGLS
#
# Figure 2
#----------------------------------


### LOAD DATA

load("./results/SUDA_A33_2017.RData" )

datetime = birdDD$days$tFirst + seconds(difftime(birdDD$days$tSecond, birdDD$days$tFirst, units = "secs")/2)

solar.angle = getElevation(birdDD$days, known.coord = birdGLS$Pos.Deployment, plot = FALSE)
coord = coord(birdDD$days, degElevation = solar.angle)[,1]
act <- birdDD$activity$mean[!is.nan(birdDD$activity$mean)]

png('./figure/Figure_3.png', width = 810, height = 550)

### select window
window = 70:160

layout(matrix(c(1,1,1,5,
                1,1,1,3,
                2,2,2,3,
                2,2,2,4,
                2,2,2,4), byrow = TRUE, ncol = 4))

par(mar = c(0,6,3,10))
plot(datetime[window], coord[window], type = "l", lwd = 3, pch = 15, ylim = c(-35.2, -28), xlab = "", ylab = "",
     xaxt = 'n', yaxt = 'n', xaxs = "i")
axis(2, at = seq(-35,-28, by = 2), labels = seq(-35,-28, by = 2), las =1,
     cex.axis = 1.2)
lines(datetime, act*(7/200)-35, col = "gray22", lty = 3)
axis(4, at = seq(-35,-28, length.out = 5), labels = seq(0,100,25), las =1, col = "gray22", 
     col.axis = "gray22", cex.axis = 1.2)
# abline(h=-40:-20, lty = 2, col = "grey")
mtext("Activity (%)", 4, line = 4, cex = 1, col = "gray22")
mtext("Longitude (°)", 2, line = 4, cex = 1)

M <- list(c(87,67), c(138, 42))
col = c("firebrick", "firebrick")

i = 1
for (m in M) {
  m_low <- max(1, m[1]-10)
  m_up <- min(wc$nc, m[1]+10)
  
  rect(datetime[m_low], -40, datetime[m_up], -20, border = NA, 
       col = adjustcolor(col[i], alpha.f = 0.2), lwd = 2)
  
  i = i+1
}
mtext("(a)", line = -1.5, adj = 1, at = datetime[M[[1]][1]-5])
mtext("(b)", line = -1.5, adj = 1, at = datetime[M[[2]][1]-5])
box(lwd = 2)

par(mar = c(6,6,0.3,6))
image2D( z = t(apply(wc$Power.xy, 2, rev))[window,], x = datetime[window], y = 1:length(wc$Period),
         col = viridis(100, begin = 0, end = 1),
         xlab = "", ylab = "", main = "",
         xaxt = 'n', yaxt = 'n', font = 2)
axis(2, at = which(wc$Period%%2==0), labels =  c(16,8, 4,2,1), las = 1,
     cex.axis = 1.2)
axis(1, at = datetime[seq(5, length(datetime), length.out = 20)],
     labels = as.Date(datetime[seq(5, length(datetime), length.out = 20)]),
     cex.axis = 1.2)
mtext("Time", 1, line = 3, cex = 1)
mtext("Periods (x 24h)", 2, line = 3, cex = 1)
mtext("Cross-Wavelet Power", 4, line = 2, cex = 1, col = "gray22")
contour(z = t(apply(wc$p.value>0.01, 2, rev))[window,], x = datetime[window], y=1:length(wc$Period), add = TRUE,
        nlevels = 1, col = 'gray90', drawlabels = FALSE)
Angle = Arg(wc$Wave.xy)
A.row = seq(1, nrow(Angle), round(nrow(Angle)/30))
A.col = seq(1, ncol(Angle), round(ncol(Angle)/30))
Angle[-A.row, ] = NA
Angle[, -A.col] = NA
for (i in 1:nrow(Angle)) {
  for (j in 1:ncol(Angle)) {
    if (wc$p.value[i,j]<0.05 & !is.na(Angle[i,j]) & j<max(window)){
      x = datetime[j]
      y = wc$nr-i
      arrow(x, y, w = 0.03, alpha = Angle[i,j], col.arrow = "gray90")
    }
  }
}

points(datetime[M[[1]][1]], M[[1]][2], pch = 4, lwd = 3, cex = 2, col =col[1])
points(datetime[M[[2]][1]], M[[2]][2], pch = 4, lwd = 3, cex = 2, col = col[2])

real.morlet.wavelet <-function(t){
  pi ^(-1/4) *  exp(-1/2 * t^2) * cos(omega0 * t)
}

# nights
nights <- birdDD$days[ window, ]

omega0 = 6
title = c("(a)", "(b)")

par(mar = c(0,1,10,3))
m <- M[[1]]
m[2] <- wc$nr - m[2]
m_low <- max(1, m[1]-10)
m_up <- min(wc$nc, m[1]+10)
plot(datetime[m_low:m_up], coord[m_low:m_up], 
     xlab = "", ylab = "", bty = 'n',
     yaxt = "n", xaxt= 'n', main = "",
     type = 'l', ylim = c(-37, -26), col = "black", lwd = 2, cex = 1.5)
lines(datetime[m_low:m_up], act[m_low:m_up]*(8/200)-36,lty = 3, col = "gray22")
scale = wc$Scale[m[2]]
time = datetime[m[1]]
xx <- seq(datetime[m_low], datetime[m_up], length.out = 100)
yy <- seq(m_low, m_up, length.out = 100)
lines(xx, mean(coord[m_low:m_up]) + 10* 1/sqrt(scale) * real.morlet.wavelet((yy - m[1])/scale),
      col = col[2], lwd = 1)
mtext(title[1], line = -2, at = datetime[m_low+2])
for (i in 1:nrow(nights)){
  if(nights$type[i] == 2){
    rect(nights$tFirst[i], -40, nights$tSecond[i], -25, border = NA, density = NA, col = adjustcolor("gray20", alpha.f = 0.3))
  }
}
box(lwd = 2)

par(mar = c(6,1,4,3))
m <- M[[2]]
m[2] <- wc$nr - m[2]
m_low <- max(1, m[1]-10)
m_up <- min(wc$nc, m[1]+10)
plot(datetime[m_low:m_up], coord[m_low:m_up], 
     xlab = "", ylab = "", bty = 'n',
     yaxt = "n", xaxt= 'n', main = "",
     type = 'l', ylim = c(-37, -26), col = "black", lwd = 2, cex = 1.5)
lines(datetime[m_low:m_up], act[m_low:m_up]*(8/200)-36,lty = 3, col = "gray22")
scale = wc$Scale[m[2]]
time = datetime[m[1]]
xx <- seq(datetime[m_low], datetime[m_up], length.out = 100)
yy <- seq(m_low, m_up, length.out = 100)
lines(xx, mean(coord[m_low:m_up]) + 10* 1/sqrt(scale) * real.morlet.wavelet((yy - m[1])/scale),
      col = col[2], lwd = 1)
mtext(title[2], line = -2, at = datetime[m_low+2])
for (i in 1:nrow(nights)){
  if(nights$type[i] == 2){
    rect(nights$tFirst[i], -40, nights$tSecond[i], -25, border = NA, density = NA, col = adjustcolor("gray20", alpha.f = 0.3))
  }
}
box(lwd = 2)


### LEGEND
par(mar = c(0,1,2.5,3))
plot(0,0, col = "white", bty = 'n', xaxt = 'n', yaxt = 'n', xlab = "", ylab = "")
legend("center", 
       c("Longitude (degrees)", "Daily Wet Time (%)", "Morlet Wavelet", "Wavelet Periods", "Nights"),
       pch = c(NA,NA,NA,4, 15),
       lwd = c(2,1,1, 2, 2),
       lty = c(1,3,1, NA, NA),
       col = c("black", "black", "firebrick", "firebrick", "grey"),
       cex = 1.2,
       bty = 'n')


dev.off()