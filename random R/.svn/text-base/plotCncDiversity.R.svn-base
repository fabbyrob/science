data = read.table("/Users/wiliarj/Desktop/temp/scaf1_8.flt.cncs", header=T)

data$minDist = ifelse(data$UP_DIST > data$DOWN_DIST, data$DOWN_DIST, data$UP_DIST)
data$piDiv = data$CNC_PI/data$CNC_DIV
data$piDiv = ifelse(is.nan(data$piDiv), NA, data$piDiv)
data$cncSize = data$CNC_END-data$CNC_START

hist(data$cncSize)

par(mfrow=c(3,1), mar=c(1,4,2,3))
plot(data$minDist, data$CNC_DIV, xlim=c(0,10000), xlab="", ylab="divergence", xaxt='n', bty='n')
abline(lm(data$CNC_DIV ~ data$minDist), col="red")
plot(data$minDist, data$CNC_PI, xlim=c(0,10000), xlab="", ylab="pi", xaxt='n', bty='n')
abline(lm(data$CNC_PI ~ data$minDist), col="red")
par(mar=c(5,4,2,3))
plot(data$minDist, data$piDiv, xlim=c(0,10000), xlab="Dist to nearest gene (bases)", ylab="pi/div", bty='n')
abline(lm(data$piDiv ~ data$minDist), col="red")

par(mfrow=c(3,1))
plot(data$UP_DIST, data$CNC_DIV, xlim=c(0,20000), xlab="", ylab="divergence", xaxt='n')
abline(lm(data$CNC_DIV ~ data$UP_DIST), col="red")
plot(data$UP_DIST, data$CNC_PI, xlim=c(0,20000), xlab="", ylab="pi", xaxt='n')
abline(lm(data$CNC_PI ~ data$UP_DIST), col="red")
plot(data$UP_DIST, data$piDiv, xlim=c(0,20000), xlab="Dist to nearest gene in - direction(bases)", ylab="pi/div")
abline(lm(data$piDiv ~ data$UP_DIST), col="red")

par(mfrow=c(3,1))
plot(data$DOWN_DIST, data$CNC_DIV, xlim=c(0,20000), xlab="", ylab="divergence", xaxt='n')
abline(lm(data$CNC_DIV ~ data$DOWN_DIST), col="red")
plot(data$DOWN_DIST, data$CNC_PI, xlim=c(0,20000), xlab="", ylab="pi", xaxt='n')
abline(lm(data$CNC_PI ~ data$DOWN_DIST), col="red")
plot(data$DOWN_DIST, data$piDiv, xlim=c(0,20000), xlab="Dist to nearest gene in + direction(bases)", ylab="pi/div")
abline(lm(data$piDiv ~ data$DOWN_DIST), col="red")

width = 500
data$cuts = cut(data$minDist,breaks=seq(0,10000,width))
boxplot(data$CNC_DIV ~ data$cuts, names = seq(width,10000,width), ylab="div", xlab="", varwidth=T,xaxt='n', bty='n')
boxplot(data$CNC_PI ~ data$cuts, names = seq(width,10000,width), ylab="pi", xlab="", varwidth=T,xaxt='n', bty='n')
boxplot(data$piDiv ~ data$cuts, names = seq(width,10000,width), ylab="pi/div", xlab="dist to nearest gene", varwidth=T, bty='n')