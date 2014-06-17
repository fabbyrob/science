lengths = read.table("/Users/wiliarj/Documents/workspaces/general/robert.williamson/trunk/sandbox/lengths", header = T)
tes = read.table("/Users/wiliarj/Documents/workspaces/general/robert.williamson/trunk/sandbox/tecounts", header = T)

par (mfrow=c(3,1))

yMin = 0
yMax = max(lengths$rubella, lengths$lyrata, lengths$thaliana)

plot(lengths$midpoint, lengths$rubella, col="green", type="l", ylim=c(yMin, yMax))
points(lengths$midpoint, lengths$lyrata, col="red", type="l")
points(lengths$midpoint, lengths$thaliana, col="blue", type="l")
abline(v=1*10^6)
abline(v=2.5*10^6)
abline(v=5.5*10^6)
abline(v=6.75*10^6)
abline(v=9.5*10^6)
abline(v=10*10^6)

yMin = 0
yMax = max(tes$rubella, tes$lyrata, tes$thaliana)
maxTE = yMax
plot(tes$midpoint, tes$rubella, col="green", type="l", ylim=c(yMin, yMax))
points(tes$midpoint, tes$lyrata, col="red", type="l")
points(tes$midpoint, tes$thaliana, col="blue", type="l")
abline(v=1*10^6)
abline(v=2.5*10^6)
abline(v=5.5*10^6)
abline(v=6.75*10^6)
abline(v=9.5*10^6)
abline(v=10*10^6)

propR = tes$rubella/lengths$rubella
propR[is.nan(propR)] = NA
propL = tes$lyrata/lengths$lyrata
propL[is.nan(propL)] = NA
propT = tes$thaliana/lengths$thaliana
propT[is.nan(propT)] = NA

yMin = 0
yMax = max(propR, propL, propT, na.rm=T)
plot(tes$midpoint, propR, col="green", type="l", ylim=c(yMin, yMax))
points(tes$midpoint, propL, col="red", type="l")
points(tes$midpoint, propT, col="blue", type="l")
abline(v=1*10^6)
abline(v=2.5*10^6)
abline(v=5.5*10^6)
abline(v=6.75*10^6)
abline(v=9.5*10^6)
abline(v=10*10^6)

maxTE = 50
x11()
par(mfrow=c(4,3))
for (i in 0:11){
	plot(tes$midpoint, tes[,i+3], col="green", type="l", ylim=c(0, maxTE))
	points(tes$midpoint, tes[,i+27], col="red", type="l", ylim=c(0, maxTE))
	points(tes$midpoint, tes[,i+51], col="blue", type="l", ylim=c(0, maxTE))
	title(main=names(tes)[i+3])
}
title(main="Number of TEs by Family", outer=T)
x11()
par(mfrow=c(4,3))
for (i in 12:22){
	plot(tes$midpoint, tes[,i+3], col="green", type="l", ylim=c(0, maxTE))
	points(tes$midpoint, tes[,i+27], col="red", type="l", ylim=c(0, maxTE))
	points(tes$midpoint, tes[,i+51], col="blue", type="l", ylim=c(0, maxTE))
	title(main=names(tes)[i+3])
}
title(main="Number of TEs by Family", outer=T)

x11()
par(mfrow=c(4,3))
for (i in 0:11){
	plot(tes$midpoint, tes[,i+3]/tes$rubella, col="green", type="l", ylim=c(0, 1))
	points(tes$midpoint, tes[,i+27]/tes$lyrata, col="red", type="l", ylim=c(0, 1))
	points(tes$midpoint, tes[,i+51]/tes$thaliana, col="blue", type="l", ylim=c(0, 1))
}
title(main="Proportion of TEs by Family", outer=T)
x11()
par(mfrow=c(4,3))
for (i in 12:22){
	plot(tes$midpoint, tes[,i+3]/tes$rubella, col="green", type="l", ylim=c(0, 1))
	points(tes$midpoint, tes[,i+27]/tes$lyrata, col="red", type="l", ylim=c(0, 1))
	points(tes$midpoint, tes[,i+51]/tes$thaliana, col="blue", type="l", ylim=c(0, 1))
	title(main=names(tes)[i+3])
}
title(main="Proportion of TEs by Family", outer=T)