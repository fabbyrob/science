type = "cncsites"
base = "/Users/wiliarj/Desktop/temp/"
data = read.table(paste(base,type,".cis", sep=""),header=T)
library(plotrix)
colors = grey.colors(length(data$type))
#colors = c("black", "darkgoldenrod", "darkgoldenrod1", "darkgoldenrod4", "darkolivegreen1", "darkolivegreen3", "cyan2", "cyan3", "cyan4", "darkorange", "darkorange3", "azure2", "azure3", "azure4")
#plot pos selection

pdf(file=paste(base,"posselection_",type,".pdf", sep=""))#, height=12, width=12)
par(mfrow=c(2,1), mar=c(1,5,2,3))
#plot alphas
plotCI(seq(data$alpha), y=data$alpha, pch=20, bty='n', li=data$alpha_min, ui=data$alpha_max, ylab=expression("" ~ alpha), xlab="Species", xaxt="n")
#axis(1, at=seq(data$alpha), labels=data$type)

#plot omegas
#plot alphas
par(mar=c(5,5,1,3))
plotCI(seq(data$alpha), y=data$omega, pch=20, bty='n', li=data$omega_min, ui=data$omega_max, ylab=expression("" ~ omega), xlab="Species", xaxt="n")
axis(1, at=seq(data$alpha), labels=data$type)

dev.off()

#plot nes
pdf(file=paste(base,"negselection_",type,".pdf", sep=""))
par(mfrow=c(1,1), mar=c(5,5,0,3))
bp = barplot(as.matrix(data[,c(9,12,15,18)]), beside = T, names.arg = c("<1", "1-10", "10-100", "100+"),xlab = expression(paste(italic('N'['e']), italic('s'), " category")), ylab = "Fraction of sites", col=colors, ylim=c(0,1.1), space = c(0.2, 1.4), axes=F)

for (i in c(0:3)){
    arrows(bp[,1+i], data[,9+i*3], bp[,1+i], data[,10+i*3], angle = 90, length = 0.05)
    arrows(bp[,1+i], data[,9+i*3], bp[,1+i], data[,8+i*3], angle = 90, length = 0.05)
}

legend("topright", legend = data$type, ncol = 2, bty = "n", fill = colors, cex=1.25)
axis(2, at=seq(0,1,0.2), labels=seq(0,1,0.2), cex.axis = 1.25)

dev.off()