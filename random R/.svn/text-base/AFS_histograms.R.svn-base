getAxis <- function(max, min = 0, tickIncrement = .1)
{
    current = min
    pos = c(min)
    labs = c(min)
    while(current < max){
        current = current+.1
        lab=current
        pos = append(pos, current)
        labs = append(labs,paste(lab,""))
    }
    list(pos,labs)
}

############

pdf(file="/Users/wiliarj/Desktop/nes_rub.pdf",width=8, height=8)
#par(mfrow=c(2,1), mar = c(3,5,1,0))
afs = read.table("/Users/wiliarj/Desktop/fig_data.csv", header=F)
#afs2 = read.table("/Users/wiliarj/Desktop/fig_data2.csv", header=F)
#colors = c("grey","blueviolet", "darkseagreen2","darkseagreen3","darkseagreen4","forestgreen", "darkolivegreen1", "darkolivegreen2")#, "darkolivegreen3", "darkolivegreen4")
#colors = c("grey", "blue", "blue1", "blue2", "blue3", "blue5", "blue6")
#colors = c("grey", "gold2", "slategray1", "skyblue1", "skyblue3", "slateblue1","slateblue3", "slateblue4")
#colors = c("grey", "gold2", "tomato", "tomato3", "red3", "violetred", "violetred2")
#colors = c("slategray1", "skyblue1", "skyblue3", "slateblue1","slateblue3", "slateblue4")
#colors = c("tomato", "tomato3", "red3", "violetred", "violetred2")#, "violetred4")
#colors = c("white","grey20","grey30","grey40","grey50","grey60")
colors = terrain.colors(7)
bcolors = colors
newColors = c()
dens = c()
ang = c()
spc = c()
colors = colors[2:7]
for (c in colors){
    dens = append(dens, c(200, 30))
    ang = append(ang,c(45, 0))
    spc = append(spc, c(.5, .1))
    #s = col2rgb(c)
    #second = s#rgb(s[1]/255, s[2]/255, s[3]/255, 0.3)
    newColors = append(newColors, c(c, c))
}
#colors = newColors
#par(mai=c(1,1,1,1))
#colors = c("white", "blueviolet","grey20","grey30","grey40","grey50")#,"grey60")
#names = c("syn", "ncoding", "CNC")#
names = afs[,1]
#data = afs[,2]
data = afs[,2:length(afs[1,])]
yMax = max(data)
#yMax = 1
yPlace = getAxis(yMax, min=-.1)
yMax = max(yPlace[[1]])
yMax = 1
#bp = barplot(as.matrix(data), beside = T, names.arg = 1:(length(afs[1,])-1),  col = colors, ylim=c(0,yMax), axes=F, space = c(0.4, 1.4))
bp = barplot(1-as.matrix(data), beside = T, names.arg = c("<1", "1-10", "10+"),  col = colors, ylim=c(0,yMax), axes=F, space = c(0.2, 1.4))
#bp = barplot(as.matrix(data), beside = T, names.arg = c("C. grandiflora", "C. rubella"),  col = colors, ylim=c(0,yMax), axes=F)
axis(2, at=yPlace[[1]], labels=yPlace[[2]], cex.axis = 1.25)
#title(main=expression(italic("C. grandiflora")), ylab="Frequency", cex.main = 1, cex.lab = 1, cex.axis = 1)
title(main=expression(italic("")), ylab="Frequency", cex.main = 1.5, cex.lab = 1.5)
#title(main="Alphas for site categories", xlab="Species", ylab="alpha", cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.5)
legend("topright", legend = names, ncol = 1, bty = "n", fill = colors, cex=2)

#par(mar = c(4,5,1,0))
#data = afs2[,2:length(afs2[1,])]
#yMax = max(data)
#yPlace = getAxis(yMax, min=-.1)
#yMax = max(yPlace[[1]])
#yMax = 1
#bp = barplot(as.matrix(data), beside = T, names.arg = 1:(length(afs2[1,])-1),  col = colors, ylim=c(0,yMax), axes=F, space = c(0, 1))
#bp = barplot(as.matrix(data), beside = T, names.arg = c("<1", "1-10", "10+"),  col = colors, ylim=c(0,yMax), axes=F)
#bp = barplot(as.matrix(data), beside = T, names.arg = c("C. grandiflora", "C. rubella"),  col = colors, ylim=c(0,yMax), axes=F)
#axis(2, at=yPlace[[1]], labels=yPlace[[2]])
#title(main=expression(italic("C. rubella")), xlab="Minor allele count", ylab="Frequency", cex.main = 1, cex.lab = 1, cex.axis = 1)
#title(main=expression(italic("C. rubella")), xlab=expression(italic("N"["e"]*"s")*" Category"), ylab="Frequency", cex.main = 1, cex.lab = 1, cex.axis = 1)
#title(main="Alphas for site categories", xlab="Species", ylab="alpha", cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.5)
#legend("topright", legend = names, ncol = 1, bty = "n", fill = colors, cex=1.4)

dev.off()
