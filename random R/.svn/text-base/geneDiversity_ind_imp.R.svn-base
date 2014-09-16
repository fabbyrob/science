####FUNCTIONS
plotMeans <- function(dataSet, classes, col = "red", lwd=1){
    mean.data = tapply(dataSet, classes, mean)
    lines(mean.data, col=col, type = "l",lwd=lwd)
}


pickRandomGenes <- function(numGenes, genes){
    indxs = sample(1:length(genes), numGenes, replace=F)
    myGenes = c()
    for (i in indxs){
        myGenes = append(myGenes, as.vector(genes[i]))
    }
    return (myGenes)
}

#########

path = '/Users/wiliarj/Desktop/temp/'
#path = 'C:/Users/Fabby/Desktop/temp/'

data = read.csv(paste(path, 'total.mod.out',sep=""),header=T)
data2 = read.csv(paste(path, 'all.mod.out',sep=""),header=T)
names(data) = c("gene", "mid", "pi", "theta", "tajd")
names(data2) = c("gene", "mid", "pi", "theta", "tajd")

flanks = 500
allcol = "black"
cols2 = "blue"
first = 1
xrange = c(-1*flanks, 2*flanks)

boxcolor = "lightgrey"
cut = 50

pdfer = 1
meanx = seq(-1*flanks+0.5*cut, 3*flanks, by=cut)

numGenes = length(unique(data$gene))
genes = unique(data$gene)
cols = rainbow(length(genes))
#newGenes = unique(data2$gene)

xats = c(seq(-1*flanks, -1, 100), 0 , 2*flanks, seq(2*flanks+100, 3*flanks, 100))#append(seq(1,(flanks/cut)+1,100/cut), seq(2*flanks/cut,3*flanks/cut,100/cut))#c(1,2,3,4,5,6,10,11,12,13,14,15)
xlabs = c(seq(-1*flanks,-1, 100), "Start", "End", seq(100, flanks, 100))#append(append(append(seq(-1*flanks, -100, 100), c("Start")), c("End")), seq(100, flanks, 100))#c(-500,-400,-300,-200,-100,"Start","End", 100,200,300,400,500)
boxleft = 0
boxright = 2*flanks
cuts = seq(-1*flanks, 3*flanks, cut)
data.class = cut(data$mid, cuts)

#all values, used to calculate plot area and box size
mean.pi = tapply(data$pi, data.class, mean)
mean.pi = mean.pi[!is.na(mean.pi)]
mean.tajd = tapply(data$tajd, data.class, mean)
mean.tajd = mean.tajd[!is.na(mean.tajd)]

data2.class = cut(data2$mid, cuts)
mean2.pi = tapply(data2$pi, data2.class, mean)
mean2.pi = mean2.pi[!is.na(mean2.pi)]
mean2.tajd = tapply(data2$tajd, data2.class, mean)
mean2.tajd = mean2.tajd[!is.na(mean2.tajd)]

#plot size and box size
ymax = max(data$pi)
box_width = ymax*.05
ymin = -1*box_width
yrange = c(ymin, ymax)


if (pdfer == 1) {pdf(paste(path, 'pi_genes_ind_imp.pdf',sep=""), width = 12, height = 8)}

#set up the plot area, but dont put any data on it
plot(data$mid, data$pi, col = cols, type = "n", ylim = yrange, main = "Pi", xaxt='n', xlab = "Position (bp relative to gene ends)", ylab = "pi", bty='n', cex.axis = 1.5, cex.lab=1.5)

rect(boxleft,-1*box_width, boxright, 0, col = boxcolor)

#plot each imp gene
for (i in c(1:numGenes)){
    sdata = data[data$gene==genes[i],]
    points(sdata$mid, sdata$pi, col=cols[i], type="l")
}

#plot the "gene box" on the area
lines(meanx, mean2.pi, col=allcol, type="l", lwd=2, lty = 2)


#add in the x-axis labels
axis(1, at = xats, lab = xlabs)
#make the legend
names = c("all genes", as.vector(genes), "Gene Position")
legend('topright', names, col = c(allcol, cols, boxcolor), bty='n', lty = c(2, rep(1, times=numGenes), 1), lwd = c(2, rep(1, times=numGenes), 10), ncol = 2, cex = 1.5)
if (pdfer == 1) {dev.off()}


ymax = max(data$tajd)
ymin = min(data$tajd)#-3.0 #min(min(mean.tajd), min(mean2.tajd))
box_width = (ymax-ymin)*.05
ymin = ymin-1*box_width/2
yrange = c(ymin, ymax) 

if (pdfer == 1) {pdf(paste(path, 'tajd_genes_ind_imp.pdf',sep=""), width = 12, height = 8)} else {quartz()}
#set up the plot area, but dont put any data on it
plot(data$mid, data$tajd, col = cols, type = "n", ylim = yrange, main = "Taj D", xaxt='n', xlab = "Position (bp relative to gene ends)", ylab = "Taj D", bty='n', cex.axis = 1.5, cex.lab=1.5)

#plot the "gene box" on the area
rect(boxleft,ymin-1*box_width/2, boxright, ymin+box_width/2, col = boxcolor)

#plot each imp gene
for (i in c(1:numGenes)){
    sdata = data[data$gene==genes[i],]
    points(sdata$mid, sdata$tajd, col=cols[i], type="l")
}

lines(meanx, mean2.tajd, col=allcol, type="l", lty = 2, lwd=2)
#plot means from file 1
#plotMeans(data2$tajd, data.class, cols, 2)

#for file 2 plot means from a random subset of genes = the number from file 1
#for (i in c(1:5)){
#    subGenes = pickRandomGenes(numGenes, newGenes)
#    # print(subGenes)
#    subdata = data2[data2$gene %in% subGenes,]
#    
#    subdata.class = cut(subdata$mid, cuts)
#    plotMeans(subdata$tajd, subdata.class, cols2)
#}

#add in the x-axis labels
axis(1, at = xats, lab = xlabs)
#make the legend
#legend('topright', names, col = c(allcol, cols, boxcolor), bty='n', lty = c(2, rep(1, times=numGenes), 1), lwd = c(2, rep(1, times=numGenes), 10), ncol = 2, cex = 1.5)
if (pdfer==1){dev.off()}
