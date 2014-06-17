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

#path = '/Users/wiliarj/Desktop/temp/'
path = 'C:/Users/Fabby/Desktop/temp/'

data = read.csv(paste(path, 'total.imp.out.csv',sep=""),header=T)
data2 = read.csv(paste(path, 'total.csv',sep=""),header=T)
names(data) = c("gene", "mid", "pi", "theta", "tajd")
names(data2) = c("gene", "mid", "pi", "theta", "tajd")
flanks = 500
genes = unique(data$gene)
allcol = "green"
cols = "red"#rainbow(length(genes))
cols2 = "blue"
first = 1
xrange = c(-1*flanks, 2*flanks)

boxcolor = "grey"
cut = 50

numGenes = length(unique(data$gene))
newGenes = unique(data2$gene)

xats = append(seq(1,(flanks/cut)+1,100/cut), seq(2*flanks/cut,3*flanks/cut,100/cut))#c(1,2,3,4,5,6,10,11,12,13,14,15)
xlabs = append(append(append(seq(-1*flanks, -100, 100), c("Start")), c("End")), seq(100, flanks, 100))#c(-500,-400,-300,-200,-100,"Start","End", 100,200,300,400,500)
boxleft = flanks/cut+1
boxright = 2*flanks/cut
cuts = seq(-1*flanks, 2*flanks, cut)
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
ymax = .15#max(max(mean.pi), max(mean2.pi))
box_width = ymax*.05
ymin = -1*box_width/2
yrange = c(ymin, ymax)

pdf(paste(path, 'pi_genes_mean.pdf',sep=""), width = 10, height = 7)

#set up the plot area, but dont put any data on it
plot(mean.pi, col = cols, type = "n", ylim = yrange, main = "Pi", xaxt='n', xlab = "Position (bp relative to gene ends)", ylab = "pi", bty='n')

#plot the "gene box" on the area
rect(boxleft,-1*box_width/2, boxright, box_width/2, col = boxcolor)
lines(mean2.pi, col=allcol, type="l", lwd=2)
#plot means from file 1
plotMeans(data$pi, data.class, cols, 2)

#for file 2 plot means from a random subset of genes = the number from file 1
for (i in c(1:5)){
    subGenes = pickRandomGenes(numGenes, newGenes)
    #print(subGenes)
    subdata = data2[data2$gene %in% subGenes,]
    
    subdata.class = cut(subdata$mid, cuts)
    plotMeans(subdata$pi, subdata.class, cols2)
}

#add in the x-axis labels
axis(1, at = xats, lab = xlabs)
#make the legend
legend('topright', c("all genes","imprinted genes", "random genes", "Gene Position"), col = c(allcol, cols, cols2, boxcolor), bty='n', lty=c(1,1,1,1), lwd = c(1,1,1, 10))
dev.off()

ymax = 1.5#max(max(mean.tajd), max(mean2.tajd))
ymin = -3.0 #min(min(mean.tajd), min(mean2.tajd))
box_width = (ymax-ymin)*.05
ymin = ymin-1*box_width/2
yrange = c(ymin, ymax) 

pdf(paste(path, 'tajd_genes_mean.pdf',sep=""), width = 10, height = 7)

#set up the plot area, but dont put any data on it
plot(mean.tajd, col = cols, type = "n", ylim = yrange, main = "Taj D", xaxt='n', xlab = "Position (bp relative to gene ends)", ylab = "Taj D", bty='n')

#plot the "gene box" on the area
rect(boxleft,ymin-1*box_width/2, boxright, ymin+box_width/2, col = boxcolor)
lines(mean2.tajd, col=allcol, type="l", lwd=2)
#plot means from file 1
plotMeans(data$tajd, data.class, cols, 2)

#for file 2 plot means from a random subset of genes = the number from file 1
for (i in c(1:5)){
    subGenes = pickRandomGenes(numGenes, newGenes)
   # print(subGenes)
    subdata = data2[data2$gene %in% subGenes,]
    
    subdata.class = cut(subdata$mid, cuts)
    plotMeans(subdata$tajd, subdata.class, cols2)
}

#add in the x-axis labels
axis(1, at = xats, lab = xlabs)
#make the legend
legend('bottomleft', c("all genes", "imprinted genes", "random genes", "Gene Position"), col = c(allcol, cols, cols2, boxcolor), bty='n', lty=c(1,1,1,1), lwd = c(1,1,1, 10))
dev.off()
