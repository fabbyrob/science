data = read.csv('/Users/wiliarj/Desktop/temp/scaf1.out.csv',header=T)

names(data) = c("gene", "mid", "pi", "theta", "tajd")

flanks = 500
genes = unique(data$gene)
cols = "red"#rainbow(length(genes))
first = 1
xrange = c(-1*flanks, 2*flanks)

boxcolor = "grey"


xats = c(-500,-400,-300,-200,-100,0,500, 600, 700, 800,900, 1000)
xlabs = c(-500,-400,-300,-200,-100,"Start","End", 100,200,300,400,500)

ymax = max(data$pi)
box_width = ymax*.05
ymin = -1*box_width/2
yrange = c(ymin, ymax)

pdf('/Users/wiliarj/Desktop/temp/pi_genes.pdf', width = 10, height = 7)
c = 1
for (gene in genes){
    geneData = data[data$gene==gene,]
    if (first == 1){
        plot(geneData$mid, geneData$pi, col = cols, type = "l", xlim = xrange, ylim = yrange, main = "Pi", xaxt='n', xlab = "Position (bp relative to gene ends)", ylab = "pi")
        axis(1, at = xats, lab = xlabs)
        rect(0,-1*box_width/2, flanks, box_width/2, col = boxcolor)
        abline(h=0)
    }else{
        points(geneData$mid, geneData$pi, col = cols, type = "l")
    }
    first = 0
    c = c + 1
}
legend('topleft', c("pi values", "Gene Position"), col = c(cols, boxcolor), bty='n', lty=c(1,1), lwd = c(1, 10))
dev.off()

ymax = max(data$tajd)
ymin = min(data$tajd)
box_width = (ymax-ymin)*.05
yrange = c(ymin, ymax) 
pdf('/Users/wiliarj/Desktop/temp/tajd_genes.pdf', width = 10, height = 7)
c = 1
first = 1
for (gene in genes){
    geneData = data[data$gene==gene,]
    if (first == 1){
        plot(geneData$mid, geneData$tajd, col = cols, type = "l", xlim = xrange, ylim = yrange, main = "Taj D", xaxt='n', xlab = "Position (bp relative to gene ends)", ylab = "Taj D")
        axis(1, at = xats, lab = xlabs)
        rect(0,-1*box_width/2, flanks, box_width/2, col = boxcolor)
        abline(h=0)
    }else{
        points(geneData$mid, geneData$tajd, col = cols, type = "l")
    }
    first = 0
    c = c + 1
}
legend('topleft', c("Taj D values", "Gene Position"), col = c(cols, boxcolor), bty='n', lty=c(1,1), lwd = c(1, 10))
dev.off()