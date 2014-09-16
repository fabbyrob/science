data = read.table("/Users/wiliarj/Desktop/temp/diversity.200kb.f5.n1.jan6.txt", header=T)
names(data) = c("CHROM", "MIDPOINT", "SAMPLE", "HOMO", "HET", "ALL", "F")
regions = read.table("/Users/wiliarj/Desktop/temp/ranges.200kb.f5.n1.jan6.txt", header=T)
names(regions) = c("CHROM", "START", "END", "SAMPLE")

#filter <30% coverage sites
#data$HET = ifelse(data$ALL/20000 < 0.3, NA, data$HET)

colors = c("blue", "grey","blue", "grey","blue", "grey","blue", "grey")#colors = rainbow(8)
offset = 20000000
data$color = ifelse(data$CHROM=="scaffold_1", colors[1], ifelse(data$CHROM=="scaffold_2", colors[2], ifelse(data$CHROM=="scaffold_3", colors[3], ifelse(data$CHROM=="scaffold_4", colors[4], ifelse(data$CHROM=="scaffold_5", colors[5], ifelse(data$CHROM=="scaffold_6", colors[6], ifelse(data$CHROM=="scaffold_7", colors[7], colors[8])))))))
data$midpoint2 = ifelse(data$CHROM=="scaffold_1", data$MIDPOINT, ifelse(data$CHROM=="scaffold_2", data$MIDPOINT+offset, ifelse(data$CHROM=="scaffold_3", data$MIDPOINT+offset*2, ifelse(data$CHROM=="scaffold_4", data$MIDPOINT+offset*3, ifelse(data$CHROM=="scaffold_5", data$MIDPOINT+offset*4, ifelse(data$CHROM=="scaffold_6", data$MIDPOINT+offset*5, ifelse(data$CHROM=="scaffold_7", data$MIDPOINT+offset*6, data$MIDPOINT+offset*7)))))))

regions$start2 = ifelse(regions$CHROM=="scaffold_1", regions$START, ifelse(regions$CHROM=="scaffold_2", regions$START+offset, ifelse(regions$CHROM=="scaffold_3", regions$START+offset*2, ifelse(regions$CHROM=="scaffold_4", regions$START+offset*3, ifelse(regions$CHROM=="scaffold_5", regions$START+offset*4, ifelse(regions$CHROM=="scaffold_6", regions$START+offset*5, ifelse(regions$CHROM=="scaffold_7", regions$START+offset*6, regions$START+offset*7)))))))
regions$end2 = ifelse(regions$CHROM=="scaffold_1", regions$END, ifelse(regions$CHROM=="scaffold_2", regions$END+offset, ifelse(regions$CHROM=="scaffold_3", regions$END+offset*2, ifelse(regions$CHROM=="scaffold_4", regions$END+offset*3, ifelse(regions$CHROM=="scaffold_5", regions$END+offset*4, ifelse(regions$CHROM=="scaffold_6", regions$END+offset*5, ifelse(regions$CHROM=="scaffold_7", regions$END+offset*6, regions$END+offset*7)))))))
regions$color = ifelse(regions$CHROM=="scaffold_1", colors[1], ifelse(regions$CHROM=="scaffold_2", colors[2], ifelse(regions$CHROM=="scaffold_3", colors[3], ifelse(regions$CHROM=="scaffold_4", colors[4], ifelse(regions$CHROM=="scaffold_5", colors[5], ifelse(regions$CHROM=="scaffold_6", colors[6], ifelse(regions$CHROM=="scaffold_7", colors[7], colors[8])))))))

samps = unique(data$SAMPLE)
height = 5*length(samps)

pdf("/Users/wiliarj/Desktop/temp/test.pdf", height = height, width = 25)
par(mfrow=c(13,1))
samp = samps[1]
for (samp in samps){
  subData = data[data$SAMPLE == samp,]
  subRegions = regions[regions$SAMPLE == samp,]
  
  ylims = c(0,0.7)#
  xlims = c(min(subData$midpoint2), max(subData$midpoint2))
  
  for (scaf in unique(subData$CHROM)){
    subScaf = subData[subData$CHROM == scaf,]
    
    if (subScaf$CHROM[1] == "scaffold_1"){
      plot(subScaf$midpoint2, subScaf$HET/subScaf$HOMO, main = samp, type="l", lwd=3, col=subScaf$color, ylim=ylims, xlim=xlims, xaxt = "n", xlab="", ylab="# Het sites / # Homo sites")
      #points(subScaf$midpoint2, subScaf$F, main = samp, lty=2, type="l", col=subScaf$color, ylim=ylims, xlim=xlims, xaxt = "n", xlab="", ylab="# Het sites / # Homo sites")
      #points(subScaf$midpoint2, subScaf$ALL/200000, type="l", col="grey")
      #abline(h=0.5)
    } else{
      points(subScaf$midpoint2, subScaf$HET/subScaf$HOMO, type="l", col=subScaf$color, lwd=3)
      #points(subScaf$midpoint2, subScaf$F, lty=2, type="l", col=subScaf$color)
      #points(subScaf$midpoint2, subScaf$ALL/200000, type="l", col="grey")
    }
  }
  
  if (length(subRegions$START) > 0){
    segments(subRegions$start2, 0, subRegions$end2, 0, col="black", lwd=3)
  }
}
dev.off()

# regions$ypos = 1
# regions$ypos = ifelse(regions$SAMPLE=="103.17", 2, regions$ypos)....ifelse..ifelse
# pdf("/Users/wiliarj/Desktop/temp/regions.pdf", width=25, height=5)
# plot(1,1, ylim=c(1,11), xlim = xlims, ylab="individual", type='n', axes=F, bty='n', xlab="")
# segments(regions$start2, regions$ypos, regions$end2, regions$ypos, col=regions$color, lwd=3)
# axis(2, c(1,2,3,4,5,6,7,8,9,10,11), labels=unique(regions$SAMPLE), las=1)
# dev.off()