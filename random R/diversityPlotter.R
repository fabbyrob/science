data = read.table("/Users/wiliarj/Desktop/temp/4fold.diversity.5ksnps.txt", header=T)

offset = 20000000

centromere = read.csv("/Users/wiliarj/Documents/data/grandiflora\ noncoding/centromere.txt", header = F)
centromere$V2 = ifelse(centromere$V1=="1", centromere$V2, ifelse(centromere$V1=="2", centromere$V2+offset, ifelse(centromere$V1=="3", centromere$V2+offset*2, ifelse(centromere$V1=="4", centromere$V2+offset*3, ifelse(centromere$V1=="5", centromere$V2+offset*4, ifelse(centromere$V1=="6", centromere$V2+offset*5, ifelse(centromere$V1=="7", centromere$V2+offset*6, centromere$V2+offset*7)))))))
centromere$V3 = ifelse(centromere$V1=="1", centromere$V3, ifelse(centromere$V1=="2", centromere$V3+offset, ifelse(centromere$V1=="3", centromere$V3+offset*2, ifelse(centromere$V1=="4", centromere$V3+offset*3, ifelse(centromere$V1=="5", centromere$V3+offset*4, ifelse(centromere$V1=="6", centromere$V3+offset*5, ifelse(centromere$V1=="7", centromere$V3+offset*6, centromere$V3+offset*7)))))))

colors = c("blue", "grey","blue", "grey","blue", "grey","blue", "grey")# rainbow(8)

data$color = ifelse(data$chromosome=="scaffold_1", colors[1], ifelse(data$chromosome=="scaffold_2", colors[2], ifelse(data$chromosome=="scaffold_3", colors[3], ifelse(data$chromosome=="scaffold_4", colors[4], ifelse(data$chromosome=="scaffold_5", colors[5], ifelse(data$chromosome=="scaffold_6", colors[6], ifelse(data$chromosome=="scaffold_7", colors[7], colors[8])))))))
data$midpoint2 = ifelse(data$chromosome=="scaffold_1", data$midpoint, ifelse(data$chromosome=="scaffold_2", data$midpoint+offset, ifelse(data$chromosome=="scaffold_3", data$midpoint+offset*2, ifelse(data$chromosome=="scaffold_4", data$midpoint+offset*3, ifelse(data$chromosome=="scaffold_5", data$midpoint+offset*4, ifelse(data$chromosome=="scaffold_6", data$midpoint+offset*5, ifelse(data$chromosome=="scaffold_7", data$midpoint+offset*6, data$midpoint+offset*7)))))))

pdf("/Users/wiliarj/Desktop/temp/div.flt.pdf", height = 8, width = 25)

par(mfrow=c(2,1))

ylims = c(min(data$pi), max(data$pi))
xlims = c(min(data$midpoint2)*0.9, max(data$midpoint2))

for (scaf in unique(data$chromosome)){
  subScaf = data[data$chromosome == scaf,]
  
  if (subScaf$chromosome[1] == "scaffold_1"){
    plot(subScaf$midpoint2, subScaf$pi, main = "Pi", type="l", lwd=2, col=subScaf$color, ylim=ylims, xlim=xlims, xaxt = "n", xlab="", ylab="Pi")
  } else{
    points(subScaf$midpoint2, subScaf$pi, type="l", lwd=2, col=subScaf$color)
  }
}

#centromeres
rect(centromere$V2, ylims[1], centromere$V3, ylims[1]*1.05, col = "grey")


#ylims = c(min(data$tajd, na.rm=T), max(data$tajd, na.rm=T))
#xlims = c(min(subData$midpoint2), max(subData$midpoint2))

#for (scaf in unique(data$chromosome)){
#  subScaf = data[data$chromosome == scaf,]
#  
#  if (subScaf$chromosome[1] == "scaffold_1"){
#    plot(subScaf$midpoint2, subScaf$tajd, main = "TajD", type="l", col=subScaf$color, ylim=ylims, xlim=xlims, xaxt = "n", xlab="", ylab="TajD")
#  } else{
#    points(subScaf$midpoint2, subScaf$tajd, type="l", col=subScaf$color)
#  }
#}

ylims = c(min(data$divergence, na.rm=T), max(data$divergence, na.rm=T))
xlims = c(min(data$midpoint2), max(data$midpoint2))

for (scaf in unique(data$chromosome)){
  subScaf = data[data$chromosome == scaf,]
  
  if (subScaf$chromosome[1] == "scaffold_1"){
    plot(subScaf$midpoint2, subScaf$divergence, main = "Divergence", type="l", lwd=2, col=subScaf$color, ylim=ylims, xlim=xlims, xaxt = "n", xlab="", ylab="Divergence")
  } else{
    points(subScaf$midpoint2, subScaf$divergence, type="l", col=subScaf$color, lwd=2)
  }
}

#centromeres
rect(centromere$V2, ylims[1], centromere$V3, ylims[1]*1.01, col = "grey")
dev.off()