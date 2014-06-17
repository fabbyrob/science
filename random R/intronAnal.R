data = read.table("/Users/wiliarj/Desktop/temp/full.div.len.txt", header=T)

usedCutOff = 1#at least 50% of the length of the intron must have been analyzable
divCutOff = 1

data$Divergence = data$Diverged/data$SitesTot
data$PropLengthUsed = data$SitesTot/data$IntronLen

subData = data[data$PropLengthUsed >= usedCutOff & data$Divergence < divCutOff,]
sitesLost = length(data$Gene)-length(subData$Gene)#number of sites that didn't pass filter
meanLen = mean(subData$IntronLen)
medianLen = median(subData$IntronLen)


onlyFirst = subData[subData$IntronNum==0,]
nonFirst = subData[subData$IntronNum!=0,]

onlyLong = subData[subData$IntronLen>=medianLen,]
onlyShort = subData[subData$IntronLen<medianLen,]

#### Plot divergence vs. Length
quartz()
plot(subData$IntronLen, subData$Divergence, xlab="Intron length (bp)", ylab="Divergence", main="Intron Length vs. Divergence")
mycor = cor(subData$IntronLen, subData$Divergence, method="spearman")
mycor.test = cor.test(subData$IntronLen, subData$Divergence, method="spearman")$p.value
text(max(subData$IntronLen)-.2*max(subData$IntronLen), max(subData$Divergence)-.1, paste("Rs =", format(mycor, digits=2), "\np =",format(mycor.test, digits=2)))

#quartz()
#plot(onlyFirst$IntronLen, onlyFirst$Divergence, xlab="Intron length (bp)", ylab="Divergence", main="Intron Length vs. Divergence First introns")
#mycor = cor(onlyFirst$IntronLen, onlyFirst$Divergence, method="spearman")
#mycor.test = cor.test(onlyFirst$IntronLen, onlyFirst$Divergence, method="spearman")$p.value
#text(max(onlyFirst$IntronLen)-.2*max(onlyFirst$IntronLen), max(onlyFirst$Divergence)-.1, paste("Rs =", format(mycor, digits=2), "\np =",format(mycor.test, digits=2)))
#
#quartz()
#plot(nonFirst$IntronLen, nonFirst$Divergence, xlab="Intron length (bp)", ylab="Divergence", main="Intron Length vs. Divergence Non-first introns")
#mycor = cor(nonFirst$IntronLen, nonFirst$Divergence, method="spearman")
#mycor.test = cor.test(nonFirst$IntronLen, nonFirst$Divergence, method="spearman")$p.value
#text(max(nonFirst$IntronLen)-.2*max(nonFirst$IntronLen), max(subData$Divergence)-.1, paste("Rs =", format(mycor, digits=2), "\np =",format(mycor.test, digits=2)))
#
#quartz()
#plot(onlyLong$IntronLen, onlyLong$Divergence, xlab="Intron length (bp)", ylab="Divergence", main="Intron Length vs. Divergence Long introns")
#mycor = cor(onlyLong$IntronLen, onlyLong$Divergence, method="spearman")
#mycor.test = cor.test(onlyLong$IntronLen, onlyFirst$Divergence, method="spearman")$p.value
#text(max(onlyLong$IntronLen)-.2*max(onlyLong$IntronLen), max(onlyLong$Divergence)-.1, paste("Rs =", format(mycor, digits=2), "\np =",format(mycor.test, digits=2)))
#
#quartz()
#plot(onlyShort$IntronLen, onlyShort$Divergence, xlab="Intron length (bp)", ylab="Divergence", main="Intron Length vs. Divergence Short introns")
#mycor = cor(onlyShort$IntronLen, onlyShort$Divergence, method="spearman")
#mycor.test = cor.test(onlyShort$IntronLen, onlyFirst$Divergence, method="spearman")$p.value
#text(max(onlyShort$IntronLen)-.2*max(onlyShort$IntronLen), max(onlyShort$Divergence)-.1, paste("Rs =", format(mycor, digits=2), "\np =",format(mycor.test, digits=2)))


#### Plot divergence vs. intron order
quartz()
plot(subData$IntronNum, subData$Divergence, xlab="Intron order", ylab="Divergence", main="Intron Order vs. Divergence")
mycor = cor(subData$IntronNum, subData$Divergence, method="spearman")
mycor.test = cor.test(subData$IntronNum, subData$Divergence, method="spearman")$p.value
text(max(subData$IntronNum)-.2*max(subData$IntronNum), max(subData$Divergence)-.1, paste("Rs =", format(mycor, digits=2), "\np =",format(mycor.test, digits=2)))


###barplot of intron order/divergence
barplot()