pollen = read.table("/Users/wiliarj/Desktop/temp/pollen.full.counts", header=F)
sporo = read.table("/Users/wiliarj/Desktop/temp/sporo.full.counts", header=F)
low = read.table("/Users/wiliarj/Desktop/temp/low.exp.full.counts", header=F)
high = read.table("/Users/wiliarj/Desktop/temp/high.exp.full.counts", header=F)
all = read.table("/Users/wiliarj/Desktop/temp/full.cnc.counts", header=F)

fpkm = read.table("/Users/wiliarj/Desktop/temp/genes.fpkm.parse", header=F)

pollen$type = "pollen"
sporo$type = "sporophyte"
low$type = "low expression"
high$type = "high expression"
all$type = "all"

full = data.frame(cnc=c(pollen[,2], sporo[,2], low[,2], high[,2], all[,2]), type = c(pollen$type, sporo$type, low$type, high$type, all$type))
boxplot(full$cnc ~ full$type, ylab = "CNCs per gene (bp)", ylim=c(0,300))

low_exp = merge(low, fpkm, by = "V1")
high_exp = merge(high, fpkm, by = "V1")

all_exp = merge(all, fpkm, by = "V1")#data.frame(gene=c(low_exp$V1, high_exp$V1), exp=c(low_exp$V2.y, high_exp$V2.y), cnc=c(low_exp$V2.x, high_exp$V2.x))

plot(log(low_exp$V2.y), low_exp$V2.x, xlab="expression", ylab="CNCs")
llm = lm(log(low_exp$V2.y) ~ low_exp$V2.x)
abline(llm, lwd=3, col="red")

plot(log(high_exp$V2.y), high_exp$V2.x, xlab="expression", ylab="CNCs")
hlm = lm(log(high_exp$V2.y) ~ high_exp$V2.x)
abline(hlm, lwd=3, col="red")

plot(log(all_exp$V2.y), all_exp$V2.x, xlab="log(expression level)", ylab="CNCs per gene (bp)", ylim=c(0, 3000))
plot(all_exp$V3, all_exp$V2.x, xlab="coeff of variation in expression", ylab="CNCs per gene (bp)", ylim=c(0, 3000))
