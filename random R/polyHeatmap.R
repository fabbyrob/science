pdf("/Users/wiliarj/Desktop/temp/polyCrCg_heatmap.pdf")
library(latticeExtra)
breaks1 = 50000
breaks2 = 50
########################
mydata = read.table("/Users/wiliarj/Desktop/temp/4fold.aligned.2.afs",header=T)
attach(mydata)

test = ifelse(Cg == 0 & Cr == 0, F, T)
freqs = Num_Sites/sum(Num_Sites[test])
col.l <- colorRampPalette(c('blue', 'green', 'purple', 'yellow', 'red'))(breaks1)

ats = c()
max = max(freqs[test])

levelplot(freqs[test]~Cg[test]*Cr[test], at=do.breaks(c(0,max), breaks1), col.regions=col.l, xlab = "C. grandiflora", ylab = "C. rubella", main = "4fold polymorphism", xlim=c(-1,27), ylim=c(-1,25))

test = ifelse(Cg == 0 | Cr == 0, F, T)
freqs = Num_Sites/sum(Num_Sites[test])
col.l <- colorRampPalette(c('blue', 'green', 'purple', 'yellow', 'red'))(breaks2)

ats = c()
max = max(freqs[test])

levelplot(freqs[test]~Cg[test]*Cr[test], at=do.breaks(c(0,max), breaks2), col.regions=col.l, xlab = "C. grandiflora", ylab = "C. rubella", main = "4fold polymorphism, neither fixed", xlim=c(0,27), ylim=c(0,25))

########################
mydata = read.table("/Users/wiliarj/Desktop/temp/0fold.aligned.2.afs",header=T)
attach(mydata)

test = ifelse(Cg == 0 & Cr == 0, F, T)
freqs = Num_Sites/sum(Num_Sites[test])
col.l <- colorRampPalette(c('blue', 'green', 'purple', 'yellow', 'red'))(breaks1)

ats = c()
max = max(freqs[test])

levelplot(freqs[test]~Cg[test]*Cr[test], at=do.breaks(c(0,max), breaks1), col.regions=col.l, xlab = "C. grandiflora", ylab = "C. rubella", main = "0fold polymorphism", xlim=c(-1,27), ylim=c(-1,25))

test = ifelse(Cg == 0 | Cr == 0, F, T)
freqs = Num_Sites/sum(Num_Sites[test])
col.l <- colorRampPalette(c('blue', 'green', 'purple', 'yellow', 'red'))(breaks2)

ats = c()
max = max(freqs[test])

levelplot(freqs[test]~Cg[test]*Cr[test], at=do.breaks(c(0,max), breaks2), col.regions=col.l, xlab = "C. grandiflora", ylab = "C. rubella", main = "0fold polymorphism, neither fixed", xlim=c(0,27), ylim=c(0,25))

########################
mydata = read.table("/Users/wiliarj/Desktop/temp/intron.aligned.2.afs",header=T)
attach(mydata)

test = ifelse(Cg == 0 & Cr == 0, F, T)
freqs = Num_Sites/sum(Num_Sites[test])
col.l <- colorRampPalette(c('blue', 'green', 'purple', 'yellow', 'red'))(breaks1)

ats = c()
max = max(freqs[test])

levelplot(freqs[test]~Cg[test]*Cr[test], at=do.breaks(c(0,max), breaks1), col.regions=col.l, xlab = "C. grandiflora", ylab = "C. rubella", main = "intron polymorphism", xlim=c(-1,27), ylim=c(-1,25))

test = ifelse(Cg == 0 | Cr == 0, F, T)
freqs = Num_Sites/sum(Num_Sites[test])
col.l <- colorRampPalette(c('blue', 'green', 'purple', 'yellow', 'red'))(breaks2)

ats = c()
max = max(freqs[test])

levelplot(freqs[test]~Cg[test]*Cr[test], at=do.breaks(c(0,max), breaks2), col.regions=col.l, xlab = "C. grandiflora", ylab = "C. rubella", main = "intron polymorphism, neither fixed", xlim=c(0,27), ylim=c(0,25))

########################
mydata = read.table("/Users/wiliarj/Desktop/temp/intergene.aligned.2.afs",header=T)
attach(mydata)

test = ifelse(Cg == 0 & Cr == 0, F, T)
freqs = Num_Sites/sum(Num_Sites[test])
col.l <- colorRampPalette(c('blue', 'green', 'purple', 'yellow', 'red'))(breaks1)

ats = c()
max = max(freqs[test])

levelplot(freqs[test]~Cg[test]*Cr[test], at=do.breaks(c(0,max), breaks1), col.regions=col.l, xlab = "C. grandiflora", ylab = "C. rubella", main = "intergene polymorphism", xlim=c(-1,27), ylim=c(-1,25))

test = ifelse(Cg == 0 | Cr == 0, F, T)
freqs = Num_Sites/sum(Num_Sites[test])
col.l <- colorRampPalette(c('blue', 'green', 'purple', 'yellow', 'red'))(breaks2)

ats = c()
max = max(freqs[test])

levelplot(freqs[test]~Cg[test]*Cr[test], at=do.breaks(c(0,max), breaks2), col.regions=col.l, xlab = "C. grandiflora", ylab = "C. rubella", main = "intergene polymorphism, neither fixed", xlim=c(0,27), ylim=c(0,25))

########################
mydata = read.table("/Users/wiliarj/Desktop/temp/3utr.aligned.2.afs",header=T)
attach(mydata)

test = ifelse(Cg == 0 & Cr == 0, F, T)
freqs = Num_Sites/sum(Num_Sites[test])
col.l <- colorRampPalette(c('blue', 'green', 'purple', 'yellow', 'red'))(breaks1)

ats = c()
max = max(freqs[test])

levelplot(freqs[test]~Cg[test]*Cr[test], at=do.breaks(c(0,max), breaks1), col.regions=col.l, xlab = "C. grandiflora", ylab = "C. rubella", main = "3utr polymorphism", xlim=c(-1,27), ylim=c(-1,25))

test = ifelse(Cg == 0 | Cr == 0, F, T)
freqs = Num_Sites/sum(Num_Sites[test])
col.l <- colorRampPalette(c('blue', 'green', 'purple', 'yellow', 'red'))(breaks2)

ats = c()
max = max(freqs[test])

levelplot(freqs[test]~Cg[test]*Cr[test], at=do.breaks(c(0,max), breaks2), col.regions=col.l, xlab = "C. grandiflora", ylab = "C. rubella", main = "3utr polymorphism, neither fixed", xlim=c(0,27), ylim=c(0,25))

########################
mydata = read.table("/Users/wiliarj/Desktop/temp/5utr.aligned.2.afs",header=T)
attach(mydata)

test = ifelse(Cg == 0 & Cr == 0, F, T)
freqs = Num_Sites/sum(Num_Sites[test])
col.l <- colorRampPalette(c('blue', 'green', 'purple', 'yellow', 'red'))(breaks1)

ats = c()
max = max(freqs[test])

levelplot(freqs[test]~Cg[test]*Cr[test], at=do.breaks(c(0,max), breaks1), col.regions=col.l, xlab = "C. grandiflora", ylab = "C. rubella", main = "5utr polymorphism", xlim=c(-1,27), ylim=c(-1,25))

test = ifelse(Cg == 0 | Cr == 0, F, T)
freqs = Num_Sites/sum(Num_Sites[test])
col.l <- colorRampPalette(c('blue', 'green', 'purple', 'yellow', 'red'))(breaks2)

ats = c()

levelplot(freqs[test]~Cg[test]*Cr[test], at=do.breaks(c(0,max), breaks2), col.regions=col.l, xlab = "C. grandiflora", ylab = "C. rubella", main = "5utr polymorphism, neither fixed", xlim=c(0,27), ylim=c(0,25))

dev.off()
