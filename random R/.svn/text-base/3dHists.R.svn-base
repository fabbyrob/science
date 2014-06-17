pdf("/Users/wiliarj/Desktop/temp/polyCrCg_hists.pdf")
library(latticeExtra)

########################
mydata = read.table("/Users/wiliarj/Desktop/temp/4fold.aligned",header=T)
attach(mydata)

test = ifelse(Cg == 0 & Cr == 0, F, T)
freqs = Num_Sites/sum(Num_Sites[test])

look = list(z=20,x=-50)
myscales = list(arrows=F)
cloud(freqs[test]~Cg[test]*Cr[test], xlim=c(27,-1), ylim=c(25,-1), panel.3d.cloud=panel.3dbars, par.settings=list(box.3d = list(col="transparent")), col.facet="grey",screen=look, xlab = "C. grandiflora", ylab = "C. rubella", zlab = "Frequency", scales=myscales, main="4fold polymorphism")


test = ifelse(Cg == 0 | Cr == 0, F, T)
freqs = Num_Sites/sum(Num_Sites[test])

look = list(z=20,x=-50)
myscales = list(arrows=F)
cloud(freqs[test]~Cg[test]*Cr[test], xlim=c(27,0), ylim=c(25,-0), panel.3d.cloud=panel.3dbars, par.settings=list(box.3d = list(col="transparent")), col.facet="grey",screen=look, xlab = "C. grandiflora", ylab = "C. rubella", zlab = "Frequency", scales=myscales, main="4fold polymorphism, neither fixed")
 

########################
mydata = read.table("/Users/wiliarj/Desktop/temp/0fold.aligned",header=T)
attach(mydata)

test = ifelse(Cg == 0 & Cr == 0, F, T)
freqs = Num_Sites/sum(Num_Sites[test])

look = list(z=20,x=-50)
myscales = list(arrows=F)
cloud(freqs[test]~Cg[test]*Cr[test], xlim=c(27,-1), ylim=c(25,-1), panel.3d.cloud=panel.3dbars, par.settings=list(box.3d = list(col="transparent")), col.facet="grey",screen=look, xlab = "C. grandiflora", ylab = "C. rubella", zlab = "Frequency", scales=myscales, main="0fold polymorphism")


test = ifelse(Cg == 0 | Cr == 0, F, T)
freqs = Num_Sites/sum(Num_Sites[test])

look = list(z=20,x=-50)
myscales = list(arrows=F)
cloud(freqs[test]~Cg[test]*Cr[test], xlim=c(27,0), ylim=c(25,-0), panel.3d.cloud=panel.3dbars, par.settings=list(box.3d = list(col="transparent")), col.facet="grey",screen=look, xlab = "C. grandiflora", ylab = "C. rubella", zlab = "Frequency", scales=myscales, main="0fold polymorphism, neither fixed")
 
########################
mydata = read.table("/Users/wiliarj/Desktop/temp/intron.aligned",header=T)
attach(mydata)

test = ifelse(Cg == 0 & Cr == 0, F, T)
freqs = Num_Sites/sum(Num_Sites[test])

look = list(z=20,x=-50)
myscales = list(arrows=F)
cloud(freqs[test]~Cg[test]*Cr[test], xlim=c(27,-1), ylim=c(25,-1), panel.3d.cloud=panel.3dbars, par.settings=list(box.3d = list(col="transparent")), col.facet="grey",screen=look, xlab = "C. grandiflora", ylab = "C. rubella", zlab = "Frequency", scales=myscales, main="intron polymorphism")


test = ifelse(Cg == 0 | Cr == 0, F, T)
freqs = Num_Sites/sum(Num_Sites[test])

look = list(z=20,x=-50)
myscales = list(arrows=F)
cloud(freqs[test]~Cg[test]*Cr[test], xlim=c(27,0), ylim=c(25,-0), panel.3d.cloud=panel.3dbars, par.settings=list(box.3d = list(col="transparent")), col.facet="grey",screen=look, xlab = "C. grandiflora", ylab = "C. rubella", zlab = "Frequency", scales=myscales, main="intron polymorphism, neither fixed")
 
########################
mydata = read.table("/Users/wiliarj/Desktop/temp/intergene.aligned",header=T)
attach(mydata)

test = ifelse(Cg == 0 & Cr == 0, F, T)
freqs = Num_Sites/sum(Num_Sites[test])

look = list(z=20,x=-50)
myscales = list(arrows=F)
cloud(freqs[test]~Cg[test]*Cr[test], xlim=c(27,-1), ylim=c(25,-1), panel.3d.cloud=panel.3dbars, par.settings=list(box.3d = list(col="transparent")), col.facet="grey",screen=look, xlab = "C. grandiflora", ylab = "C. rubella", zlab = "Frequency", scales=myscales, main="intergene polymorphism")


test = ifelse(Cg == 0 | Cr == 0, F, T)
freqs = Num_Sites/sum(Num_Sites[test])

look = list(z=20,x=-50)
myscales = list(arrows=F)
cloud(freqs[test]~Cg[test]*Cr[test], xlim=c(27,0), ylim=c(25,-0), panel.3d.cloud=panel.3dbars, par.settings=list(box.3d = list(col="transparent")), col.facet="grey",screen=look, xlab = "C. grandiflora", ylab = "C. rubella", zlab = "Frequency", scales=myscales, main="intergene polymorphism, neither fixed")
 
########################
mydata = read.table("/Users/wiliarj/Desktop/temp/3utr.aligned",header=T)
attach(mydata)

test = ifelse(Cg == 0 & Cr == 0, F, T)
freqs = Num_Sites/sum(Num_Sites[test])

look = list(z=20,x=-50)
myscales = list(arrows=F)
cloud(freqs[test]~Cg[test]*Cr[test], xlim=c(27,-1), ylim=c(25,-1), panel.3d.cloud=panel.3dbars, par.settings=list(box.3d = list(col="transparent")), col.facet="grey",screen=look, xlab = "C. grandiflora", ylab = "C. rubella", zlab = "Frequency", scales=myscales, main="3utr polymorphism")


test = ifelse(Cg == 0 | Cr == 0, F, T)
freqs = Num_Sites/sum(Num_Sites[test])

look = list(z=20,x=-50)
myscales = list(arrows=F)
cloud(freqs[test]~Cg[test]*Cr[test], xlim=c(27,0), ylim=c(25,-0), panel.3d.cloud=panel.3dbars, par.settings=list(box.3d = list(col="transparent")), col.facet="grey",screen=look, xlab = "C. grandiflora", ylab = "C. rubella", zlab = "Frequency", scales=myscales, main="3utr polymorphism, neither fixed")
 
########################
mydata = read.table("/Users/wiliarj/Desktop/temp/5utr.aligned",header=T)
attach(mydata)

test = ifelse(Cg == 0 & Cr == 0, F, T)
freqs = Num_Sites/sum(Num_Sites[test])

look = list(z=20,x=-50)
myscales = list(arrows=F)
cloud(freqs[test]~Cg[test]*Cr[test], xlim=c(27,-1), ylim=c(25,-1), panel.3d.cloud=panel.3dbars, par.settings=list(box.3d = list(col="transparent")), col.facet="grey",screen=look, xlab = "C. grandiflora", ylab = "C. rubella", zlab = "Frequency", scales=myscales, main="5utr polymorphism")


test = ifelse(Cg == 0 | Cr == 0, F, T)
freqs = Num_Sites/sum(Num_Sites[test])

look = list(z=20,x=-50)
myscales = list(arrows=F)
cloud(freqs[test]~Cg[test]*Cr[test], xlim=c(27,0), ylim=c(25,-0), panel.3d.cloud=panel.3dbars, par.settings=list(box.3d = list(col="transparent")), col.facet="grey",screen=look, xlab = "C. grandiflora", ylab = "C. rubella", zlab = "Frequency", scales=myscales, main="5utr polymorphism, neither fixed")
 
dev.off()
