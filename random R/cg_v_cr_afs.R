mydata = read.table("/Users/wiliarj/Documents/workspaces/general/robert.williamson/trunk/random R/io/aligned.afs.txt", header=T)
cloud(freqs ~ mydata$CgAF * mydata$CrAF, xlab = "Cg", ylab = "Cr", xlim=c(-1,2), ylim=c(-1,2), type="h",xbase=.7, ybase=.7, panel.3d.cloud = panel.3dbars, scales=list(arrows=F), screen = list(z=-60,x=-50), col.facet="grey")
