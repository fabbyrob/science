myNames = c("midpoint","pi","theta","tajd")
path = "/data/robert.williamson/finalGenome_analysis/windows/"
#path = '/Users/wiliarj/Desktop/temp/'

options(scipen=15)#make non-scientific notation
pdf(file=paste(path, 'pi_windows.pdf', sep=""),width=15, height=11)
par(mfrow=c(4,2), oma=c(5,0,5,0))

for (i in 1:8){
    scafdata = read.csv(paste(path, "scaf",i,".1010SNP1.csv", sep = ""), header=T)
    names(scafdata) = myNames
    
    plot(scafdata$midpoint, scafdata$pi, main=paste("chromosome ",i, sep = ""), type="l", ylab = "pi", xlab = "position", ylim=c(0,0.025))
    abline(mean(scafdata$pi),0,col="red", lty=2)
}

title(main="Pi", outer=T, cex.main=5, font.main=3)
dev.off()

#dev.new(width=12, height=12)
pdf(file=paste(path, 'tajd_windows.pdf', sep=""),width=15, height=11)
par(mfrow=c(4,2), oma=c(5,0,5,0))

for (i in 1:8){
    scafdata = read.csv(paste(path, "scaf",i,".1010SNP1.csv", sep = ""), header=T)
    names(scafdata) = myNames
    
    plot(scafdata$midpoint, as.vector(scafdata$tajd), main=paste("chromosome ",i), type="l", ylab = "Taj D", xlab = "position", ylim=c(-2,1))
    abline(mean(scafdata$tajd),0,col="red", lty=2)
    abline(0,0)
}

title(main="Tajima's D", outer=T, cex.main=5, font.main=3)
dev.off()

 ########## FUNCTIONS ###############
 plotGenes <- function(genes)
 {
     for (i in 1:length(genes[,1])){
        segments(x0=genes[i,2], y0=.005, x1=genes[i,3],y1=0.005, lwd = 3, col= "red")
     }
 }
 
 plotNC <- function(genes, end)
 {
     segments(x0=0, y0=.006, x1=genes[1,2],y1=0.006, lwd = 1, col= "blue")
     for (i in 2:length(genes[,1])-1){
        segments(x0=genes[i,2], y0=.006, x1=genes[i,3],y1=0.006, lwd = 1, col= "blue")
     }
     segments(x0=genes[i+1,3], y0=.006, x1=end,y1=0.006, lwd = 1, col= "blue")
 }
 
 getAxis <- function(max, tickIncrement = 5)
 {
     current = 0
     pos = c(0)
     labs = c(0)
     while(current < max){
         current = current+2.5*10^5
         lab = current/(10^6)
         if (current %% (tickIncrement*10^6) == 0){
             pos = append(pos, current)
             labs = append(labs,paste(lab,"Mb"))
         }
         if (current > max){
             if (current - pos[length(pos)] > 2*10^6){
                 pos = append(pos, current)
                 labs = append(labs,paste(lab,"Mb"))
             }
             else{
                 pos[length(pos)] = current
                 labs[length(labs)] = paste(lab,"Mb")
             }
         }
     }
     list(pos,labs)
 }
 ##################
 
 #options(scipen=15)#make non-scientific notation
#path = '/Users/wiliarj/Desktop/temp/'
path = "/data/robert.williamson/finalGenome_analysis/windows/"
 for (i in 1:8){
     myNames = c("midpoint","pi","theta","tajd")
     zfolddata = read.csv(paste(path, "scaf",i,".0fold.1010SNP1.csv",sep=""), header=T)
     names(zfolddata) = myNames
     ffolddata = read.csv(paste(path, "scaf",i,".4fold.1010SNP1.csv",sep=""), header=T)
     names(ffolddata) = myNames
     introndata = read.csv(paste(path, "scaf",i,".intron.1010SNP1.csv",sep=""), header=T)
     names(introndata) = myNames
     intergenedata = read.csv(paste(path, "scaf",i,".intergene.1010SNP1.csv",sep=""), header=T)
     names(intergenedata) = myNames
     tutrdata = read.csv(paste(path, "scaf",i,".3utr.1010SNP1.csv",sep=""), header=T)
     names(tutrdata) = myNames
     futrdata = read.csv(paste(path, "scaf",i,".5utr.1010SNP1.csv",sep=""), header=T)
     names(futrdata) = myNames
 
     #annot = read.table(paste("/data/robert.williamson/finalGenome_analysis/annot/scaf",i,".annot.sorted",sep=""),header=F)
 
     #----Pi plots----
     #dev.new(width=12, height=12)
     maxes=c(tail(ffolddata$midpoint,n=1),tail(zfolddata$midpoint,n=1),tail(introndata$midpoint,n=1),tail(intergenedata$midpoint,n=1))
     xMax = max(maxes)
     piMin = 0
     piMax = 0.052
     xPlaces = getAxis(xMax)
     xMax = tail(xPlaces[[1]], n=1)
 
     pdf(file=paste(path, 'pi_windows_',i,'.pdf', sep = ""),width=15, height=11)
     par(mfrow=c(2,2), oma=c(5,2,5,2), mar = c(1,4,1,2))
 
     par(xaxs='i')
     plot(ffolddata$midpoint, ffolddata$pi, main="4fold", type="l", ylab = "", xlab = "", ylim=c(piMin,piMax), xlim=c(0,xMax), frame=F, axes=F, col=139)
     axis(2)
     title(ylab="pi")
     abline(mean(ffolddata$pi),0,col="red", lty=2)
     par(mar = c(4,1,1,0))
     plot(zfolddata$midpoint, zfolddata$pi, main="0fold", type="l", ylab = "", xlab = "", ylim=c(piMin,piMax), xlim=c(0,xMax), frame=F, axes=F, col=139)
     axis(2)
     abline(mean(zfolddata$pi),0,col="red", lty=2)
     par(mar = c(4,4,1,2))
     plot(introndata$midpoint, introndata$pi, main="intron", type="l", ylab = "",xlab="", ylim=c(piMin,piMax), xlim=c(0,xMax), frame=F, axes=F, col=139)
     abline(mean(introndata$pi),0,col="red", lty=2)
     axis(2)
     axis(1, at=xPlaces[[1]], labels=xPlaces[[2]])
     title(xlab="position",ylab="pi")
     par(mar = c(4,1,1,0))
     #plotGenes(annot)
     #plotNC(annot, max(introndata$midpoint))
     plot(intergenedata$midpoint, intergenedata$pi, main="intergene", type="l", ylab = "",xlab="", ylim=c(piMin,piMax), xlim=c(0,xMax), frame=F, axes=F, col=139)
     abline(mean(intergenedata$pi),0,col="red", lty=2)
     axis(2)
     axis(1, at=xPlaces[[1]], labels=xPlaces[[2]])
     title(xlab="position")
     #plotGenes(annot)
     title(main=paste("Pi - Chromosome ",i,sep=""), outer=T, cex.main=5, font.main=3)
     dev.off()
 
     #----TajD plots----
     #dev.new(width=12, height=12)
     dMin = -3
     dMax = 0.5
     pdf(file=paste(path, 'tajd_windows_',i,'.pdf', sep = ""),width=15, height=11)
     par(mfrow=c(2,2), oma=c(5,2,5,2), mar = c(1,4,1,2))
 
     par(xaxs='i')
     plot(ffolddata$midpoint, ffolddata$tajd, main="4fold", type="l", ylab = "", xlab = "",xlim=c(0,xMax), ylim=c(dMin,dMax), frame=F, axes=F, col=139)
     axis(2)
     title(ylab="Taj D")
     abline(mean(ffolddata$tajd),0,col="red", lty=2)
     abline(0,0)
     par(mar = c(4,1,1,0))
     plot(zfolddata$midpoint, zfolddata$tajd, main="0fold", type="l", ylab = "", xlab = "",xlim=c(0,xMax), ylim=c(dMin,dMax), frame=F, axes=F, col=139)
     axis(2)
     abline(mean(zfolddata$tajd),0,col="red", lty=2)
     abline(0,0)
     par(mar = c(4,4,1,2))
     plot(introndata$midpoint, introndata$tajd, main="intron", type="l", ylab = "",xlab="", xlim=c(0,xMax), ylim=c(dMin,dMax), frame=F, axes=F, col=139)
     axis(2)
     axis(1, at=xPlaces[[1]], labels=xPlaces[[2]])
     title(xlab="position",ylab="Taj D")
     abline(mean(introndata$tajd),0,col="red", lty=2)
     abline(0,0)
     par(mar = c(4,1,1,0))
     #plotGenes(annot)
     #plotNC(annot, max(introndata$midpoint))
     plot(intergenedata$midpoint, intergenedata$tajd, main="intergene", type="l", ylab = "",xlab="", xlim=c(0,xMax), ylim=c(dMin,dMax), frame=F, axes=F, col=139)
     axis(2)
     axis(1, at=xPlaces[[1]], labels=xPlaces[[2]])
     title(xlab="position")
     abline(mean(intergenedata$tajd),0,col="red", lty=2)
     abline(0,0)
     #plotGenes(annot)
     title(main=paste("Tajima\'s D - Chromosome ",i,sep=""), outer=T, cex.main=5, font.main=3)
     dev.off()
 
     print(paste("Finished scaffold - ",i,sep=""))
 }
 
 print("Done")
