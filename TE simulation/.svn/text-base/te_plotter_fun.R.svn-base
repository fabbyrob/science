args = commandArgs(trailingOnly = TRUE)
if (length(args) > 0){
  base = ""
  myFile = args[1]
  prefix = args[2]
} else{
  base = "/Users/wiliarj/Desktop/temp/"
  myFile = paste(base,"allSims2.csv", sep="")
  prefix = ""
}

dataS = read.csv(myFile, header=F)
names(dataS) = c("ID", "teNum", "l", "N", "G", "u", "x", "calc", "t", "cross", "s", "r", "K", "sel", "c", "b", "e", "z", "v", "generation", "class", "fitness", "meanteNum", "silence", "silenceDiff", "ne")


########################remake wright and schoen 1999
#library(car)
data = dataS
last = max(data$generation)

#all here is a misnomer, it's only the ones I varied when I made this (i.e. the ones stephen varied)
#big matrix of all input values
#scatterplotMatrix(~r + u + teNum + N + x + calc + v, data = data, subset=generation==last, legend.plot=FALSE, reg.line=FALSE, spread=FALSE)

#big matrix of all inputs with all outputs
#scatterplotMatrix(~r + u + teNum + N + x + calc + v + fitness + meanteNum, data = data, subset=generation==last, legend.plot=FALSE, reg.line=FALSE, spread=FALSE)

samp = sample(unique(data$ID), size = 15, replace=FALSE)
cols=rainbow(length(unique(data$u)))
#scatterplot(meanteNum ~ generation | u, data=data, subset=ID %in% samp, grid=FALSE, reg.line=FALSE, smoother.args=list(span=0.8), cex=0.1, spread=TRUE,legend.plot=FALSE,col=cols, xlab="Generation", ylab="TE num", boxplot="y")
#legend("top",legend=unique(data$u),col=cols,lty=1, title="u")

subData= data[which((data$generation == 9999 & (data$calc == "ALL" | data$calc == "HOMO")) | (data$generation == 9999 & data$calc == "HET")),]
#N=1000
#scatterplot(meanteNum ~ r, data=subData, grid=FALSE, legend.plot=FALSE, subset=N==1000 & calc=="ALL", main="ALL", reg.line=FALSE, spread=TRUE, boxplot="")
#scatterplot(meanteNum ~ r, data=subData, grid=FALSE, legend.plot=FALSE, subset=N==1000 & calc=="HOMO", main="HOMO", reg.line=FALSE, spread=TRUE, boxplot="")
#scatterplot(meanteNum ~ r, data=subData, grid=FALSE, legend.plot=FALSE, subset=N==1000 & calc=="HET", main="HET", reg.line=FALSE, spread=TRUE, boxplot="")

cols = c("black", "red")
#scatterplot(meanteNum ~ r | v, data=subData, grid=FALSE, legend.plot=FALSE, subset=N==100 & calc=="ALL", main="ALL", reg.line=FALSE, spread=TRUE, boxplot="")
#legend("topleft", legend=c("v = 0", "v = 0.005"), col=cols, lty=1)
#scatterplot(meanteNum ~ r | v, data=subData, grid=FALSE, legend.plot=FALSE, subset=N==100 & calc=="HOMO", main="HOMO", reg.line=FALSE, spread=TRUE, boxplot="")
#legend("topleft", legend=c("v = 0", "v = 0.005"), col=cols, lty=1)
#scatterplot(meanteNum ~ r | v, data=subData, grid=FALSE, legend.plot=FALSE, subset=N==100 & calc=="HET", main="HET", reg.line=FALSE, spread=TRUE, boxplot="")
#legend("topleft", legend=c("v = 0", "v = 0.005"), col=cols, lty=1)
for (j in unique(dataS$b)){
  sD = dataS[dataS$b == j,]
  for (selection in unique(sD$sel)){
    sD1 = sD[sD$sel == selection,]
    for (N in unique(sD1$N)){
        sD2 = sD1[sD1$N == N,]
      for (mag in unique(sD2$e)){
          sD3 = sD2[sD2$e == mag,]
        for (silsel in unique(sD3$z)){
            sD4 = sD3[sD3$z == silsel,]
          for (silmut in unique(sD4$s)){
              sD5 = sD4[sD4$s == silmut,]
            for (u in unique(sD5$u)){
                sD6 = sD5[sD5$u == u,]
              for (l in unique(sD6$l)){
                subData2 = sD6[sD6$l == l,]
                subData2 = subData2[which(subData2$meanteNum > 0),]
                agData = aggregate(subData2, by=list(calc=subData2$calc, v=subData2$v, r=subData2$r), FUN=mean)
                agData2 = aggregate(subData2, by=list(calc=subData2$calc, v=subData2$v, r=subData2$r), FUN=sd)
                
                pdf(paste(base, prefix, "_", "TEs_N", N, "_e", mag, "_z", silsel, "_s", silmut, "_u", u,"_l",l, ".pdf", sep = ""), width=5, height=8)
                par(mfrow=c(3,1), oma=c(3, 3, 5, 3))
                for (calc in unique(agData$calc)){
                  subPoints = agData[which(agData$calc==calc),]
                  subError = agData2[which(agData2$calc==calc),]
                  plot(subPoints$r, subPoints$meanteNum, main=calc, xlab="selfing rate", ylab="TE #", pch=subPoints$col)
                  segments(subPoints$r, subPoints$meanteNum-subError$meanteNum, subPoints$r, subPoints$meanteNum+subError$meanteNum)
                  
                  mylm = lm(subPoints$meanteNum ~ subPoints$r)
                  abline(mylm)
                  text(par('usr')[1]+0.1, (par('usr')[4]+par('usr')[3])/2, summary(mylm)$r.squared)
                }
                title(paste("N = ", N, "e = ", mag, " z = ", silsel, " s = ", silmut, " u = ", u), outer=TRUE)
                dev.off()
              }
            }
          }
        }
      }
    }
  }
}

########################remake wright and schoen 1999

