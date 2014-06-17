args = commandArgs(trailingOnly = TRUE)
if (length(args) > 0){
  base = ""
  myFile = args[1]
} else{
  base = "/Users/wiliarj/Desktop/temp/"
  myFile = paste(base,"allSims2.csv", sep="")
}

dataS = read.csv(myFile, header=F)
names(dataS) = c("ID", "teNum", "l", "N", "G", "u", "x", "calc", "t", "cross", "s", "r", "K", "sel", "c", "b", "e", "z", "v", "generation", "class", "fitness", "meanteNum", "silence", "silenceDiff", "ne")

#names(data) = c("ID", "teNum", "l", "N", "G", "u", "x", "calc", "t", "cross", "s", "r", "K", "sel", "c", "b", "e", "z", "v", "generation", "freq")


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




for (j in unique(data$b)){
  for (selection in unique(data$sel)){
    for (N in unique(data$N)){
      for (mag in unique(data$e)){
        for (silsel in unique(data$z)){
          for (silmut in unique(data$s)){
            for (u in unique(data$u)){
              subData2 = subData[subData$N==N & subData$e == mag & subData$z == silsel & subData$s == silmut & subData$u == u,]
              #subData2 = subData2[which(subData2$meanteNum > 0),]
              agData = aggregate(subData2, by=list(calc=subData2$calc, v=subData2$v, r=subData2$r), FUN=mean)
              agData2 = aggregate(subData2, by=list(calc=subData2$calc, v=subData2$v, r=subData2$r), FUN=sd)
              
              ylims = c(min(subData2$silence), max(subData2$silence))
              
              pdf(paste(base, "Silencing_N", N, "_e", mag, "_z", silsel, "_s", silmut, "_u", u, ".pdf", sep = ""), width=5, height=8)
              par(mfrow=c(3,1), oma=c(3, 3, 5, 3))
              for (calc in c("ALL","HOMO","HET")){
                subPoints = agData[which(agData$calc==calc),]
                subError = agData2[which(agData2$calc==calc),]
                mylm = lm(subPoints$silence ~ subPoints$r)
                plot(subPoints$r, subPoints$silence, main=calc, ylim=ylims, xlab="selfing rate", ylab="Silencing", pch=subPoints$col)
                abline(mylm)
                segments(subPoints$r, subPoints$silence-subError$silence, subPoints$r, subPoints$silence+subError$silence)
                text(par('usr')[1]+0.1, par('usr')[3]*2, format(summary(mylm)$r.squared, digits=2))
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

#dataS = read.csv("/Users/wiliarj/Desktop/temp/headSims.csv", header=F)
#names(dataS) = c("ID", "teNum", "l", "N", "G", "u", "x", "calc", "t", "cross", "s", "r", "K", "sel", "c", "b", "e", "z", "v", "generation", "class", "fitness", "meanteNum", "silence", "silenceDiff", "ne")

##plot the rates
for (j in unique(dataS$b)){
  for (selection in unique(dataS$sel)){
    for (N in unique(dataS$N)){
      for (mag in unique(dataS$e)){
        for (silsel in unique(dataS$z)){
          for (silmut in unique(dataS$s)){
            for (u in unique(dataS$u)){
            
              subData2 = dataS[dataS$N==N & dataS$e == mag & dataS$z == silsel & dataS$s == silmut & dataS$u == u,]
              #subData2 = subData2[which(subData2$meanteNum > 0),]
              #agData = aggregate(subData2, by=list(calc=subData2$calc, v=subData2$v, r=subData2$r), FUN=mean)
              #agData2 = aggregate(subData2, by=list(calc=subData2$calc, v=subData2$v, r=subData2$r), FUN=sd)
              
              
              pdf(paste(base, "Silencing_Slope_N", N, "_e", mag, "_z", silsel, "_s", silmut, "_u", u, ".pdf", sep = ""), width=5, height=8)
              par(mfrow=c(3,1), oma=c(3, 3, 5, 3))
              for (calc in c("ALL","HOMO","HET")){
                rs = c()
                coeffs = c()
                for (r in unique(subData2$r)){
                  subsubData = subData2[subData2$calc == calc & subData2$r == r & subData2$generation < 1000,]
                  reg = lm(subsubData$silence ~ subsubData$generation)
                  rs = append(rs, r)
                  coeffs = append(coeffs, coef(reg)[2])
                }
                
                mylm = lm (coeffs ~ rs)
                plot(rs, coeffs, main=calc, xlab="selfing rate", ylab="Silencing Slope")
                abline(mylm)
                text(par('usr')[1]+abs(par('usr')[1]*2), par('usr')[3]+abs(par('usr')[3]*0.5), format(summary(mylm)$r.squared, digits=2))
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


last1000 = subset(data, generation>9000)
agData = aggregate(last1000, by=list(ID = last1000$ID, calc = last1000$calc), FUN=mean)
agData$earliest = 9000

for (id in unique(data$ID)){
  eql = agData[agData$ID == id,]$silence
  gens = min(data[data$ID == id & data$silence > eql,]$generation)
  agData[agData$ID == id,]$earliest = gens
}


#plot means and stdev of replicates time to eq
for (j in unique(dataS$b)){
  for (selection in unique(dataS$sel)){
    for (N in unique(dataS$N)){
      for (mag in unique(dataS$e)){
        for (silsel in unique(dataS$z)){
          for (silmut in unique(dataS$s)){
            for (u in unique(dataS$u)){
              subData2 = agData[agData$N==N & agData$e == mag & agData$z == silsel & agData$s == silmut & agData$u == u,]
              agData2 = aggregate(agData, by=list(calc=agData$calc, r=agData$r), FUN=mean)
              agData2std = aggregate(agData, by=list(calc=agData$calc, r=agData$r), FUN=sd)
              
              pdf(paste(base, "Time_to_EQ_N", N, "_e", mag, "_z", silsel, "_s", silmut, "_u", u, ".pdf", sep = ""), width=5, height=8)
              par(mfrow=c(3,1))
              for (calc in c("ALL","HOMO","HET")){
                subPoints = agData2[which(agData2$calc==calc),]
                subError = agData2std[which(agData2std$calc==calc),]
                mylm = lm(subPoints$earliest ~ subPoints$r)
                ylims=c(min(subPoints$earliest-subError$earliest), max(subPoints$earliest-subError$earliest))
                plot(subPoints$r, subPoints$earliest, main=calc, ylims=ylims, xlab="selfing rate", ylab="Earliest generation \npast equilibrium", pch=20, col="black")
                abline(mylm, lty=2)
                segments(subPoints$r, subPoints$earliest-subError$earliest, subPoints$r, subPoints$earliest+subError$earliest)
                text(0, par('usr')[3]*1.2, format(summary(mylm)$r.squared, digits=2))
              }
              dev.off()

#plot means and stdev of replicates equilibrium
              pdf(paste(base, "Silence_EQ_N", N, "_e", mag, "_z", silsel, "_s", silmut, "_u", u, ".pdf", sep = ""), width=5, height=8)
              par(mfrow=c(3,1))
              for (calc in c("ALL","HOMO","HET")){
                subPoints = agData2[which(agData2$calc==calc),]
                subError = agData2std[which(agData2std$calc==calc),]
                mylm = lm(subPoints$silence ~ subPoints$r)
                ylims=c(min(subPoints$earliest-subError$earliest), max(subPoints$earliest-subError$earliest))
                plot(subPoints$r, subPoints$silence, main=calc, ylims=ylims, xlab="selfing rate", ylab="Equilibrium silencing", pch=20, col="black")
                abline(mylm, lty=2)
                segments(subPoints$r, subPoints$silence-subError$silence, subPoints$r, subPoints$silence+subError$silence)
                text(0, par('usr')[3]*1.2, format(summary(mylm)$r.squared, digits=2))
              }
              dev.off()
            }
          }
        }
      }
    }
  }
}

#plot silencing vs TE num

pdf(paste(base, "silence_v_tes2.pdf", sep=""), width=8, height=11)
par(mfrow=c(3,1))
for (calc in c("ALL","HOMO","HET")){
  subData = data[data$calc == calc,]
  plot(subData$meanteNum, subData$silence)
  lines(lowess(subData$meanteNum, subData$silence, f=1/3), col="blue")
  title(calc)
}
dev.off()