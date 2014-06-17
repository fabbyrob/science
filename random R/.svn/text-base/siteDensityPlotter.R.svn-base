data = read.table('/Users/wiliarj/Desktop/temp/codingDensity.200kb.txt', header=T)

data$newMid = data$MIDPOINT

cols = rainbow(length(unique(data$CHROM)))
data$cols = cols[1]
offset = 0
ctr = 1
for (c in unique(data$CHROM)){
  data$newMid = ifelse(data$CHROM==c, data$newMid + offset, data$newMid)
  data$cols = ifelse(data$CHROM==c, cols[ctr], data$cols)
  offset = offset + 25000000
  ctr = ctr + 1
}

data$prop = data$COUNT/200000
par(mfrow=c(1,1))
for (t in unique(data$TYPE)){
  sData = subset(data, data$TYPE == t)
  
  ylims = c(min(sData$prop, na.rm=T), max(sData$prop, na.rm=T))
  xlims = c(min(sData$newMid, na.rm=T), max(sData$newMid, na.rm=T))
  
  for (c in unique(sData$CHROM)){
    ssData = subset(sData, sData$CHROM==c)
    if (c == "scaffold_1"){
      plot(ssData$newMid, ssData$prop, bty='n', xlab='position', xaxt='n', ylab='proportion of sites', ylim=ylims, xlim=xlims, type="l", col=ssData$cols)
      title(main=t)
    }
    else{
      points(ssData$newMid, ssData$prop, type="l", col=ssData$cols)
    }
  }
}

# for (c in unique(data$CHROM)){
#   sData = subset(data, data$CHROM==c & data$TYPE == "ALL")
#   ntData = subset(newTable, newTable$CHROM==c& newTable$TYPE == "ALL")
#   if (c == "scaffold_1"){
#     plot(sData$newMid, sData$prop, bty='n', xlab='position', xaxt='n', ylab='proportion of sites', ylim=ylims, xlim=xlims, type="l", col=sData$cols)
#     points(ntData$newMid, ntData$prop, type="l", col=ntData$cols, lty=2)
#     legend("bottomleft", legend=c("coding","intergenic"),lty=c(1,2),bty='n',ncol=1, cex=0.75, x.intersp=0.25, y.intersp = 2)
#   }
#   else{
#     points(sData$newMid, sData$prop, type="l", col=sData$cols)
#     points(ntData$newMid, ntData$prop, type="l", col=ntData$cols, lty=2)
#   }
# }