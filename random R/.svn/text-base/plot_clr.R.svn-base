##### FUNCTIONS #####
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

#####################

path = '/Users/wiliarj/Desktop/temp/'
#path = 'C:/Users/Fabby/Desktop/temp/'

grand = "blue"
rub = "red"

imprinted = read.table(paste(path, "imprinted_loci.txt", sep = ""), header = T)
icolor = "green"

pdfer = 1

if (pdfer == 1){pdf(file=paste(path, "clr.pdf",sep=""),width=11, height=8)}
par(mfrow=c(4,2))
for (i in c(1:8)){
    data = read.table(paste(path, 'scaf', i, '.clr.out',sep=""),header=T)
    data2 = read.table(paste(path, 'scaf', i, '.rub.clr.out',sep=""),header=T)

    ymax = max(max(data$LR), max(data2$LR))
    xmax = max(max(data$midpoint), max(data2$midpoint))
    
    plot(data$midpoint, data$LR, type = "l", col=grand,  xlim = c(0, xmax), xlab = "", ylab = "", xaxt='n', lwd=1)
    par(new=T)
    plot(data2$midpoint, data2$LR, type="l", col=rub, xlab = "", ylab = "", xaxt='n', yaxt='n', lwd=1)
    
    axis(4)
    xPlaces = getAxis(xmax)
    axis(1, at=xPlaces[[1]], labels=xPlaces[[2]])
    
    title(main=paste("Chromosome", i), xlab = "Position", ylab = "LR")
    #mtext("LR C. rubella", side=4, line=3)
    
    simprinted = imprinted[imprinted$scaf==i,]
    if (length(simprinted$gene) > 0){
        for (j in c(1:length(simprinted$gene))){
            #segments(x0=simprinted[j, 3], y0=.005, x1=simprinted[j, 4],y1=0.005, lwd = 5, col= icolor)
            mid = (simprinted[j,3]+simprinted[j,4])/2
            arrows(mid, ymax, mid, ymax/2, col = icolor)
        }
    }
    if (i == 7){
        rect(8000000,0, 8700000, ymax+.1, lwd = 1.5, border = "grey")
        legend('topleft', c("C. grandiflora", "C. rubella", "imprinted loci"), col=c(grand, rub, icolor),  bty='n', lty=rep(1, 3), lwd = c(1, 1, 1))
    }
    
}
if (pdfer == 1){dev.off()}
