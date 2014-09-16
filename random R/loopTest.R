plotGenes <- function(genes)
{
    for (i in 1:length(genes[,1])){
       segments(x0=genes[i,2], y0=.005, x1=genes[i,3],y1=0.005, lwd = 1, col= "red")
    
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

getAxis <- function(max)
{
    current = 0
    pos = c(0)
    labs = c(0)
    while(current < max){
        current = current+5*10^6
        pos = append(pos, current)
        lab = current/(10^6)
        labs = append(labs,paste(lab," Mb"))
    }
    pos = append(pos, max)
    labs = append(labs, "")
    list(pos,labs)
}

data = read.csv("/Users/Fabby/Desktop/scaf6.4fold.10050SNP.csv",header=T)
genes = read.table("/Users/Fabby/Desktop/scaf6.annot.sorted")
attach(data)

plot(MidPoint, pi, type="l",ylim=c(0.005,.04))
abline(avg,0,col="red",lty=2)
segments(10*10^6, .007, 14*10^6, .007, col="blue", lwd=5)
