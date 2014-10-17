data = read.table("/Users/williarj/Desktop/temp/2n.coverage.txt", header=T)

name = "fijiensis.2n.coverage.mean"
pdfer = 1
ratio = 2#0 = %homo and % het sep; 1 = ratio of het:homo; 2 = coverage
names = unique(data$ind)
if(pdfer == 1){
    pdf(paste("/Users/williarj/Desktop/temp/",name,".pdf", sep = ""), width = 12, height = 15)
}
if (pdfer == 2){
    jpeg(paste("/Users/williarj/Desktop/temp/",name,".jpg", sep = ""))
}
par(mfrow=c(10, 1), mar=c(3,3,1,1))
data$total = data$num_homo+data$num_het+data$num_homo_ref
gap = 2000000
win = 5000

maxs = c()
for (i in unique(data$pos)){
    sdata = data[data$pos==i,]
    maxs = append(maxs, max(sdata$midpoint))
}

scafs = unique(data$pos)
i = 0
tp = "p"
for (n in names){
    i = i + 1
    sdata = data[data$ind==n,]
    if (ratio == 0){
        ymax = max(sdata$num_homo/sdata$total, sdata$num_het/sdata$total, na.rm=T)
        #ymax = 1
    } else if (ratio == 1){
        temp = sdata$num_het/sdata$num_homo
        temp = temp[!(temp == Inf)]
        ymax = max(temp, na.rm=T)
        #ymax = 10
    } else{
        ymax = max(data$total/win)
    }
    xmax = max(sum(maxs)+3*gap)
    j = 0
    hets = 0
    hom = 0
    tot = 0
    for (s in scafs){
        ssdata = sdata[sdata$pos==s,]
        if (j == 0){
            if (ratio == 0){
                plot(ssdata$midpoint, ssdata$num_het/(ssdata$total), col="red", xlim=c(0, xmax), ylim=c(0, ymax), ylab = "", xaxt='n', bty='n', type=tp)
                points(ssdata$midpoint, ssdata$num_homo/(ssdata$total), col="blue", type=tp)
            }else if (ratio == 1){
                plot(ssdata$midpoint, ssdata$num_het/ssdata$num_homo, col="red", xlim=c(0, xmax), ylim=c(0, ymax), ylab = "", xaxt='n', bty='n', type=tp)
            }else {
                plot(ssdata$midpoint, ssdata$total/win, col="red", xlim=c(0, xmax), ylim=c(0, ymax), ylab = "", xaxt='n', bty='n', type=tp)    
            }
        } else{
            if (s == "scaffold_2"){
                const = maxs[1]+gap
            }
            if (s == "scaffold_3"){
                const = sum(maxs[1:2])+2*gap
            }
            if (s == "scaffold_4"){
                const = sum(maxs[1:3])+3*gap
            }
            if (s == "scaffold_5"){
                const = sum(maxs[1:4])+4*gap
            }
            if (s == "scaffold_6"){
                const = sum(maxs[1:5])+5*gap
            }
            if (s == "scaffold_7"){
                const = sum(maxs[1:6])+6*gap
            }
            if (s == "scaffold_8"){
                const = sum(maxs[1:7])+7*gap
            }
            if (ratio == 0){
                points(ssdata$midpoint+const, ssdata$num_het/(ssdata$total), col="red", xlim=c(0, xmax), ylim=c(0, ymax), ylab = "", xaxt='n', bty='n', type=tp)
                points(ssdata$midpoint+const, ssdata$num_homo/(ssdata$total), col="blue", type=tp)
            }else if (ratio == 1){
                points(ssdata$midpoint+const, ssdata$num_het/ssdata$num_homo, col="red", xlim=c(0, xmax), ylim=c(0, ymax), ylab = "", xaxt='n', bty='n', type=tp)
            } else {
                points(ssdata$midpoint+const, ssdata$total/win, col="red", xlim=c(0, xmax), ylim=c(0, ymax), ylab = "", xaxt='n', bty='n', type=tp)    
            }
        }
        j = 1
    }
    hets = sum(sdata$num_het)
    hom = sum(sdata$num_homo)
    tot = sum(sdata$total)
    text(x=(sum(maxs)+7*gap)/2, y=ymax*0.9, paste(n, "homos:", hom, "hets:", hets, "total:", tot), cex=1.2)
    #if (i %% 2 == 1){
   #     title(ylab="Counts")
    #}
}

if (ratio == 0){
    legend('topright', c('Homo', 'Het'), col=c("blue","red"), lty=c(1,1), bty='n')
} else if (ratio == 1){
    legend('topright', c('Het:Homo'), col=c("red"), lty=c(1,1), bty='n')
} else{
    legend('topright', c('Total'), col=c("red"), lty=c(1,1), bty='n')  
}
if (pdfer != 0){
    dev.off()
}
#
#if(pdfer == 1){
#    pdf(paste("/Users/wiliarj/Desktop/temp/",name,"2.pdf", sep = ""), width = 12, height = 15)
#}
#if (pdfer == 2){
#    jpeg(paste("/Users/wiliarj/Desktop/temp/",name,"2.jpg", sep = ""))
#}
#if (pdfer == 0){
#    quartz()
#}
#par(mfrow=c(6, 1), mar=c(3,3,1,1))
#for (n in names[8:13]){
#    i = i + 1
#    sdata = data[data$ind==n,]
#    if (ratio == 0){
#        ymax = max(sdata$num_homo/sdata$total, sdata$num_het/sdata$total, na.rm=T)
#        #ymax = .03
#    } else if (ratio == 1){
#        temp = sdata$num_het/sdata$num_homo
#        temp = temp[!(temp == Inf)]
#        ymax = max(temp, na.rm=T)
#        #ymax = 10
#    } else{
#        ymax = 1
#    }
#    xmax = max(sum(maxs)+7*gap)
#    j = 0
#    hets = 0
#    hom = 0
#    tot = 0
#    for (s in scafs){
#        ssdata = sdata[sdata$pos==s,]
#        if (j == 0){
#            if (ratio == 0){
#                plot(ssdata$midpoint, ssdata$num_het/(ssdata$total), col="red", xlim=c(0, xmax), ylim=c(0, ymax), ylab = "", xaxt='n', bty='n', type=tp)
#                points(ssdata$midpoint, ssdata$num_homo/(ssdata$total), col="blue", type=tp)
#            }else if (ratio == 1){
#                plot(ssdata$midpoint, ssdata$num_het/ssdata$num_homo, col="red", xlim=c(0, xmax), ylim=c(0, ymax), ylab = "", xaxt='n', bty='n', type=tp)
#            }else {
#                plot(ssdata$midpoint, ssdata$total/win, col="red", xlim=c(0, xmax), ylim=c(0, ymax), ylab = "", xaxt='n', bty='n', type=tp)    
#            }
#        } else{
#            if (s == "scaffold_2"){
#                const = maxs[1]+gap
#            }
#            if (s == "scaffold_3"){
#                const = sum(maxs[1:2])+2*gap
#            }
#            if (s == "scaffold_4"){
#                const = sum(maxs[1:3])+3*gap
#            }
#            if (s == "scaffold_5"){
#                const = sum(maxs[1:4])+4*gap
#            }
#            if (s == "scaffold_6"){
#                const = sum(maxs[1:5])+5*gap
#            }
#            if (s == "scaffold_7"){
#                const = sum(maxs[1:6])+6*gap
#            }
#            if (s == "scaffold_8"){
#                const = sum(maxs[1:7])+7*gap
#            }
#            if (ratio == 0){
#                points(ssdata$midpoint+const, ssdata$num_het/(ssdata$total), col="red", xlim=c(0, xmax), ylim=c(0, ymax), ylab = "", xaxt='n', bty='n', type=tp)
#                points(ssdata$midpoint+const, ssdata$num_homo/(ssdata$total), col="blue", type=tp)
#            }else if (ratio == 1){
#                points(ssdata$midpoint+const, ssdata$num_het/ssdata$num_homo, col="red", xlim=c(0, xmax), ylim=c(0, ymax), ylab = "", xaxt='n', bty='n', type=tp)
#            }else {
#                points(ssdata$midpoint+const, ssdata$total/win, col="red", xlim=c(0, xmax), ylim=c(0, ymax), ylab = "", xaxt='n', bty='n', type=tp)    
#            }
#        }
#        j = 1
#    }
#    hets = sum(sdata$num_het)
#    hom = sum(sdata$num_homo)
#    tot = sum(sdata$total)
#    text(x=(sum(maxs)+7*gap)/2, y=ymax*0.9, paste(n, "homos:", hom, "hets:", hets, "total:", tot), cex=1.2)
#    #if (i %% 2 == 1){
#    #     title(ylab="Counts")
#    #}
#}
#
#if (ratio == 0){
#    legend('topright', c('Homo', 'Het'), col=c("blue","red"), lty=c(1,1), bty='n')
#} else if (ratio == 1){
#    legend('topright', c('Het:Homo'), col=c("red"), lty=c(1,1), bty='n')
#} else{
#    legend('topright', c('Total'), col=c("red"), lty=c(1,1), bty='n')  
#}
#if (pdfer != 0){
#    dev.off()
#}