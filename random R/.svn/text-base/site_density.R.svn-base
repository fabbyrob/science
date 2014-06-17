
#path = "/Users/wiliarj/Desktop/temp/"
path = "/data/robert.williamson/finalGenome_analysis/rubella/cr_cg_mapping/"

sizes = read.table(paste(path, 'chrom_sizes.txt', sep= ""), header = F)
windows = c(50000)
for (window in windows){
    #window = 100000
    fileNameFinal = paste(path, "ratio_dist_",window,".pdf", sep = "")
    
    
    pdf(fileNameFinal, width = 15, height = 11)
    
    #quartz()
    
    par(mfrow=c(4,2))#, oma=c(5,0,5,0))
    
    for(j in 1:8){
        fixedsites = read.table(paste(path, 'scaf',j,'.26cg.0cr.sites', sep = ""), header=F)
        polysites = read.table(paste(path, 'scaf',j,'.sharedpoly.sites', sep = ""), header=F)
        
        names(fixedsites) = c('scaf','base','cg','cr')
        names(polysites) = c('scaf','base','cg','cr')
        
    
        
        fixed = c()
        poly = c()
        x = c()
        i = 1
        top5 = c(0,0,0,0,0)
        top5Pos = c(0,0,0,0,0)
        
        while (i < sizes[j,]){
            x1 = i+window/2
            i = i + window
            
            fixed1 = length(fixedsites[which(fixedsites$base >= i & fixedsites$base < i + window),][,1])
            poly1 = length(polysites[which(polysites$base >= i & polysites$base < i + window),][,1])
            
            fixed = append(fixed, fixed1/window)
            
            for (k in 1:5){
                if (fixed1 > top5[k]){
                    top5[k] = fixed1
                    top5pos = x1
                    break
                }
            }
            
            poly = append(poly, poly1/window)
            x = append(x, x1)
        }
        
        #ymax = max(fixed, poly)
        
        #plot(x, fixed, type = "l", lty = 2, ylim = c(0, 0.005), col = "blue", xlab = "position", ylab = paste("Fraction of sites in window (size=", window,")", sep = ""), cex=0.1)
        plot(x, fixed/poly, type = "l", lty = 1, col = "red", ylim=c(0,9), xlab = "position", ylab = paste("Ratio of fixed:polymorphic sites (window=",j,")", sep = ""), cex=0.1)
        #points(x, poly, type = "l", col = "red")
        #legend('topleft', c("shared polymorphism","fixed differences"), lty=c(1,2), col=c("red", "blue"),bty='n')
        #legend('topleft', c("fixed differences"), lty=c(2), col=c("blue"),bty='n')
        
        
        title(main = paste("Distribution of sites on chromosome", j))
        
        print(paste("Scaffold", j))
        print(top5)
        print(top5Pos)
        
    }
    dev.off()
}