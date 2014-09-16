#######FUNCTIONS##########

#TODO combine these all into one function with separate functionality for each plot type

plotSilenceDiff <- function(IDs, classes, data, yMax, xAxis = FALSE, myTitle = FALSE, myylab = ""){
    for (i in 1:length(IDs)){
        currentID = IDs[i]
        for (j in 1:length(classes)){
            currentClass = classes[j]
            subData = data[data$ID == currentID & data$class == currentClass,]
            if (j == 1 & i == 1){
                myMar = c()
                if (xAxis == TRUE){
                    myMar = append(myMar, myXLabSpace)
                } else{
                    myMar = append(myMar, myXSpace)
                }
                if (myylab != ""){
                    myMar = append(myMar, myYLabSpace)
                } else{
                    myMar = append(myMar, myYSpace)
                }
                if (myTitle == TRUE){
                    myMar = append(myMar, myTitleSpace)
                } else{
                    myMar = append(myMar, myTSpace)
                }
                myMar = append(myMar, myRSpace)
                par(mar = myMar)
                
                
                plot(subData$generation, abs(subData$silenceDiff), col=subData$colors, ylim=c(0,yMax), xaxt = 'n', yaxt = 'n', xlab = '', ylab = myylab, cex.lab = 1.25, font.lab = 2)
                if (myylab != ""){
                    axis(2)
                } else{
                    axis(2, at = NULL, labels = FALSE)
                }
                if (xAxis == TRUE){
                    axis(1)
                    #title(xlab = 'Generation')
                } else{
                    axis(1, at = NULL, labels = FALSE)
                }
                if (myTitle == TRUE){
                    title("Silence Diff")
                }
                box()
            } else{
                points(subData$generation, subData$silenceDiff, col=subData$colors)
            }
        }
    }
    abline(h=mean(subData$silenceDiff), col = "yellow", lty=2, lwd=2)
}

plotSilenceLines <- function(IDs, classes, data, yMax, xAxis = FALSE, myTitle = FALSE, myylab = ""){
    totalSlope = 0
    counter = 0
    cols = rainbow(length(IDs))
    for (i in 1:length(IDs)){
        color = cols[i]
        currentID = IDs[i]
        for (j in 1:length(classes)){
            currentClass = classes[j]
            subData = data[data$ID == currentID & data$class == currentClass,]
            if (j == 1 & i == 1){
                myMar = c()
                if (xAxis == TRUE){
                    myMar = append(myMar, myXLabSpace)
                } else{
                    myMar = append(myMar, myXSpace)
                }
                if (myylab != ""){
                    myMar = append(myMar, myYLabSpace)
                } else{
                    myMar = append(myMar, myYSpace)
                }
                if (myTitle == TRUE){
                    myMar = append(myMar, myTitleSpace)
                } else{
                    myMar = append(myMar, myTSpace)
                }
                myMar = append(myMar, myRSpace)
                par(mar = myMar)
                
                
                plot(subData$generation, subData$silence, col=color, type="l", ylim=c(0, yMax), xaxt = 'n', yaxt = 'n', xlab = '', ylab = myylab, cex.lab = 1.25, font.lab = 2)
                if (myylab != ""){
                    axis(2)
                } else{
                    axis(2, at = NULL, labels = FALSE)
                }
                if (xAxis == TRUE){
                    axis(1)
                    #title(xlab = 'Generation')
                } else{
                    axis(1, at = NULL, labels = FALSE)
                }
                if (myTitle == TRUE){
                    title("Silencing")
                }
                box()
                
                totalSlope = totalSlope + lm(subData$silence[0:50]~subData$generation[0:50])$coefficients[2]
                counter = counter + 1
            } else{
                points(subData$generation, subData$silence, col=color, type="l")
            }
        }
    }
    meanSlope = totalSlope/counter
    #text(5, .2, paste("mS=", format(meanSlope, digits=2)))
}

plotTEpropLines <- function(IDs, classes, data, yMax, xAxis = FALSE, myTitle = FALSE, myylab = ""){
    yMax = max(data$meanteNum)
    cols = rainbow(length(IDs))
    for (i in 1:length(IDs)){
        currentID = IDs[i]
        color = cols[i]
        for (j in 1:length(classes)){
            currentClass = classes[j]
            subData = data[data$ID == currentID & data$class == currentClass,]
            if (j == 1 & i == 1){
                myMar = c()
                if (xAxis == TRUE){
                    myMar = append(myMar, myXLabSpace)
                } else{
                    myMar = append(myMar, myXSpace)
                }
                if (myylab != ""){
                    myMar = append(myMar, myYLabSpace)
                } else{
                    myMar = append(myMar, myYSpace)
                }
                if (myTitle == TRUE){
                    myMar = append(myMar, myTitleSpace)
                } else{
                    myMar = append(myMar, myTSpace)
                }
                myMar = append(myMar, myRSpace)
                par(mar = myMar)
                
                plot(subData$generation, subData$meanteNum, col=color, type="l", xaxt = 'n', yaxt = 'n', xlab = '', ylab = myylab, cex.lab = 1.25, font.lab = 2, ylim=c(0, yMax))
                if (myylab != ""){
                    axis(2)
                } else{
                    axis(2, at = NULL, labels = FALSE)
                }
                if (xAxis == TRUE){
                    axis(1)
                    #title(xlab = 'Generation')
                } else{
                    axis(1, at = NULL, labels = FALSE)
                }
                if (myTitle == TRUE){
                    title("Number of TEs")
                }
                box()
            } else{
                points(subData$generation, subData$meanteNum, col=color, type="l")
            }
        }
        #plot the total prop over time in black
        #subData = data[data$ID == currentID,]
        #agdata = aggregate(subData$meanteNum, FUN=sum, list(subData$generation))
        #print(max(agdata$x))
        #points(agdata$Group.1, agdata$x, col="black", type="l")
    }
}

plotFitnessLines <- function(IDs, classes, data, yMin, yMax, xAxis = FALSE, myTitle = FALSE, myylab = ""){
    for (i in 1:length(IDs)){
        currentID = IDs[i]
        for (j in 1:length(classes)){
            currentClass = classes[j]
            subData = data[data$ID == currentID & data$class == currentClass,]
            if (j == 1 & i == 1){
                myMar = c()
                if (xAxis == TRUE){
                    myMar = append(myMar, myXLabSpace)
                } else{
                    myMar = append(myMar, myXSpace)
                }
                if (myylab != ""){
                    myMar = append(myMar, myYLabSpace)
                } else{
                    myMar = append(myMar, myYSpace)
                }
                if (myTitle == TRUE){
                    myMar = append(myMar, myTitleSpace)
                } else{
                    myMar = append(myMar, myTSpace)
                }
                myMar = append(myMar, myRSpace)
                par(mar = myMar)
                
                plot(subData$generation, subData$fitness, type="l", ylim=c(yMin, yMax), xaxt = 'n', yaxt = 'n', xlab = '', ylab = myylab, cex.lab = 1.25, font.lab = 2)
                
                if (myylab != ""){
                    axis(2)
                } else{
                    axis(2, at = NULL, labels = FALSE)
                }
                if (xAxis == TRUE){
                    axis(1)
                    #title(xlab = 'Generation')
                } else{
                    axis(1, at = NULL, labels = FALSE)
                }
                if (myTitle == TRUE){
                    title(paste("Selfing rate =", data$r[1]))
                }
                #box()
            } else{
                if (j == 1){###TODO make this more efficient (hint: remove a for loop)
                    points(subData$generation, subData$fitness, type="l")
                }
            }
        }
    }
}

getClassNames <- function(numclasses){
    ClassList = c()
    for (i in 1:numclasses){
        name = paste("Class", i)
        ClassList = append(ClassList, name) 

    
    }
    return(ClassList)
    
}
    

summaryStats <- function(Low, Med, High){
    Lend = Low[Low$generation==999,]
    Mend = Med[Med$generation==999,]
    Hend = High[High$generation==999,]
    
    
#fitness
    min(Low$fitness)
    min(Med$fitness)
    min(High$fitness)
    
    mean(Lend$fitness)
    mean(Mend$fitness)
    mean(Hend$fitness)
    
    
    Lag = aggregate(Low$fitness, FUN=mean, list(Low$generation, Low$ID))
    Mag = aggregate(Med$fitness, FUN=mean, list(Med$generation, Med$ID))
    Hag = aggregate(High$fitness, FUN=mean, list(High$generation, High$ID))
    
    Ls = Lag[Lag$Group.1==999,]
    Ms = Mag[Mag$Group.1==999,]
    Hs = Hag[Hag$Group.1==999,]
    
    t = t.test(Ls$x, Hs$x)
    print(t$p.value)
    
#prop TEs
    
#total of each class, each generation, for each sim
    Lag = aggregate(Low$meanteNum, FUN=sum, list(Low$class, Low$generation, Low$ID))
    Ltot = aggregate(Low$meanteNum, FUN=sum, list(Low$generation, Low$ID))
    Mag = aggregate(Med$meanteNum, FUN=sum, list(Med$class, Med$generation, Med$ID))
    Mtot = aggregate(Med$meanteNum, FUN=sum, list(Med$generation, Med$ID))
    Hag = aggregate(High$meanteNum, FUN=sum, list(High$class, High$generation, High$ID))
    Htot = aggregate(High$meanteNum, FUN=sum, list(High$generation, High$ID))
    
#max of each class for each type
    for (i in c(0:4)){
        sub = Lag[Lag$Group.1==i,]
        print(i)
        print(max(sub$x))
        sub = sub[sub$Group.2==999,]
        print(mean(sub$x))
    }
    print("total")
    print(max(Ltot$x))
    sub = Ltot[Ltot$Group.1==999,]
    print("mean number of TEs low")
    print(mean(sub$x))
    Ls = sub
    
#max of each class for each type
    for (i in c(0:4)){
        sub = Mag[Mag$Group.1==i,]
        print(i)
        print(max(sub$x))
        sub = sub[sub$Group.2==999,]
        print(mean(sub$x))
    }
    print("total")
    print(max(Mtot$x))
    sub = Mtot[Mtot$Group.1==999,]
    print("mean number of TEs med")
    print(mean(sub$x))
    Ms = sub
    
#max of each class for each type
    for (i in c(0:4)){
        sub = Hag[Hag$Group.1==i,]
        print(i)
        print(max(sub$x))
        sub = sub[sub$Group.2==999,]
        print(mean(sub$x))
    }
    print("total")
    print(max(Htot$x))
    sub = Htot[Htot$Group.1==999,]
    print("mean number of TEs high")
    print(mean(sub$x))
    Hs = sub
    
    t = t.test(Ls$x, Hs$x)
    print(t$p.value)
    
#silencing
    Lag = aggregate(Low$silence, FUN=sum, list(Low$class, Low$generation, Low$ID))
    Ltot = aggregate(Low$silence, FUN=sum, list(Low$generation, Low$ID))
    Mag = aggregate(Med$silence, FUN=sum, list(Med$class, Med$generation, Med$ID))
    Mtot = aggregate(Med$silence, FUN=sum, list(Med$generation, Med$ID))
    Hag = aggregate(High$silence, FUN=sum, list(High$class, High$generation, High$ID))
    Htot = aggregate(High$silence, FUN=sum, list(High$generation, High$ID))
    
#max of each class for each type
    for (i in c(0:max(Low$b))){
        sub = Lag[Lag$Group.1==i,]
        print(i)
        print(max(sub$x))
        sub = sub[sub$Group.2==max(Lag$generation),]
        print(mean(sub$x))
    }
    print("total")
    print(max(Ltot$x))
    sub = Ltot[Ltot$Group.1==max(Lag$generation),]
    print(mean(sub$x))
    Ls = sub
    
#max of each class for each type
    for (i in c(0:max(Low$b))){
        sub = Mag[Mag$Group.1==i,]
        print(i)
        print(max(sub$x))
        sub = sub[sub$Group.2==max(Mag$generation),]
        print(mean(sub$x))
    }
    print("total")
    print(max(Mtot$x))
    sub = Mtot[Mtot$Group.1==max(Mag$generation),]
    print(mean(sub$x))
    Ms = sub
    
#max of each class for each type
    for (i in c(0:max(Low$b))){
        sub = Hag[Hag$Group.1==i,]
        print(i)
        print(max(sub$x))
        sub = sub[sub$Group.2==max(Hag$generation),]
        print(mean(sub$x))
    }
    print("total")
    print(max(Htot$x))
    sub = Htot[Htot$Group.1==max(Hag$generation),]
    print(mean(sub$x))
    Hs = sub
    
    t = t.test(Ls$x, Hs$x)
    print(t$p.value)
}
                
plotSelfingSummary <- function(data){

    par(mfrow=c(1,3))
    
    for (calc in c("ALL", "HET", "HOMO")){
        subData = data[data$calc == calc,]
        subData = subData[subData$generation==max(subData$generation),]
        agData = aggregate(subData$meanteNum, FUN=mean, list(subData$calc, subData$r))
        sdData = aggregate(subData$meanteNum, FUN=sd, list(subData$calc, subData$r))
        y = agData[agData$Group.1==calc,]
        s = sdData[sdData$Group.1==calc,]
        line = lm(y$x~y$Group.2)
        ymin = min(y$x-s$x)
        ymax = max(y$x+s$x)
        plot(y$Group.2, y$x, main = calc, xlab = "selfing rate", ylab = "mean end TE number", xlim = c(0,1), ylim = c(ymin, ymax))
        #abline(line, col = "red")
        for (i in c(1:length(s$Group.1))){
            segments(s$Group.2[i], y$x[i]-s$x[i], s$Group.2[i], y$x[i]+s$x[i])
        }
    }
    
}

plotNe <- function(data, ymax, ymin){
    agData = aggregate(data$ne, FUN=mean, by=list(data$generation))
    plot(agData$Group.1, agData$x, ylim = c(ymin, ymax), ylab="Ne", type="l")
}

plotSilTE <- function(IDs, classes, data, yMin, yMax, xAxis = FALSE, myTitle = FALSE, myylab = ""){
    for (i in 1:length(IDs)){
        currentID = IDs[i]
        for (j in 1:length(classes)){
            currentClass = classes[j]
            subData = data[data$ID == currentID & data$class == currentClass,]
            if (j == 1 & i == 1){
                myMar = c()
                if (xAxis == TRUE){
                    myMar = append(myMar, myXLabSpace)
                } else{
                    myMar = append(myMar, myXSpace)
                }
                if (myylab != ""){
                    myMar = append(myMar, myYLabSpace)
                } else{
                    myMar = append(myMar, myYSpace)
                }
                if (myTitle == TRUE){
                    myMar = append(myMar, myTitleSpace)
                } else{
                    myMar = append(myMar, myTSpace)
                }
                myMar = append(myMar, myRSpace)
                par(mar = myMar)
                
                plot(subData$silence, subData$meantenum, type="l", ylim=c(yMin, yMax), xaxt = 'n', yaxt = 'n', xlab = '', ylab = myylab, cex.lab = 1.25, font.lab = 2)
                
                if (myylab != ""){
                    axis(2)
                } else{
                    axis(2, at = NULL, labels = FALSE)
                }
                if (xAxis == TRUE){
                    axis(1)
                    title(xlab = 'Silencing rate')
                } else{
                    axis(1, at = NULL, labels = FALSE)
                }
                if (myTitle == TRUE){
                    title(paste("Selfing rate =", data$r[1]))
                }
                #box()
            } else{
                if (j == 1){###TODO make this more efficient (hint: remove a for loop)
                    points(subData$generation, subData$fitness, type="l")
                }
            }
        }
    }
}

#########

### MARGIN CONSTANTS ####

myTitleSpace = 2
myYLabSpace = 4
myXLabSpace = 3
myRFarSpace = 0

myYSpace = 4
myXSpace = 1
myTSpace = 1
myRSpace = 1


#######

### READ IN DATA #######
args = commandArgs(trailingOnly = TRUE)
if (length(args) > 0){
    base = ""
    myFile = args[1]
} else{
    base = "/Users/wiliarj/Desktop/temp/"
    myFile = paste(base,"allSims.csv", sep="")
}

if (! exists("dataS")){
    #base = "/data/robert.williamson/programs/basicSims.csv"

    #base = "/cap1/robert.williamson/TEsim/data/2013.07.05/allSims.csv"

    #base = "/Users/fabby/Desktop/allSims2.csv"

    dataS = read.csv(myFile,header=F)
    names(dataS) = c("ID", "teNum", "l", "N", "G", "u", "x", "calc", "t", "cross", "s", "r", "K", "sel", "c", "b", "e", "z", "v", "generation", "class", "fitness", "meanteNum", "silence", "silenceDiff", "ne")
}

data = dataS

attach(data)

ids = unique(data$ID)
colors = rainbow(length(unique(data$ID)))
data$colors = colors[1]
for (i in c(2:10)){
  data$colors = ifelse(data$ID == ids[i], colors[i], data$colors)
}

#TODO: generate a list of colors of the correct length based on the number of TE classes
#colors = c("red","orange","yellow","green","blue")
classes = unique(data$class[!is.na(data$class)])
colors = rainbow(length(classes))
#ClassNames = getClassNames(length(classes))
data$colors2 = colors[data$class]
data$TEprop = data$meanteNum/(2*data$G)

#########

###### PLOT #########

#values for B (max # classes)
for (j in unique(data$b)){
    for (selection in unique(data$sel)){
        for (trans in unique(data$u)){
            for (mag in unique(data$e)){
                for (silsel in unique(data$z)){
                    for (silmut in unique(data$s)){
                      for (l in unique(data$l)){
                          selectionName = ifelse(selection == "T", "Tournament", ifelse(selection == "O", "Roulette", ifelse(selection =="R", "Random", "Truncation")))
                          #fileName = ifelse(selection == "T", "Tournament", ifelse(selection == "O", "Roulette", "Truncation"))
                          #fileName = paste(fileName, "epsilon", me)
                          ClassNames = getClassNames(j)
                          subColors = colors[1:j]
                      
                          bSubset = data[data$l == l & data$b == j & sel == selection & u == trans & e == mag & z== silsel & cross == cross & s == silmut,]
                          
                          #base = "/data/robert.williamson/programs/"
  
                          #base = "/cap1/robert.williamson/TEsim/data/2013.07.05/"
                          #pdf(paste(base,"summary",j,selection,trans,mag,silsel,silmut, ".pdf", sep = ""), width=11, height=8)
                          #plotSelfingSummary(data)
                          #dev.off()
  
                          #base = "/cap1/robert.williamson/TEsim/data/2013.07.05/"
                          #pdf(paste(base,"summary",j,selection,trans,mag,silsel,silmut, ".pdf", sep = ""), width=11, height=8)
                          #plotSelfingSummary(data)
                          #dev.off()
  
                          
                          SilMax = max(bSubset$silence[!is.na(bSubset$silence)])
                          TEMax = max(bSubset$TEprop[!is.na(bSubset$TEprop)])
                          TEMax = max(bSubset$meanteNum[!is.na(bSubset$meanteNum)])
                          FitnessMax = max(bSubset$fitness[!is.na(bSubset$fitness)])
                          FitnessMin = min(bSubset$fitness[!is.na(bSubset$fitness)])
                          DiffMax = max(bSubset$silenceDiff[!is.na(bSubset$silenceDiff)])
                      
                          NeMax = max(bSubset$ne)
                          NeMin = min(bSubset$ne)
                          
                          #different selection types
                          for (mySel in unique(data$calc)){
                              fileName = selectionName
                              fileName = paste(paste(fileName, mySel, sep="_"), "b", j, "u", trans,"x", mag, "z", silsel, "s", silmut, "l", l, sep="_")
                              fileNameFinal = paste(base,"timeSeries_",fileName, ".pdf", sep="")
                          
                              print(fileNameFinal)
                              
                              pdf(fileNameFinal)
                              mySubSet = bSubset[bSubset$calc==mySel,]
                              Low = mySubSet[mySubSet$r==0.1,]
                              Med = mySubSet[mySubSet$r==0.5,]
                              High = mySubSet[mySubSet$r==0.9,]
                              
                      
                              
                              par(mfcol=c(3,3), oma=c(1,1,1,1))
                              
                              ##TODO make these plotting things into a function that gets called three times
                              #~~~~~ Low selfing plots
                              #plot(Low$generation, Low$fitness, main="Low selfing (0.1)")
                              #agdata = aggregate(Low$meanteNum/(2*Low$G), FUN=sum, list(Low$generation, Low$ID))
                              #plot(agdata$Group.1, agdata$x)
                              #plot(Low$generation, Low$meanteNum/(2*Low$G), col=Low$colors)
                              #plot(Low$generation, Low$silence, col=Low$colors)
                              LowIDs = unique(Low$ID)
                              plotFitnessLines(LowIDs, classes, Low, FitnessMin, FitnessMax, myylab = "Host Fitness", myTitle = TRUE)
                              
                              #TODO: add a function that plots the legend correctly based on the number of TE classes, generating the names and taking in the colors generated earlier (SEE other TODO). 
                              
                              plotTEpropLines(LowIDs, classes, Low, TEMax, myylab = "Number TEs")
                              plotSilenceLines(LowIDs, classes, Low, SilMax, myylab = "Silencing", xAxis = T)
                              #plotNe(Low, NeMax, NeMin)
                              #plotSilTE(LowIDs, classes, Low, 0, TEMax, myylab = "TE num", xAxis = T)
                              #plotSilenceDiff(LowIDs, classes, Low, DiffMax, xAxis = TRUE, myylab = "Silencing Diff")
                              
                              #~~~~~ medium selfing plots
                              #plot(Med$generation, High$fitness, main="Moderate selfing (0.5)")
                              #agdata = aggregate(Med$meanteNum/(2*Med$G), FUN=sum, list(Med$generation, Med$ID))
                              #plot(agdata$Group.1, agdata$x)
                              #plot(Med$generation, Med$meanteNum/(2*Med$G), col=Med$colors)
                              #plot(Med$generation, Med$silence, col=Med$colors)
                              MedIDs = unique(Med$ID)
                              plotFitnessLines(MedIDs, classes, Med, FitnessMin, FitnessMax, myTitle = TRUE)
                              plotTEpropLines(MedIDs, classes, Med, TEMax)
                              plotSilenceLines(MedIDs, classes, Med, SilMax, xAxis = T)
                              #plotNe(Med, NeMax, NeMin)
                              #plotSilTE(MedIDs, classes, Med, 0, TEMax, xAxis = T)
                              #plotSilenceDiff(MedIDs, classes, Med, DiffMax, xAxis = TRUE)
                              
                              #~~~~~ High selfing plots
                              #plot(High$generation, High$fitness, main="High selfing (0.9)")
                              #agdata = aggregate(High$meanteNum/(2*High$G), FUN=sum, list(High$generation, High$ID))
                              #plot(agdata$Group.1, agdata$x)
                              #plot(High$generation, High$meanteNum/(2*High$G), col=High$colors)
                              #plot(High$generation, High$silence, col=High$colors)
                              HighIDs = unique(High$ID)
                              plotFitnessLines(HighIDs, classes, High, FitnessMin, FitnessMax, myTitle = TRUE)
                              plotTEpropLines(HighIDs, classes, High, TEMax)
                              plotSilenceLines(HighIDs, classes, High, SilMax, xAxis = T)
                              #plotNe(High, NeMax, NeMin)
                              #plotSilTE(HighIDs, classes, High, 0, TEMax, myylab = "TE num", xAxis = T)
                              #plotSilenceDiff(HighIDs, classes, High, DiffMax, xAxis = TRUE)
                              
                              legend('topright', c(ClassNames, "total"), col = c(subColors, "black"), bty='n', lty=c(1,1,1,1,1), ncol = 2, cex = 0.7)
                              
                              mtext("Generation", cex = 1, font = 2, side = 1, outer = TRUE)
                      
                              if (mySel == "HOMO"){           
                                  name = 'homozygous'
                              }else if (mySel == "HET"){
                                  name = 'heterozygous'
                              } else{
                                  name = 'all'
                              }
                                  
                                  
                              #summaryStats(Low, Med, High)
                              #title(main=paste("Selection against", name, "TEs with", selectionName, "selection", "u =",trans, "x =",mag, " v=", silsel), font.main = 1, cex.main = 0.9, outer=T)
                              dev.off()
                            }
                        }
                    }
                }
            }
        }
    }

}

##########