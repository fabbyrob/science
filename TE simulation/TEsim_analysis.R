.mode csv
.output allSims2.csv
select ID, a.teNum, l, N, G, u, x, calc, t, cross, s, r, K, sel, c, b, e, z, v, generation, class, fitness, output.teNum, silence, silencediff, ne from (select * from simulation) as a left join output on ID = simID;

.mode csv
.output lastGen.csv
select ID, a.teNum, l, N, G, u, x, calc, t, cross, s, r, K, sel, c, b, e, z, v, generation, class, fitness, output.teNum, silence, silencediff, ne from (select * from simulation) as a left join output on ID = simID where (generation = 999 or generation = 4999);

.mode csv
.output finalFreqs.csv
select ID, a.teNum, l, N, G, u, x, calc, t, cross, s, r, K, sel, c, b, e, z, v, generation, freq from (select * from simulation) as a left join finalFreqs on ID = simID ;


#select ID, a.teNum, l, N, G, u, x, calc, t, cross, s, r, K, sel, c, b, generation, class, fitness, output.teNum, silence from (select * from simulation where calc = "ALL") as a left join output on ID = simID;
#.output homo.csv
#select ID, a.teNum, l, N, G, u, x, calc, t, cross, s, r, K, sel, c, b, generation, class, fitness, output.teNum, silence from (select * from simulation where calc = "HOMO") as a left join output on ID = simID;
#.output het.csv
#select ID, a.teNum, l, N, G, u, x, calc, t, cross, s, r, K, sel, c, b, generation, class, fitness, output.teNum, silence from (select * from simulation where calc = "HET") as a left join output on ID = simID;

for (k in 1:3){
if (k == 1){
	data = read.csv('/Users/wiliarj/Desktop/temp/allSims.csv',header=F)
	type = 'HET'
	pdf (file = paste('/Users/wiliarj/Desktop/temp/',type,'.pdf'))
	}
if (k == 2){
	data = read.csv('/Users/wiliarj/Desktop/temp/all.csv',header=F)
	type = 'ALL'
	pdf (file = paste('/Users/wiliarj/Desktop/temp/',type,'.pdf'))
	}
if (k == 3){
	data = read.csv('/Users/wiliarj/Desktop/temp/homo.csv',header=F)
	type = 'HOMO'
	pdf (file = paste('/Users/wiliarj/Desktop/temp/',type,'.pdf'))
	}

#for (j in 1:2){
names(data) = c("ID", "teNum", "l", "N", "G", "u", "x", "calc", "t", "cross", "s", "r", "K", "sel", "c", "b", "generation", "fitness", "meanteNum", "silence")
attach(data)
meanteNumRel = meanteNum/G



#x11()
#shapes = c(1)#ifelse(data$t==j, 1, 5)
#colors = ifelse(data$G == 100, "green", ifelse(data$G == 200, "blue", "black"))
#colors = ifelse(data$cross == 0, "black", ifelse(data$cross == 1, "purple", ifelse(data$cross == 2, "blue", ifelse(data$cross == 5, "darkgreen", ifelse(data$cross == 10, "green", ifelse(data$cross == 20, "yellow", ifelse(data$cross == 25, "orange", "red")))))))
colors = ifelse(data$l == 0.1, "red", ifelse(data$l == 0.5, "blue", "darkgreen"))
shapes = ifelse(data$l == 0.1, 1, ifelse(data$l == 0.5, 2, 3))
par(mfrow=c(2,3),oma=c(2,2,2,2))
lowr = data[data$r==.1,]
lowcolors = colors[data$r==.1]
lowshapes = shapes[data$r==.1]
highr = data[data$r==.9,]
highcolors = colors[data$r==.9]
highshapes = shapes[data$r==.9]
lowids = unique(lowr$ID)
highids = unique(highr$ID)
#plot(lowr$generation, lowr$fitness, col= lowcolors, pch = lowshapes,type="l")
#plot(highr$generation, highr$fitness, col = highcolors, pch = highshapes,type="l")
#plot(lowr$generation, lowr$meanteNum, col = lowcolors, pch = lowshapes,type="l")
#plot(highr$generation, highr$meanteNum, col  = highcolors, pch = highshapes,type="l")
#plot(lowr$generation, lowr$silence, col=lowcolors, pch = lowshapes,type="l")
#plot(highr$generation, highr$silence, col  = highcolors, pch = highshapes,type="l")

for (i in 1:length(lowids)){
	myIndexes = lowr$ID == lowids[i]
	sub = lowr[myIndexes,]
	subcol = lowcolors[myIndexes]
	subshape = lowshapes[myIndexes]
	if (i == 1){
		plot(sub$generation, sub$fitness, type="l", col=subcol, lty=subshape,ylim=c(.999,1))
	}else{
		points(sub$generation, sub$fitness, type="l", col=subcol, lty=subshape)
	}
}
title(main = "outcrosser, r = 0.1")

legend('bottomright', c("l=0.1", "l=0.5", "l=0.8"),lty=c(1,2, 3), col=c("red", "blue", "darkgreen"), bty='n')

for (i in 1:length(lowids)){
    myIndexes = lowr$ID == lowids[i]
    sub = lowr[myIndexes,]
	subcol = lowcolors[myIndexes]
	subshape = lowshapes[myIndexes]
    if (i == 1){
        plot(sub$generation, sub$meanteNum/(2*sub$G), type="l", col=subcol, lty=subshape,ylim=c(0,1))
    }else{
        points(sub$generation, sub$meanteNum/(2*sub$G), type="l", col=subcol, lty=subshape)
    }
}

slopes = c()
earlyslopes1 = c()
lateslopes1 = c()
earlyslopes2 = c()
lateslopes2 = c()
earlyslopes3 = c()
lateslopes3 = c()
for (i in 1:length(lowids)){
    myIndexes = lowr$ID == lowids[i]
    sub = lowr[myIndexes,]
	subcol = lowcolors[myIndexes]
	subshape = lowshapes[myIndexes]
    if (i == 1){
        plot(sub$generation, sub$silence, type="l", col=subcol, lty=subshape,ylim=c(0,1))
    }else{
        points(sub$generation, sub$silence, type="l", col=subcol, lty=subshape)
    }
	mod = coef(lm(sub$silence~sub$generation))[2]
	slopes = append(slopes, mod)
	
	subearly = sub[generation<20,]
	mod1 = coef(lm(subearly$silence~subearly$generation))[2]
	sublate = sub[generation>=20,]
	mod2 = coef(lm(sublate$silence~sublate$generation))[2]
	
	if (sub$l == 0.1){
		earlyslopes1 = append(earlyslopes1, mod1)
		lateslopes1 = append(lateslopes1, mod2)
	} else if(sub$l == 0.5){
		earlyslopes2 = append(earlyslopes2, mod1)
		lateslopes2 = append(lateslopes2, mod2)	
	} else if(sub$l == 0.8){
		earlyslopes3 = append(earlyslopes3, mod1)
		lateslopes3 = append(lateslopes3, mod2)	
	}
}
earlies = paste(format(mean(earlyslopes1), digits = 2), " ", format(mean(lateslopes1), digits = 2),"\n",format(mean(earlyslopes2), digits = 2), " ", format(mean(lateslopes2), digits = 2),"\n",format(mean(earlyslopes3), digits = 2), " ", format(mean(lateslopes2), digits = 2))
lates = paste(format(mean(lateslopes1), digits = 2),"\n",format(mean(lateslopes2), digits = 2),"\n",format(mean(lateslopes3), digits = 2))
title(main = earlies)
#abline(h=.1, col= "pink", lwd = 3)

for (i in 1:length(highids)){
	myIndexes = highr$ID == highids[i]
	sub = highr[myIndexes,]
	subcol = lowcolors[myIndexes]
	subshape = lowshapes[myIndexes]
	if (i == 1){
		plot(sub$generation, sub$fitness, type="l", col=subcol, lty=subshape,ylim=c(.999,1))
	}else{
		points(sub$generation, sub$fitness, type="l", col=subcol, lty=subshape)
	}
}
title(main = "selfer, r = 0.9")
for (i in 1:length(highids)){
    myIndexes = highr$ID == highids[i]
    sub = highr[myIndexes,]
	subcol = lowcolors[myIndexes]
	subshape = lowshapes[myIndexes]
    if (i == 1){
        plot(sub$generation, sub$meanteNum/(2*sub$G), type="l", col=subcol, lty=subshape,ylim=c(0,1))
    }else{
        points(sub$generation, sub$meanteNum/(2*sub$G), type="l", col=subcol, lty=subshape)
    }
}

slopes = c()
earlyslopes1 = c()
lateslopes1 = c()
earlyslopes2 = c()
lateslopes2 = c()
earlyslopes3 = c()
lateslopes3 = c()
for (i in 1:length(highids)){
    myIndexes = highr$ID == highids[i]
    sub = highr[myIndexes,]
	subcol = lowcolors[myIndexes]
	subshape = lowshapes[myIndexes]
    if (i == 1){
        plot(sub$generation, sub$silence, type="l", col=subcol, lty=subshape,ylim=c(0,1))
    }else{
        points(sub$generation, sub$silence, type="l", col=subcol, lty=subshape)
    }
	mod = coef(lm(sub$silence~sub$generation))[2]
	slopes = append(slopes, mod)
	
	subearly = sub[generation<20,]
	mod1 = coef(lm(subearly$silence~subearly$generation))[2]
	sublate = sub[generation>=20,]
	mod2 = coef(lm(sublate$silence~sublate$generation))[2]
	
	if (sub$l == 0.1){
		earlyslopes1 = append(earlyslopes1, mod1)
		lateslopes1 = append(lateslopes1, mod2)
	} else if(sub$l == 0.5){
		earlyslopes2 = append(earlyslopes2, mod1)
		lateslopes2 = append(lateslopes2, mod2)	
	} else if(sub$l == 0.8){
		earlyslopes3 = append(earlyslopes3, mod1)
		lateslopes3 = append(lateslopes3, mod2)	
	}
}
earlies = paste(format(mean(earlyslopes1), digits = 2), " ", format(mean(lateslopes1), digits = 2),"\n",format(mean(earlyslopes2), digits = 2), " ", format(mean(lateslopes2), digits = 2),"\n",format(mean(earlyslopes3), digits = 2), " ", format(mean(lateslopes2), digits = 2))
lates = paste(format(mean(lateslopes1), digits = 2),"\n",format(mean(lateslopes2), digits = 2),"\n",format(mean(lateslopes3), digits = 2))
title(main = earlies)#abline(h=.1, col= "pink", lwd = 3)
title(main = paste(type, ""), sub = "", outer=T)

#########
#
##x11()
#shapes = c(1)#ifelse(data$t==j, 1, 5)
##colors = ifelse(data$G == 100, "green", ifelse(data$G == 200, "blue", "black"))
#colors = ifelse(data$cross == 0, "black", ifelse(data$cross == 1, "purple", ifelse(data$cross == 2, "blue", ifelse(data$cross == 5, "darkgreen", ifelse(data$cross == 10, "green", ifelse(data$cross == 20, "yellow", ifelse(data$cross == 25, "orange", "red")))))))
#par(mfrow=c(2,3),oma=c(2,2,2,2))
#lowr = data[data$r==.1 &  data$s == .01 & data$l == 0.1 & data$u == .001 & data$t==j,]
#lowcolors = colors[data$r==.1 &  data$s == .01 & data$l == 0.1& data$u == .001 & data$t==j]
#lowshapes = shapes[data$r==.1 &  data$s == .01 & data$l == 0.1& data$u == .001 & data$t==j]
#highr = data[data$r==.9 &  data$s == .01 & data$l == 0.1& data$u == .001 & data$t==j,]
#highcolors = colors[data$r==.9 &  data$s == .01 & data$l == 0.1& data$u == .001 & data$t==j]
#highshapes = shapes[data$r==.9 &  data$s == .01 & data$l == 0.1& data$u == .001 & data$t==j]
#lowids = unique(lowr$ID)
#highids = unique(highr$ID)
##plot(lowr$generation, lowr$fitness, col= lowcolors, pch = lowshapes,type="l")
##plot(highr$generation, highr$fitness, col = highcolors, pch = highshapes,type="l")
##plot(lowr$generation, lowr$meanteNum, col = lowcolors, pch = lowshapes,type="l")
##plot(highr$generation, highr$meanteNum, col  = highcolors, pch = highshapes,type="l")
##plot(lowr$generation, lowr$silence, col=lowcolors, pch = lowshapes,type="l")
##plot(highr$generation, highr$silence, col  = highcolors, pch = highshapes,type="l")
#
#
#for (i in 1:length(lowids)){
#	myIndexes = lowr$ID == lowids[i]
#	sub = lowr[myIndexes,]
#	subcol = lowcolors[myIndexes]
#	subshape = lowshapes[myIndexes]
#	if (i == 1){
#		plot(sub$generation, sub$fitness, col = subcol,type="l",ylim=c(.999,1))
#	}else{
#		points(sub$generation, sub$fitness, col = subcol,type="l")
#	}
#}
#title(main = "outcrosser, r = 0.1")
#
#for (i in 1:length(lowids)){
#	myIndexes = lowr$ID == lowids[i]
#	sub = lowr[myIndexes,]
#	subcol = lowcolors[myIndexes]
#	subshape = lowshapes[myIndexes]
#	if (i == 1){
#		plot(sub$generation, sub$meanteNum, col = subcol,type="l",ylim=c(0,1000))
#	}else{
#		points(sub$generation, sub$meanteNum, col = subcol,type="l")
#	}
#}
#slopes = c()
#for (i in 1:length(lowids)){
#	myIndexes = lowr$ID == lowids[i]
#	sub = lowr[myIndexes,]
#	subcol = lowcolors[myIndexes]
#	subshape = lowshapes[myIndexes]
#	if (i == 1){
#		plot(sub$generation, sub$silence, col = subcol,type="l",ylim=c(0,1))
#	}else{
#		points(sub$generation, sub$silence, col = subcol,type="l")
#	}
#	mod = coef(lm(sub$silence~sub$generation))[2]
#	slopes = append(slopes, mod)
#}
#title(main = paste("Mean slope = ",format(mean(slopes), digits = 2)))
##abline(h=.1, col= "pink", lwd = 3)
#
#for (i in 1:length(highids)){
#	myIndexes = highr$ID == highids[i]
#	sub = highr[myIndexes,]
#	subcol = highcolors[myIndexes]
#	subshape = highshapes[myIndexes]
#	if (i == 1){
#		plot(sub$generation, sub$fitness, col = subcol,type="l",ylim=c(.999,1))
#	}else{
#		points(sub$generation, sub$fitness, col = subcol,type="l")
#	}
#}
#title(main = "selfer, r = 0.9")
#for (i in 1:length(highids)){
#	myIndexes = highr$ID == highids[i]
#	sub = highr[myIndexes,]
#	subcol = highcolors[myIndexes]
#	subshape = highshapes[myIndexes]
#	if (i == 1){
#		plot(sub$generation, sub$meanteNum, col = subcol,type="l",ylim=c(0,1000))
#	}else{
#		points(sub$generation, sub$meanteNum, col = subcol,type="l")
#	}
#}
#
#slopes = c()
#for (i in 1:length(highids)){
#	myIndexes = highr$ID == highids[i]
#	sub = highr[myIndexes,]
#	subcol = highcolors[myIndexes]
#	subshape = highshapes[myIndexes]
#	if (i == 1){
#		plot(sub$generation, sub$silence, col = subcol,type="l",ylim=c(0,1))
#	}else{
#		points(sub$generation, sub$silence, col = subcol,type="l")
#	}
#	mod = coef(lm(sub$silence~sub$generation))[2]
#	slopes = append(slopes, mod)
#}
#title(main = paste("Mean slope = ",format(mean(slopes), digits = 2)))
##abline(h=.1, col= "pink", lwd = 3)
#title(main = paste(type, "s = .01, l = .1, u = .001, t =", j), sub = "more red = more crossovers. ", outer=T)
#}
dev.off()
}
#####
##x11()
#shapes = ifelse(data$x == .001, 1, ifelse(data$x == .01, 22, ifelse(data$x == .05, 2, 3)))
#colors = ifelse(data$G == 100, "green", ifelse(data$G == 200, "blue", "black"))
#par(mfrow=c(2,2))
#lowr = data[data$r==.1 &  data$s == .01 & data$l == 0.1 & data$u == .001 ,]
#lowcolors = colors[data$r==.1 &  data$s == .01 & data$l == 0.1& data$u == .001 ]
#lowshapes = shapes[data$r==.1 &  data$s == .01 & data$l == 0.1& data$u == .001 ]
#highr = data[data$r==.9 &  data$s == .01 & data$l == 0.1& data$u == .001 ,]
#highcolors = colors[data$r==.9 &  data$s == .01 & data$l == 0.1& data$u == .001 ]
#highshapes = shapes[data$r==.9 &  data$s == .01 & data$l == 0.1& data$u == .001 ]
##plot(lowr$generation, lowr$fitness, col= lowcolors, pch = lowshapes)
##plot(highr$generation, highr$fitness, col = highcolors, pch = highshapes)
#plot(lowr$generation, lowr$meanteNum, col = lowcolors, pch = lowshapes)
#plot(highr$generation, highr$meanteNum, col  = highcolors, pch = highshapes)
#plot(lowr$generation, lowr$silence, col=lowcolors, pch = lowshapes)
#plot(highr$generation, highr$silence, col  = highcolors, pch = highshapes)
#
##x11()
#shapes = ifelse(data$G == 20, 1, ifelse(data$G == 100, 3, 15))
#colors = ifelse(data$G == 20, "green", ifelse(data$G == 100, "blue", "black"))
#par(mfrow=c(3,2))
#lowr = data[data$r==.1 &  data$s == .01 & data$l == 0.01 & data$u == .01 ,]
#lowcolors = colors[data$r==.1 &  data$s == .01 & data$l == 0.01& data$u == .01 ]
#lowshapes = shapes[data$r==.1 &  data$s == .01 & data$l == 0.01& data$u == .01 ]
#highr = data[data$r==.9 &  data$s == .01 & data$l == 0.01& data$u == .01 ,]
#highcolors = colors[data$r==.9 &  data$s == .01 & data$l == 0.01& data$u == .01 ]
#highshapes = shapes[data$r==.9 &  data$s == .01 & data$l == 0.01& data$u == .01 ]
#plot(lowr$generation, lowr$fitness, col= lowcolors, pch = lowshapes)
#plot(highr$generation, highr$fitness, col = highcolors, pch = highshapes)
#plot(lowr$generation, lowr$meanteNum, col = lowcolors, pch = lowshapes)
#plot(highr$generation, highr$meanteNum, col  = highcolors, pch = highshapes)
#plot(lowr$generation, lowr$silence, col=lowcolors, pch = lowshapes)
#plot(highr$generation, highr$silence, col  = highcolors, pch = highshapes)
#
##x11()
#
#par(mfrow=c(3,2))
#lowr = data[data$r==.1 &  data$s == .1 & data$l == 0.01,]
#lowcolors = colors[data$r==.1 &  data$s == .1 & data$l == 0.01]
#lowshapes = shapes[data$r==.1 &  data$s == .1 & data$l == 0.01]
#highr = data[data$r==.9 &  data$s == .1 & data$l == 0.01,]
#highcolors = colors[data$r==.9 &  data$s == .1 & data$l == 0.01]
#highshapes = shapes[data$r==.9 &  data$s == .1 & data$l == 0.01]
#plot(lowr$generation, lowr$fitness, col= lowcolors, pch = lowshapes)
#plot(highr$generation, highr$fitness, col = highcolors, pch = highshapes)
#plot(lowr$generation, lowr$meanteNum, col = lowcolors, pch = lowshapes)
#plot(highr$generation, highr$meanteNum, col  = highcolors, pch = highshapes)
#plot(lowr$generation, lowr$silence, col=lowcolors, pch = lowshapes)
#plot(highr$generation, highr$silence, col  = highcolors, pch = highshapes)
#
##x11()
#shapes = ifelse(data$G == 20, 1, ifelse(data$G == 100, 3, 15))
#colors = ifelse(data$G == 20, "green", ifelse(data$G == 100, "blue", "black"))
#par(mfrow=c(3,2))
#lowr = data[data$r==.1 &  data$s == .01 & data$l == 0.01.6,]
#lowcolors = colors[data$r==.1 &  data$s == .01 & data$l == 0.01.6]
#lowshapes = shapes[data$r==.1 &  data$s == .01 & data$l == 0.01.6]
#highr = data[data$r==.9 &  data$s == .01 & data$l == 0.01.6,]
#highcolors = colors[data$r==.9 &  data$s == .01 & data$l == 0.01.6]
#highshapes = shapes[data$r==.9 &  data$s == .01 & data$l == 0.01.6]
#plot(lowr$generation, lowr$fitness, col= lowcolors, pch = lowshapes)
#plot(highr$generation, highr$fitness, col = highcolors, pch = highshapes)
#plot(lowr$generation, lowr$meanteNum, col = lowcolors, pch = lowshapes)
#plot(highr$generation, highr$meanteNum, col  = highcolors, pch = highshapes)
#plot(lowr$generation, lowr$silence, col=lowcolors, pch = lowshapes)
#plot(highr$generation, highr$silence, col  = highcolors, pch = highshapes)
