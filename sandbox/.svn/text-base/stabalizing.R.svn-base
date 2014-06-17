name = "N1000_L100_u0.15_E300_K5000_o0.5_t3_replicate0"
max = 0.5
data = read.csv(paste("/Users/wiliarj/Desktop/temp/", name,".txt", sep=""), header=T)
initial = data[data$Population == "I",]
after = data[data$Population == "A",]

pdf(paste("/Users/wiliarj/Desktop/temp/",name,".pdf"))
par(mfrow=c(3,1))
plot(data$Generation, data$meanFitness, type="l", xaxt="n", xlab="", ylab="Mean Fitness")
plot(data$Generation, data$meanPhenotype, type="l", xlab="Generation", ylab="Phenotype")
abline(h=max)
abline(v=max(initial$Generation), lty=2, col="grey")
points(data$Generation, data$meanNeut, type="l", col="grey")
legend("topleft", c("selected", "neutral"), col=c("black", "grey"), lty=c(1,1), bty="n", cex = 0.75, x.intersp = 0.25, y.intersp=0.5)

plot(data$Generation, data$alpha, xlab="Generation", ylab = "alpha")
dev.off()

plotFitness = function(opt, loci){
  t = c(0.2,0.3,0.5,0.7,1,1.5,2,3,5)
  x = seq(-2,2,0.01)
  cols  = rainbow(length(t))
  
  maxD = max(abs(opt-loci*0.02), abs(opt-loci*(-0.02)))
  
  plot(x, (1-abs(opt-x)/maxD)^t[1], type="l", col=cols[1], xlab="Phenotype", ylab="Fitness")
  for (i in 2:length(t)){
    points(x, (1-abs(opt-x)/maxD)^t[i], type="l", col=cols[i])
  }
  abline(v=opt,lty=2)
  legend("topleft", t, t, lty=1, col=cols, cex=0.5, bty='n')
}