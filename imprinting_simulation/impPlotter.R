path = "/Users/wiliarj/Desktop/temp/"

out = read.csv(paste(path, "out.txt", sep = ""), header=T)
self = read.csv(paste(path, "self.txt", sep = ""), header=T)
mix = read.csv(paste(path, "mix.txt", sep = ""), header=T)

#plot maternal and paternal values for each

matCol = "red"
patCol = "blue"

ymin = min(out[,2], mix[,2], self[,2], out[,3], mix[,3], self[,3])
ymax = max(out[,2], mix[,2], self[,2], out[,3], mix[,3], self[,3])

plot(out$Generation, out$meanMaternal, ylim=c(ymin, ymax), col = matCol, lty=1, type = "l")
points(out$Generation, out$meanPaternal, ylim=c(ymin, ymax), col = patCol, lty=1, type = "l")

#points(mix$Generation, mix$meanMaternal, ylim=c(ymin, ymax), col = matCol, lty=2, type = "l")
#points(mix$Generation, mix$meanPaternal, ylim=c(ymin, ymax), col = patCol, lty=2, type = "l")

points(self$Generation, self$meanMaternal, ylim=c(ymin, ymax), col = matCol, lty=3, type = "l")
points(self$Generation, self$meanPaternal, ylim=c(ymin, ymax), col = patCol, lty=3, type = "l")

legend('bottomright', c("out","mix","self", "pat", "mat"), lty=c(1,2,3,1,1), col=c("black","black","black","blue", "red"), bty='n')
quartz()

ymin = min(out$meanPaternal-out$meanMaternal, self$meanPaternal-self$meanMaternal)
ymin = min(out$meanPaternal-out$meanMaternal, self$meanPaternal-self$meanMaternal)

plot(out$Generation, out$meanPaternal-out$meanMaternal, ylim=c(ymin, ymax), lty=1, type = "l")
points(self$Generation, self$meanPaternal-self$meanMaternal, ylim=c(ymin, ymax), lty=3, type = "l")

quartz()


ymin = min(out[,4], mix[,4], self[,4], out[,5], mix[,5], self[,5])
ymax = max(out[,4], mix[,4], self[,4], out[,5], mix[,5], self[,5])

plot(out$Generation, out$meanMhet, ylim=c(ymin, ymax), col = matCol, lty=1, type = "l")
points(out$Generation, out$meanPhet, ylim=c(ymin, ymax), col = patCol, lty=1, type = "l")

#points(mix$Generation, mix$meanMhet, ylim=c(ymin, ymax), col = matCol, lty=2, type = "l")
#points(mix$Generation, mix$meanPhet, ylim=c(ymin, ymax), col = patCol, lty=2, type = "l")

points(self$Generation, self$meanMhet, ylim=c(ymin, ymax), col = matCol, lty=3, type = "l")
points(self$Generation, self$meanPhet, ylim=c(ymin, ymax), col = patCol, lty=3, type = "l")

legend('bottomright', c("out","mix","self", "pat", "mat"), lty=c(1,2,3,1,1), col=c("black","black","black","blue", "red"), bty='n')

quartz()

plot(out$Generation, out$meanBirths, lty=1, type="l")
#points(mix$Generation, mix$meanBirths, lty=2, type="l")
points(self$Generation, self$meanBirths, lty=3, type="l")
