data = read.table("/Users/wiliarj/Desktop/temp/diversity.boot", header=T)
attach(data)

ymin = min(data$lowCI)
ymax = max(data$highCI)

plot(starting_dist, point_estimate, ylim = c(ymin, ymax), xlab = "distance from nearest non intergenic site", ylab="fraction of polymorphic sites")


for (i in c(1:length(starting_dist))){
    segments(starting_dist[i], lowCI[i], starting_dist[i], highCI[i])
}