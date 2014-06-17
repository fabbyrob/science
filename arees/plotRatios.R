real = read.table("/Users/wiliarj/Desktop/temp/real_data.out", header=T)
avg = read.table("/Users/wiliarj/Desktop/temp/average_syn_ratios.txt",header=T)
mu = read.table("/Users/wiliarj/Desktop/temp/mu_syn_ratios.txt",header=T)

pdf(file="/Users/wiliarj/Desktop/temp/theta_ks_ratios.pdf")
par(mfrow=c(3,1))

xs = c(0,2)
hist(real$ratio_syn, xlab="theta:Ks", xlim=xs, ylab="Number of loci", main="C. grandiflora")
hist(avg$ratio, xlab="theta:Ks", xlim=xs, ylab="Number of loci", main="Generated using average theta (equal mu)")
hist(mu$ratio, xlab="theta:Ks", xlim=xs, ylab="Number of loci", main="Generating using locus theta (different mu)")

xs = c(0,0.3)
hist(real$theta_syn, xlab="theta", xlim=xs, ylab="Number of loci", main="C. grandiflora")
hist(avg$theta, xlab="theta", xlim=xs, ylab="Number of loci", main="Generated using average theta (equal mu)")
hist(mu$theta, xlab="theta", xlim=xs, ylab="Number of loci", main="Generating using locus theta (different mu)")

xs = c(0,0.7)
hist(real$divergence_syn, xlab="ds", xlim=xs, ylab="Number of loci", main="C. grandiflora")
hist(avg$divergence, xlab="ds", xlim=xs, ylab="Number of loci", main="Generated using average theta (equal mu)")
hist(mu$divergence, xlab="ds", xlim=xs, ylab="Number of loci", main="Generating using locus theta (different mu)")
dev.off()