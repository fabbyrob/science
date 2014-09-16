library(ape)
data = read.table("/Users/wiliarj/Desktop/temp/scaf1.diff", header=T)
M <- matrix(0, 13, 13)
M[row(M) > col(M)] <- data$DIV
names = unique(unlist(list(unique(data$SAMP1), unique(data$SAMP2))))
rownames(M) <- colnames(M) <- names
tr = nj(M)
plot(tr, "u")