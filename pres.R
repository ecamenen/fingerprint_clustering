#data_test
set.seed(7)
x=ceiling(runif(9,0,10))
y=ceiling(runif(9,0,10))
d_test=cbind(x,y)

names=paste("M",seq(0,9),sep="")

#Plotting data
plot(d_test,xlim=c(0,10), ylim=c(0,10), font.lab=3, pch=21, col="blue",type='n',axes=F)
title(main="Metabolite network\n (in nb. of nodes)")
axis(1,seq(0,10),lwd=3,font.axis=3); axis(2,seq(0,10),lwd=3,font.axis=3)
abline(h=seq(0,10), v=seq(0,10), col="grey")
text(d_test, labels=names)
points(0,0,pch=21, col="red", cex=1.5, bg="yellow")

d_test = rbind(c(0,0),d_test)

#Distance matrix
dist_matrix = matrix(NA, 10,10)
rownames(dist_matrix) <- names -> colnames(dist_matrix)
for (i in 1:nrow(d_test)){
  for (j in 1:nrow(d_test)){ 
    dist_matrix[i, j] = abs((d_test[i, 1] - d_test[j, 1])  + (d_test[i, 2] - d_test[j, 2]))
  }
}
(dist_matrix = as.dist(dist_matrix))

