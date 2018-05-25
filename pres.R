rm(list=ls())
nb_met = 5
max_node = 10

#data_test
x = c(0, 3, 5, 8)
y = c(2, 0, 2, 6)
d_test = cbind(x, y)



#Plotting data
plotNetwork = function (d, m, n){
  names = paste("M", seq(1,n-1),sep="")
  plot(d_test, xlim=c(0,m), ylim=c(0,m), font.lab=3, type='n',axes=F)
  title(main="Metabolite network")
  axis(1, seq(0,m), lwd=3, font.axis=3); axis(2, seq(0,m), lwd=3, font.axis=3);
  abline(h=seq(0,m), v=seq(0,m), col="grey")
  text(d_test, labels=names)
  points(0,0,pch=21, col="red", cex=1.5, bg="yellow")

}

pdf("cou.pdf", width=7,height = 5)
plotNetwork(d_test, max_node, nb_met)
text(0.3,0.3, "M0")
dev.off()

d_test = rbind(c(0,0), d_test)
names = paste("M", seq(0, nb_met-1), sep="")

pdf("cou2.pdf", width=7,height = 5)
par(mar=c(2, 5, 4, 2) + 0.1 )
plot(0:14, 0:14, axes=F, type="n", ylab="")
title(ylab = "Distance between metabolites", font.lab=3, cex.lab=1.4)
axis(2,seq(0,14, 2), lwd=3, font.axis=3)
axis(2,seq(0,14, 1), lwd=3, font.axis=3, labels=rep("",15))
axis(1,seq(0,14, 3), lwd=3, font.axis=3,labels=names)
abline(h=seq(0,14, 2),col="grey")
dev.off()


#Distance matrix
dist_matrix = matrix(NA, nb_met, nb_met)
rownames(dist_matrix) <- names -> colnames(dist_matrix)
for (i in 1:nrow(d_test)){
  for (j in 1:nrow(d_test)){ 
    dist_matrix[i, j] = abs((d_test[i, 1] - d_test[j, 1]))  + abs((d_test[i, 2] - d_test[j, 2]))
  }
}
(dist_matrix = as.dist(dist_matrix))

#dendogramme
plotRect = function (x1, y1, x2, y2, c="red") {
  rect(x1, y1, x2, y2, border = c, lwd=2)
  rect(x1+0.03, y1, x2-0.03, y2-0.1, col = "white", border="white")
}

plot(d_test,xlim=c(0,nb_met), xaxt = "n", ylim=c(0,max(dist_matrix)+1), font.lab=3, xlab="Metabolites", ylab="Nb. nodes between metabolites", type='n',axes=F)
# par("usr")[1] - 0.25 as the vertical placement
# srt = 45 as text rotation angle
# adj = 1 to place right end of text at tick mark
# xpd = TRUE to allow for text outside the plot region
text(seq(0,nb_met-1), par("usr")[1] - 1.5, srt = 45, adj = 1, labels = names, xpd = TRUE)
#plotRect(0.5, 0, 2, 3, "orange")
plotRect(0, 0, 1, 2)
rect(3, 0, 2, 4, border="orange", lwd=2)
rect(2.98, 0, 2.02, 3.9, col = "white", border="white")
abline(h=0, col = "white", lwd=4)
axis(1, seq(0,nb_met-1), labels=FALSE, lwd=3, font.axis=3); axis(2, seq(0,max(dist_matrix)+1,5), lwd=3, font.axis=3)

