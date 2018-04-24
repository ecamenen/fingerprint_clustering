rm(list=ls())
nb_met = 5
max_node = 10

#data_test
x = c(0, 3, 5, 8)
y = c(2, 0, 2, 6)
d_test = cbind(x, y)

names = paste("M", seq(1,nb_met-1),sep="")

#Plotting data
plot(d_test, xlim=c(0,max_node), ylim=c(0,max_node), font.lab=3, type='n',axes=F)
title(main="Metabolite network")
axis(1, seq(0,max_node), lwd=3, font.axis=3); axis(2, seq(0,max_node), lwd=3, font.axis=3);
abline(h=seq(0,max_node), v=seq(0,max_node), col="grey")
text(d_test, labels=names)
points(0,0,pch=21, col="red", cex=1.5, bg="yellow")

d_test = rbind(c(0,0), d_test)
names = paste("M", seq(0, nb_met-1), sep="")

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
plotRect = function (x1, y1, x2, y2) {
  rect(x1, y1, x2, y2, border = "red", lwd=2)
  rect(x1+0.03, y1, x2-0.03, y2-0.1, col = "white", border="white")
}


plot(d_test,xlim=c(0,nb_met), xaxt = "n", ylim=c(0,max(dist_matrix)+1), font.lab=3, type='n',axes=F)
# par("usr")[1] - 0.25 as the vertical placement
# srt = 45 as text rotation angle
# adj = 1 to place right end of text at tick mark
# xpd = TRUE to allow for text outside the plot region
text(seq(0,nb_met-1), par("usr")[1] - 1.5, srt = 45, adj = 1, labels = names, xpd = TRUE)
plotRect(0.5, 0, 1.5, 3)
plotRect(0, 0, 1, 2)
abline(h=0, col = "white", lwd=4)
axis(1, seq(0,nb_met-1), labels=FALSE, lwd=3, font.axis=3); axis(2, seq(0,max(dist_matrix)+1,5), lwd=3, font.axis=3)

