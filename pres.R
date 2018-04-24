set.seed(1)
x=ceiling(runif(9,0,10))
y=ceiling(runif(9,0,10))
d_test=cbind(x,y)
d_test=rbind(d_test, c(0,0))
plot(d_test,xlim=c(0,10), ylim=c(0,10))

dist_t = dist(d_test, method = "euclidian")
