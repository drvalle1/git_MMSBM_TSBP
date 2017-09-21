boxplot(pi1)

plot(NA,NA,xlim=c(0,nloc),ylim=c(0,1))
for (i in 1:ngroups){
  lines(1:nloc,pi1[,i],col=i)
}

rango=range(c(true.lambda,lambda[1:6,1:6]))
plot(true.lambda,lambda[1:6,1:6],xlim=rango,ylim=rango)
lines(rango,rango)