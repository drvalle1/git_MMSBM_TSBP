rm(list=ls(all=TRUE))
set.seed(1)
library('MCMCpack')

nnodes=100
ngroups=6
lambda=matrix(rgamma(ngroups*ngroups,10,1),ngroups,ngroups)

#introduce some sparseness
ind1=sample(1:ngroups,size=10,replace=T)
ind2=sample(1:ngroups,size=10,replace=T)
for (i in 1:10) lambda[ind1[i],ind2[i]]=0.1
image(lambda)
true.lambda=lambda

#create community membership
pi1=matrix(NA,nnodes,ngroups)
seq1=seq(from=0,to=ngroups,length.out=nnodes)
for (i in 1:ngroups) pi1[,i]=dnorm(seq1,mean=i,sd=0.2)

#normalize to sum to 1
for (i in 1:nnodes) pi1[i,]=pi1[i,]/sum(pi1[i,])
plot(NA,NA,xlim=c(0,nnodes),ylim=c(0,1))
for (i in 1:ngroups) lines(pi1[,i],col=i)
true.pi1=pi1
#------------------------------------------------
z=matrix(NA,nnodes,nnodes)
y=matrix(NA,nnodes,nnodes)
for (i in 1:nnodes){
  print(i)
  for (j in i:nnodes){
    if (i!=j){
      zi=rmultinom(1,size=1,prob=pi1[i,])
      ind.i=apply(zi==1,2,which)
      z[i,j]=ind.i
      
      zj=rmultinom(1,size=1,prob=pi1[j,])
      ind.j=apply(zj==1,2,which)
      z[j,i]=ind.j
      
      lambda1=lambda[ind.i,ind.j]
      y[i,j]=rpois(1,lambda1)
      lambda1=lambda[ind.j,ind.i]
      y[j,i]=rpois(1,lambda1)
    }
    if (i==j){
      zi=rmultinom(1,size=1,prob=pi1[i,])
      ind.i=apply(zi==1,2,which)
      z[i,j]=ind.i
      lambda1=lambda[ind.i,ind.i]
      y[i,j]=rpois(1,lambda1)
    }
  }
}
image(y)
z.true=z

setwd('U:\\independent studies\\MM SBM\\TSBP')
write.csv(y,'fake data.csv',row.names=F)
