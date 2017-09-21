# rm(list=ls(all=TRUE))
library('MCMCpack')
library('Rcpp')
set.seed(1)

setwd('U:\\independent studies\\MM SBM\\git_MMSBM_TSBP')
source('gibbs functions.R')
sourceCpp('aux1.cpp')

dat=read.csv('fake data.csv',as.is=T)
nloc=nrow(dat)
ngroups=6

#priors
gamma1=0.1

#initial values
lambda=matrix(rgamma(ngroups*ngroups,10,2),ngroups,ngroups)
pi1=rdirichlet(nloc,rep(1,ngroups))
z=matrix(sample(1:ngroups,size=nloc*nloc,replace=T),nloc,nloc)
combo=expand.grid(gr2=1:ngroups,gr1=1:ngroups)
a.lambda=b.lambda=0.1

ngibbs=1000
vec.lambda=matrix(NA,ngibbs,ngroups*ngroups)
vec.pi=matrix(NA,ngibbs,nloc*ngroups)

options(warn=2)
for (i in 1:ngibbs){
  print(i)
  lambda=sample.lambda(ngroups,nloc,z,dat,a.lambda,b.lambda)
  v=sample.v(z,nloc,ngroups,gamma1)
  pi1=v.to.pi(v,nloc,ngroups)
  z=sample.z(nloc,ngroups,pi1,dat,lambda,combo)
  
  #store results
  vec.lambda[i,]=lambda
  vec.pi[i,]=pi1
}
