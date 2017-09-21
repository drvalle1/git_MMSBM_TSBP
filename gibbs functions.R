sample.lambda=function(ngroups,nloc,z,dat,a.lambda,b.lambda){
  nitems=soma=matrix(0,ngroups,ngroups)
  for (i in 1:nloc){
    for (j in 1:nloc){
      ind1=z[i,j]  
      ind2=z[j,i]
      nitems[ind1,ind2]=nitems[ind1,ind2]+1
      soma[ind1,ind2]=soma[ind1,ind2]+dat[i,j]
    }
  }
  matrix(rgamma(ngroups*ngroups,soma+a.lambda,nitems+b.lambda),ngroups,ngroups)
}
#-----------------------------------
sample.v=function(z,nloc,ngroups,gamma){
  vpk=matrix(NA,nloc,ngroups)
  for (i in 1:(ngroups-1)){
    npk=rowSums(z==i)
    np.gk=rowSums(z>i)
    vpk[,i]=rbeta(nloc,npk+1,np.gk+gamma)
  }
  vpk[,ngroups]=1
  vpk
}
#-----------------------------------
v.to.pi=function(v,nloc,ngroups){
  pi=matrix(NA,nloc,ngroups)
  pi[,1]=v[,1]
  pi[,2]=v[,2]*(1-v[,1])
  for (i in 3:ngroups){
    pi[,i]=v[,i]*apply(1-v[,1:(i-1)],1,prod)
  }
  pi
}
#-----------------------------------
sample.z=function(nloc,ngroups,pi1,dat,lambda,combo){
  z=matrix(NA,nloc,nloc)
  lpi=log(pi1)
  llambda=log(lambda)
  diag1=diag(lambda)
  ldiag1=diag(llambda)

  for (i in 1:nloc){
    for (j in i:nloc){
      if (i==j) prob=lpi[i,]+dat[i,j]*ldiag1-diag1
      if (i!=j) prob=getprob(lpi1=lpi[i,],lpi2=lpi[j,],llambda=llambda,lambda=lambda,
                             nloc=nloc,ngroups=ngroups,dat12=dat[i,j],dat21=dat[j,i])
      prob=prob-max(prob)
      prob=exp(prob)
      prob=prob/sum(prob)
      ind=rmultinom(1,size=1,prob=prob)
      ind1=which(ind==1)
      
      if (i==j) z[i,j]=ind1
      if (i!=j) {
        z[i,j]=combo$gr1[ind1]
        z[j,i]=combo$gr2[ind1]
      }
    }
  }
  z
}
