#include <Rcpp.h>
#include <iostream>
#include <ctime>
#include <fstream>
using namespace Rcpp;

/***************************************************************************************************************************/
/*********************************                      UTILS          *****************************************************/
/***************************************************************************************************************************/

// [[Rcpp::export]]
NumericVector getprob(NumericVector lpi1,
                      NumericVector lpi2,
                      NumericMatrix llambda,
                      NumericMatrix lambda,
                      int nloc, int ngroups, int dat12, int dat21) {
  
  NumericVector prob(ngroups*ngroups);
  int oo=0;
  
  for(int i=0; i<ngroups;i++){
    for (int j=0; j<ngroups; j++){
      prob[oo]=lpi1[i]+lpi2[j]+dat12*llambda(i,j)-lambda(i,j)+dat21*llambda(j,i)-lambda(j,i);
      oo=oo+1;
    }
  }

  return (prob);
}
