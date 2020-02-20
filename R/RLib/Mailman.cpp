#include <Rcpp.h>
using namespace Rcpp;
#include <string>
#include <vector>
#include <cmath>
using namespace std;

// [[Rcpp::export]]
Rcpp::List MailmanProductL(NumericMatrix G1, NumericMatrix G2, NumericVector sg) {
  
  NumericMatrix mP(G1.nrow(), G2.ncol());
  NumericVector bs(G1.nrow());
  NumericVector Pcompact(G1.ncol());
  NumericMatrix Px(G2.nrow(), G2.ncol());
  
  for(int i = bs.length()-1; i>=0; i--) {
    if (i == (bs.length()-1) ) {
      bs[i] = 1;
    } else {
      bs[i] = bs[i+1] * sg.length();
    }
  }
  
  for(int i = 0; i < G1.ncol(); i++) {
    int idx = 0;
    for(int j = 0; j < G1.nrow(); j++) {
      idx += bs[j] * G1(j,i);
    }
    Pcompact[i]=idx;
  }
  
  for(int i=0; i < G2.ncol(); i++) {
    for(int j = 0; j < Pcompact.length(); j++) {
      Px(Pcompact[j],i) += G2(j,i); 
    }
    
    NumericVector z=Px(_,i);
    for(int j = 0; j < G1.nrow(); j++) {
      NumericVector z1(bs[j]);
      
      for(int l = 0; l < bs[j]; l++) {
        for(int k = 0; k < sg.length(); k++) {
          mP(j,i) += sg[k] * z[k*bs[j]+l];
          z1[l] +=z[k*bs[j]+l];
        }
      }
      z=z1;
    }
  }
  return Rcpp::List::create(Rcpp::Named("P")=Pcompact, 
                            Rcpp::Named("Px")=Px,
                            Rcpp::Named("mP")=mP,
                            Rcpp::Named("bs")=bs);
}

// [[Rcpp::export]]
Rcpp::List MailmanProductR(NumericMatrix G1, NumericMatrix G2, NumericVector sg) {
  
  NumericMatrix mP(G1.nrow(), G2.ncol());
  NumericVector bs(G1.ncol());
  NumericVector Pcompact(G1.nrow());
  NumericMatrix Ux(G1.nrow(), G2.ncol());
  
  // 
  //cal bs
  for(int i = bs.length()-1; i>=0; i--) {
    if (i == (bs.length()-1) ) {
      bs[i] = 1;
    } else {
      bs[i] = bs[i+1] * sg.length();
    }
  }
  
  //make pcompact
  for(int i = 0; i < G1.nrow(); i++) {
    int idx = 0;
    for(int j = 0; j < G1.ncol(); j++) {
      idx += bs[j] * G1(i,j);
    }
    Pcompact[i] = idx;
  }
  
  //cal Ux
  for(int i = 0; i < Ux.ncol(); i++) {
    
    NumericVector z(bs[bs.length()-1]);
    for(int j = bs.length() -1; j >=0; j--) {
      double v=G2(j,i);
      int len=bs[j]*sg.length();
      NumericVector z1(len);
      for(int k = 0; k < sg.length(); k++) {
        for(int l = 0; l < bs[j]; l++) {
          z1[k*bs[j] + l]=sg[k]*v + z[l];
        }
      }
      z=z1;
    }
    Ux(_,i)=z;
  }
  
  //cal mP
  for(int i = 0; i < Pcompact.length(); i++) {
    for(int j = 0; j < Ux.ncol(); j++) {
      mP(i,j) = Ux(Pcompact[i],j);
    }
  }
  
  return Rcpp::List::create(Rcpp::Named("P")=Pcompact, 
                            Rcpp::Named("Px")=Ux,
                            Rcpp::Named("mP")=mP,
                            Rcpp::Named("bs")=bs);
}
