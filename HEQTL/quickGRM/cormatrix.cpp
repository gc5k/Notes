#include <Rcpp.h>
using namespace Rcpp;
#include <string>
#include <vector>
using namespace std;

// [[Rcpp::export]] 
NumericMatrix CorMatrix(NumericMatrix Mat)
{
  int i=0,j=0,k=0;
  int num;
  double sum,temp;
  int c=Mat.ncol();
  int r=Mat.nrow();
  double* m=(double*) malloc (sizeof(double) * c*r);
  memset(m, 0, sizeof(double)* c*r);
  for(i=0; i<r;i++)
  {
    for(j=0; j<c; j++)
    {
      m[i*c+j]=Mat(i,j);
    }
  }
  NumericMatrix exCor(r,r);
  for(i=0; i<r;i++)
  {
    for(j=0; j<=i; j++)
    {
      sum=0.0;
      num=0;
      for(k=0; k<c; k++)
      {
        //temp = Mat(i,k)*Mat(j,k);
        temp = m[i*c+k]*m[j*c+k];
        if(abs(temp)>1.0e-8) {num++;sum += temp;}
      }
      exCor(j,i) = exCor(i,j) = sum/(double)num;
    }
  }
  return wrap(exCor);
}

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


// [[Rcpp::export]]
NumericMatrix GRM(NumericMatrix G1, NumericMatrix G2) {
  NumericMatrix gM(G1.nrow(), G2.ncol());
  int len=G1.ncol();

  for(int i=0; i < G1.nrow(); i++) {
    for(int j = 0; j < G2.ncol(); j++) {
      for(int k = 0; k < len; k++) {
        gM(i,j) += G1(i,k) * G2(k,j);
      }
    }
  }

  return wrap(gM);
}
