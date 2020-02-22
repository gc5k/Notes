#include <Rcpp.h>
using namespace Rcpp;
#include <string>
#include <vector>
#include <cmath>
using namespace std;

//[[Rcpp::export(CorMatrixRcpp)]]
NumericMatrix CorMatrix(NumericMatrix Mat)
{
  int num;
  double sum,temp;
  int c=Mat.ncol();
  int r=Mat.nrow();
  NumericMatrix exCor(r,r);
  for(int i=0; i<r;i++)
  {
    for(int j=0; j<=i; j++)
    {
      sum=0.0;
      num=0;
      for(int k=0; k<c; k++)
      {
        temp = Mat(i,k)*Mat(j,k);
        if(abs(temp)>1.0e-8) {num++;sum += temp;}
      }
      exCor(j,i) = exCor(i,j) = sum/(double)num;
    }
  }
  return wrap(exCor);
}

//[[Rcpp::export(CorMatrix2Rcpp)]]
NumericMatrix CorMatrix2(NumericMatrix Mat1, NumericMatrix Mat2)
{
  try {
    if (Mat1.ncol() != Mat2.ncol()) {
      throw std::range_error("unequal size of freq and Dprime."); 
    }
  } catch (std::exception &ex) {
    forward_exception_to_r(ex);
  } catch (...) {
    ::Rf_error("Unequal ncol for matrix 1 and 2.");
  }

  int num;
  double sum,temp;
  int r1=Mat1.nrow();
  int r2=Mat2.nrow();
  int c=Mat1.ncol();
  NumericMatrix exCor(r1,r2);
  for(int i=0; i<r1;i++)
  {
    for(int j=0; j<r2; j++)
    {
      sum=0.0;
      num=0;
      for(int k=0; k<c; k++)
      {
        temp = Mat1(i,k)*Mat2(j,k);
        if(abs(temp)>1.0e-8) {num++;sum += temp;}
      }
      exCor(i,j) = sum/(double)num;
    }
  }
  return wrap(exCor);
}


inline static bool check(int l1, int l2, int diff) {
  try {
    if (l1 != (l2+diff)) {
      throw std::range_error("unequal size of freq and Dprime."); 
    }
  } catch (std::exception &ex) {
    forward_exception_to_r(ex);
  } catch (...) {
    ::Rf_error("length(freq) -length (Dprime) != 1, or too small freq.length (at least 2).");
  }
  return true;
}

// [[Rcpp::export(Dprime2LDRcpp)]]
NumericVector Dprime2LD(NumericVector freq, NumericVector Dp) {

  if(!check(freq.length(), Dp.length(), 1)) Rcpp::stop("failed");

  NumericVector ld(Dp.length());
  
  for(int i = 0; i < ld.length(); i++) {
    if (Dp[i]>0) {
      double d1=freq[i] * (1-freq[i+1]);
      double d2=(1-freq[i])*freq[i+1];
      ld[i] = Dp[i] * R::fmin2(d1, d2);
    } else {
      double d1=freq[i] * freq[i+1];
      double d2=(1-freq[i])* (1-freq[i+1]);
      ld[i] = Dp[i] * R::fmin2(d1, d2);
    }
  }
  return ld;
}

// [[Rcpp::export(Dprime2CorRcpp)]]
NumericVector Dprime2Cor(NumericVector freq, NumericVector Dp) {
  if(!check(freq.length(), Dp.length(), 1)) Rcpp::stop("failed");
  
  NumericVector ld=Dprime2LD(freq, Dp);
  NumericVector corl=NumericVector(freq.length()-1);
  for(int i = 0; i < corl.length(); i++) {
    corl[i] = ld[i]/sqrt(freq[i]*(1-freq[i])*freq[i+1]*(1-freq[i+1]));
  }
  return(corl);
}

// [[Rcpp::export(LD2CorRcpp)]]
NumericVector LD2Cor (NumericVector freq, NumericVector ld) {
  if(!check(freq.length(), ld.length(), 1)) Rcpp::stop("failed");

  NumericVector corl= NumericVector(freq.length()-1);
  for(int i = 0; i < corl.length(); i++) {
    corl[i] = ld[i]/sqrt(freq[i]*(1-freq[i])*freq[i+1]*(1-freq[i+1]));
  }
  return corl;
}

// [[Rcpp::export(GenerateGenoDprimeRcpp)]] 
NumericMatrix GenerateGenoDprime(NumericVector freq, NumericVector Dp, int N) {
  if(!check(freq.length(), Dp.length(), 1)) Rcpp::stop("failed");

  NumericMatrix G(N, freq.length());
  NumericVector ld = Dprime2LD(freq, Dp);

  for(int h = 0; h < N; h++) {
    NumericMatrix gMat = NumericMatrix(freq.length(), 2);
    
    for(int i = 0; i < gMat.nrow(); i++) {
      
      for(int j = 0; j < gMat.ncol(); j++) {
        if(i == 0) {
          gMat(i,j) = (int) Rcpp::rbinom(1, 1, 1-freq[i])[0];
        } else {
          double ff = freq[i-1] * freq[i];
          double hap_pro = (ff +ld[i-1])/freq[i-1];
          if (gMat(i-1,j)!=0) {
            ff = (1-freq[i-1])*(1-freq[i]);
            hap_pro = (ff +ld[i-1])/(1-freq[i-1]);
          }

          gMat(i,j) = gMat(i-1,j);
          if (Rcpp::runif(1, 0, 1)[0] > hap_pro) {
            gMat(i,j) = 1-gMat(i-1,j);
          }
        }
      }
    }

    for(int i = 0; i < freq.length(); i++) {
      G(h,i) = 2 - (gMat(i,0) + gMat(i,1));
    }
  }

  return G;
}

// [[Rcpp::export(GenerateGenoLDRcpp)]]
NumericMatrix GenerateGenoLD(NumericVector freq, NumericVector ld, int N) {
  if(!check(freq.length(), ld.length(), 1)) Rcpp::stop("failed");

  NumericMatrix G(N, freq.length());

  for(int h = 0; h < N; h++) {
    NumericMatrix gMat = NumericMatrix(freq.length(), 2);
    
    for(int i = 0; i < gMat.nrow(); i++) {
      
      for(int j = 0; j < gMat.ncol(); j++) {
        if(i == 0) {
          gMat(i,j) = (int) Rcpp::rbinom(1, 1, 1-freq[i])[0];
        } else {
          double ff = freq[i-1] * freq[i];
          double hap_pro = (ff +ld[i-1])/freq[i-1];
          if (gMat(i-1,j)!=0) {
            ff = (1-freq[i-1])*(1-freq[i]);
            hap_pro = (ff +ld[i-1])/(1-freq[i-1]);
          }
          
          gMat(i,j) = gMat(i-1,j);
          if (Rcpp::runif(1, 0, 1)[0] > hap_pro) {
            gMat(i,j) = 1-gMat(i-1,j);
          }
        }
      }
    }
    for(int i = 0; i < freq.length(); i++) {
      G(h,i) = 2- (gMat(i,0) + gMat(i,1));
    }
  }
  
  return G;
}

// [[Rcpp::export(GenerateHapLDRcpp)]]
NumericMatrix GenerateHapLD(NumericVector freq, NumericVector ld, int N) {

  if(!check(freq.length(), ld.length(), 1)) Rcpp::stop("failed");
  
  NumericMatrix G(N*2, freq.length());
  
  for(int h = 0; h < N; h++) {
    NumericMatrix gMat = NumericMatrix(freq.length(), 2);
    
    for(int i = 0; i < gMat.nrow(); i++) {
      
      for(int j = 0; j < gMat.ncol(); j++) {
        if(i == 0) {
          gMat(i,j) = (int) Rcpp::rbinom(1, 1, 1-freq[i])[0];
        } else {
          double ff = freq[i-1] * freq[i];
          double hap_pro = (ff +ld[i-1])/freq[i-1];
          if (gMat(i-1,j)!=0) {
            ff = (1-freq[i-1])*(1-freq[i]);
            hap_pro = (ff +ld[i-1])/(1-freq[i-1]);
          }
          
          gMat(i,j) = gMat(i-1,j);
          if (Rcpp::runif(1, 0, 1)[0] > hap_pro) {
            gMat(i,j) = 1-gMat(i-1,j);
          }
        }
      }
    }
    for(int i = 0; i < freq.length(); i++) {
      G(h*2,i) = 1 - gMat(i,0);
      G(h*2+1,i) = 1 - gMat(i,1);
    }
  }
  
  return G;
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
