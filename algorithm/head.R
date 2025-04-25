library(statmod)
library(MASS)
library(expm)
library(Rcpp)

ngrid = 4
quad = gauss.quad(ngrid, "hermite")
b_list = quad$nodes
w_list = quad$weights

expit <- function(x) 1/(1+exp(-x))

mprod0 <- function(x,y) x%*%y

cppFunction('
  NumericVector prod(NumericMatrix x, NumericVector y){
    int l = x(_,0).size();
    NumericVector z (l);
    for (int i=0; i<l; i++){
      z[i] = sum(x(i,_)*y);
    }
    return z;
  }
')
cppFunction('
  NumericMatrix mprod(NumericMatrix x, NumericMatrix y){
    Function mprod0("mprod0");
    NumericMatrix z = mprod0(x,y);
    return z;
  }
')

bdiag0 <- function(...) as.matrix(bdiag(...))

cppFunction('
  NumericMatrix madd(NumericMatrix x, NumericMatrix y){
    Function add("add");
    NumericMatrix z = add(x,y);
    return z;
  }
')
cppFunction('
  NumericMatrix mneg(NumericMatrix x){
    Function mneg0("mneg0");
    NumericMatrix z = mneg0(x);
    return z;
  }
')

cppFunction('
  NumericMatrix diag0(NumericVector x){
    int l = x.size();
    NumericMatrix z (l);
    for (int i=0; i<l; i++){
      z(i,i) = x[i];
    }
    return z;
  }
')
