#include <Rcpp.h>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]
double cpp_med2(Rcpp::NumericVector xx) {
  Rcpp::NumericVector x = Rcpp::clone(xx);
  std::size_t n = x.size() / 2;
  std::nth_element(x.begin(), x.begin() + n, x.end());

  if (x.size() % 2) return x[n];
  return (x[n] + *std::max_element(x.begin(), x.begin() + n)) / 2.;
}

// inverse of empirical distribution (type 1 quantile definition in R)
// [[Rcpp::export]]
double cpp_quantile(Rcpp::NumericVector xx, double q) {
  Rcpp::NumericVector x = Rcpp::clone(xx);
  std::size_t n = std::floor(x.size() * q);
  if (n - x.size()*q==0){
    std::nth_element(x.begin(), x.begin() + n - 1, x.end());
    return x[n-1];
  } else {
    std::nth_element(x.begin(), x.begin() + n, x.end());
    return x[n];
  }
}

typedef std::vector<double> ints;
void insert( ints &cont, double value ) {
  ints::iterator it = std::lower_bound( cont.begin(), cont.end(), value, std::greater<double>() ); // find proper position in descending order
  cont.insert( it, value ); // insert before iterator it
}

// inverse of empirical distribution (type 1 quantile definition in R)
// [[Rcpp::export]]
NumericVector cpp_cumquantile(Rcpp::NumericVector xx, double q) {
  Rcpp::NumericVector x = Rcpp::clone(xx);
  double n = x.size();
  NumericVector quant(n);
  ints sort_x(1);
  sort_x[0] = x[0];
  quant[0] = x[0];
  for(int i = 1; i < n; ++i){
    double value = x[i];
    insert(sort_x, value); //x is sorted in decsending order!!!
    int pos = i+1-std::floor((i+1)*q);
    double pos1 = i+1-(i+1)*q;
    if (std::floor(pos1)==std::ceil(pos1)){
      quant[i] = sort_x[pos];
    } else {
      quant[i] = sort_x[pos-1];
    }
  }
  return quant;
}

// [[Rcpp::export]]
NumericVector cumsum_median_constrast_Cpp(NumericVector ts, String type){
  int n = ts.size();
  NumericVector result(n);
  int x = (type=="L");
  if(x!=0){
    for(int i = 0; i < n-1; ++i) {
      result[i] = cpp_med2(ts[Rcpp::Range(0, i)])-cpp_med2(ts[Rcpp::Range(i+1, n-1)]);
    }
    result[n-1] = 0;
  }
  if(x==0){
    for(int i = 1; i < n; ++i) {
      result[i] = cpp_med2(ts[Rcpp::Range(i, n-1)])-cpp_med2(ts[Rcpp::Range(0, i-1)]);
    }
    result[0] = 0;
  }
  return result;
}

// [[Rcpp::export]]
NumericVector cumsum_quantile_constrast_Cpp_check(NumericVector ts, String type, double q){
  int n = ts.size();
  NumericVector result(n);
  int x = (type=="L");
  if(x!=0){
    for(int i = 0; i < n-1; ++i) {
      result[i] = cpp_quantile(ts[Rcpp::Range(0, i)],q)-cpp_quantile(ts[Rcpp::Range(i+1, n-1)],q);
    }
    result[n-1] = 0;
  }
  if(x==0){
    for(int i = 1; i < n; ++i) {
      result[i] = cpp_quantile(ts[Rcpp::Range(i, n-1)],q)-cpp_quantile(ts[Rcpp::Range(0, i-1)],q);
    }
    result[0] = 0;
  }
  return result;
}

// [[Rcpp::export]]
NumericVector cumsum_quantile_constrast_Cpp(NumericVector ts, String type, double q){
  int n = ts.size();
  NumericVector result(n);
  int x = (type=="L");
  if(x!=0){
    //for(int i = 0; i < n-1; ++i) {
    //  result[i] = cpp_quantile(ts[Rcpp::Range(0, i)],q)-cpp_quantile(ts[Rcpp::Range(i+1, n-1)],q);
    //}
    //result[n-1] = 0;

    NumericVector result1 = cpp_cumquantile(ts, q);
    NumericVector ts2 = Rcpp::clone(ts);
    std::reverse(ts2.begin(),ts2.end());
    NumericVector result2 = cpp_cumquantile(ts2, q);
    std::reverse(result2.begin(), result2.end());
    result[Rcpp::Range(0, n-2)] = result1[Rcpp::Range(0, n-2)] - result2[Rcpp::Range(1, n-1)];
    result[n-1] = 0;
  }
  if(x==0){
    //for(int i = 1; i < n; ++i) {
    //  result[i] = cpp_quantile(ts[Rcpp::Range(i, n-1)],q)-cpp_quantile(ts[Rcpp::Range(0, i-1)],q);
    //}
    //result[0] = 0;
    NumericVector result1 = cpp_cumquantile(ts, q);
    NumericVector ts2 = Rcpp::clone(ts);
    std::reverse(ts2.begin(),ts2.end());
    NumericVector result2 = cpp_cumquantile(ts2, q);
    std::reverse(result2.begin(), result2.end());
    result[Rcpp::Range(1, n-1)] = result2[Rcpp::Range(1, n-1)] - result1[Rcpp::Range(0, n-2)];
    result[0] = 0;
  }
  return result;
}

// [[Rcpp::export]]
double cpp_acf(NumericVector ts) {
  NumericVector ydm = ts - mean(ts);
  double n = ydm.size();
  double acv0 = sum(pow(ydm, 2.0))/n;
  double acf1 = sum(ydm[Rcpp::Range(1, n-1)]*ydm[Rcpp::Range(0, n-2)])/n/acv0;
  return acf1;
}

// [[Rcpp::export]]
NumericVector cumsum_acf_constrast_Cpp(NumericVector ts, String type){
  int n = ts.size();
  NumericVector result(n);
  int x = (type=="L");
  if(x!=0){
    for(int i = 1; i < n-2; ++i) {
      result[i] = cpp_acf(ts[Rcpp::Range(0, i)])-cpp_acf(ts[Rcpp::Range(i+1, n-1)]);
    }
    result[0] = 0;
    result[n-2] = 0;
    result[n-1] = 0;
  }
  if(x==0){
    for(int i = 2; i < n-1; ++i) {
      result[i] = cpp_acf(ts[Rcpp::Range(i, n-1)])-cpp_acf(ts[Rcpp::Range(0, i-1)]);
    }
    result[0] = 0;
    result[1] = 0;
    result[n-1] = 0;
  }
  return result;
}


