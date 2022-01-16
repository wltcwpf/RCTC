#include <Rcpp.h>
#include <Rmath.h>

using namespace Rcpp;

//' Subset selection
//'
//' Select subsets of two components of ground motions.
//'
//' @param data1 One horizontal component of ground motions
//' @param data2 Another horizontal component of ground motions
//' @param fraction A factor to specify how much subset will be selected.
//' 0 means select all dataset.
//' @param length The minimum length of two input time sereis, integer type
//' @param time_dt Time step, double type
//' @param type_return Either 1 or 2. It specifies what does this \code{subset_select} return.
//' If set as 1, it will return a matirx (length is the length of selected subset,
//' width is 3) for subset. The first row is subset selection for \code{data1}, the second row is
//' subset selection for \code{data2}, while the third row is subset selection for time.
//' If set as 2, it will return only a double type value, which is threshold value for subset selection
//' @return \code{subset_select} returns two subsets of \code{data1} and \code{data2}, or a threshold value.
//' It depends on \code{type_return}.
//' @export
// [[Rcpp::export]]


Rcpp::NumericMatrix subset_select(Rcpp::NumericVector data1, Rcpp::NumericVector data2, double fraction, int length, double time_dt, int type_return){
  // it is used to select subset of sinc interpolated dataset, and return subset data points and corresponding time
  // input parameters are: two datasets, tne min length of datasets, and time step, type_return: 1 for data [first two rows are data, row 3 is time]; 2 for alevel
  double alevel, peak1, peak2, peak_min;
  int count = 0; // count is used to record the selected data points
  Rcpp::NumericMatrix subset(3,length);
  peak1 = 0;
  peak2 = 0;
  for(int i = 0; i < length; i++){
    if(peak1 < fabs(data1[i])){
      peak1 = fabs(data1[i]);
    }
    if(peak2 < fabs(data2[i])){
      peak2 = fabs(data2[i]);
    }
    if(peak1 < peak2){
      peak_min = peak1;
    } else{
      peak_min = peak2;
    }
  }
  alevel = fraction*peak_min;
  for(int i = 0; i < length; i++){
    if(sqrt(data1[i]*data1[i] + data2[i]*data2[i]) >= alevel){
      subset(0,count) = data1[i];
      subset(1,count) = data2[i];
      subset(2,count) = i*time_dt;
      count = count + 1;
    }
  }
  Rcpp::NumericMatrix result(3, count);
  for(int i = 0; i < count; i++){
    result(0,i) = subset(0,i);
    result(1,i) = subset(1,i);
    result(2,i) = subset(2,i);
  }
  if(type_return == 1){
    return result;
  } else{
    Rcpp::NumericMatrix level(1,1);
    level(0,0) = alevel;
    return level;
  }
}

