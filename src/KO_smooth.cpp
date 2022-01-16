#include <Rcpp.h>
#include <Rmath.h>

using namespace Rcpp;

//' Konno Ohmachi smoothing
//'
//' Apply Konno Ohmachi smoothing to the Fourier Amplitude Spectra (FAS)
//'
//' @useDynLib RCTC
//' @importFrom Rcpp sourceCpp
//' @param freq An array of the frequency
//' @param amp An array of the corresponding FAS
//' @param b The coefficient of bandwidth. A smaller value will lead to more smoothing
//' @param rate The truncation rate of KO smoothing. A larger value will be more accurate but more time consuming
//' @return The smoothed FAS
//' @references Konno, K., Ohmachi, T. Ground-motion characteristics estimated from spectral ratio between horizontal and vertical components of microtremor.
//' Bulletin of the Seismological Society of America. 1998.88 (1): 228-241.
//' @export
// [[Rcpp::export]]
NumericVector ko_smooth(NumericVector freq, NumericVector amp, float b=20, float rate=2.5) {

  double fmul, fdiv, f1, f2, a1, a2, c1;
  int i, j;
  fmul = pow(10.0, rate / b);
  fdiv = 1.0 / fmul;
  int len_f = freq.length();
  NumericVector smoothed(len_f);
  // f/fc is the same as the ratio of frequency indexes j and i
  // Store log i, sin i, and cos i of indexes to speed up internal loops
  NumericVector sins(len_f);
  NumericVector coss(len_f);
  NumericVector logs(len_f);
  for (i = 0; i < len_f; i++) {
    double val = b*log10(freq[i]);
    logs[i] = val;
    sins[i] = sin(val);
    coss[i] = cos(val);
  }
  for (i = 0; i < len_f; i++) {
    f1 = freq[i] * fdiv;
    f2 = freq[i] * fmul;
    a1 = 0.0;
    a2 = 0.0;
    double sinfc = sins[i];
    double cosfc = coss[i];
    double logfc = logs[i];
    for (j = 0; j < len_f; j++) {
      if((freq[j] >= f1) & (freq[j] <= f2)){
        if(freq[j] == freq[i]){
          a1 += amp[j];
          a2 += 1;
        }else{
          c1 = (sins[j] * cosfc - sinfc * coss[j]) / (logs[j] - logfc);
          c1 *= c1; // square
          c1 *= c1; // square again: total c1^4
          a1 += c1 * amp[j];
          a2 += c1;
        }
      }
    }
    smoothed[i] = a1 / a2;
  }
  return smoothed;
}
