#include <Rcpp.h>
#include <Rmath.h>

using namespace Rcpp;

//' Response computation
//'
//' Compute response acceleration, pseudo-spectral acceleration (PSA), and spectral displacement
//'
//' @param data An acceleration time series array
//' @param period_t An array of oscillator periods
//' @param damping Damping ratio, expressed as a decimal (default value is 0.05)
//' @param time_dt Time step, in sec
//' @param type_return A dummy variable controlling output type, set to either 1 or 2.
//' If 1 is selected, it returns a two-row matrix with actual spectral acceleration
//' in the first row and PSA in the second row at each of the specified periods in \code{period_t}. If 2 is selected,
//' it returns a row vector of spectral displacements
//' @return The response acceleration and PSA, or spectral displacement. It depends on \code{type_return}.
//' @export
// [[Rcpp::export]]

Rcpp::NumericMatrix PS_cal_cpp(Rcpp::NumericVector data, Rcpp::NumericVector period_t, double damping, double time_dt, int type_return){
  // int type_return represents which type of results returns: 1 for PSA [2 by 111 periods]; 2 for relative disp [111 periods by data points]
  int len_per = period_t.size();
  int len_data = data.size();
  double period, omega_n, omega_D, k, A, B, C, D, AA, BB, CC, DD, Ra, max;
  Rcpp::NumericMatrix PS_disp(len_per,len_data);
  Rcpp::NumericMatrix PS_vel(len_per,len_data);
  Rcpp::NumericVector ACC(len_per);
  Rcpp::NumericVector PSA(len_per);
  for(int i=0; i<len_per; i++){
    PS_disp(i,0) = 0;
    PS_vel(i,0) = 0;
    period = period_t[i];
    omega_n = 2*4*atan(1.0)/period;
    omega_D = omega_n*sqrt(1-damping*damping);
    k = pow(omega_n,2);
    A = exp(-damping*omega_n*time_dt)*(damping/sqrt(1-pow(damping,2))*sin(omega_D*time_dt) + cos(omega_D*time_dt));
    B = exp(-damping*omega_n*time_dt)*(1/omega_D*sin(omega_D*time_dt));
    C = 1/k*(2*damping/omega_n/time_dt + exp(-damping*omega_n*time_dt)*(((1-2*pow(damping,2))/omega_D/time_dt-damping/sqrt(1-pow(damping,2)))*sin(omega_D*time_dt)-(1+2*damping/omega_n/time_dt)*cos(omega_D*time_dt)));
    D = 1/k*(1-2*damping/omega_n/time_dt+exp(-damping*omega_n*time_dt)*(((2*pow(damping,2)-1)/omega_D/time_dt)*sin(omega_D*time_dt)+2*damping/omega_n/time_dt*cos(omega_D*time_dt)));
    AA = -exp(-damping*omega_n*time_dt)*(omega_n/sqrt(1-pow(damping,2))*sin(omega_D*time_dt));
    BB = exp(-damping*omega_n*time_dt)*(cos(omega_D*time_dt)-damping/sqrt(1-pow(damping,2))*sin(omega_D*time_dt));
    CC = 1/k*(-1/time_dt+exp(-damping*omega_n*time_dt)*((omega_n/sqrt(1-pow(damping,2))+damping/time_dt/sqrt(1-pow(damping,2)))*sin(omega_D*time_dt)+1/time_dt*cos(omega_D*time_dt)));
    DD = 1/k/time_dt*(1-exp(-damping*omega_n*time_dt)*(damping/sqrt(1-pow(damping,2))*sin(omega_D*time_dt)+cos(omega_D*time_dt)));
    for (int j=1; j<len_data; j++){
      PS_disp(i,j) = A*PS_disp(i,j-1) + B*PS_vel(i,j-1) - C*data[j-1] - D*data[j];
      PS_vel(i,j) = AA*PS_disp(i,j-1) + BB*PS_vel(i,j-1) - CC*data[j-1] - DD*data[j];
      Ra = -2*damping*omega_n*PS_vel(i,j) - pow(omega_n,2)*PS_disp(i,j);
      if (ACC[i] < fabs(Ra)){
        ACC[i] = fabs(Ra);
      }
    }
    max=0;
    for(int n=0; n<len_data; n++){
      if(max < fabs(PS_disp(i,n))){
        max = fabs(PS_disp(i,n));
      }
    }
    PSA[i] = pow(omega_n,2)*max;
  }
  Rcpp::NumericMatrix result(2, PSA.size());
  for(int j=0;j<PSA.size();j++){
    for(int i=0;i<2;i++){
      if(i==0){
        result(i,j)=ACC[j];
      } else{
        result(i,j)=PSA[j];
      }
    }
  }
  if(type_return == 1){
    return result;
  } else {
    return PS_disp;
  }
}


