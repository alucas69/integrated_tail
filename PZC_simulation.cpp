#include <RcppArmadillo.h>

//[[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

//[[Rcpp::export]]
vec LFZ0_obj_cpp(
    const vec &ydata, const vec &v_values, const vec &e_values, 
    const double alpha , const double smoothing = 0
) {
  int dim_n = ydata.size();
  vec obj = (ydata - v_values);
  if ((min(v_values) >= 0) || (min(e_values) >= 0)) {
    Rcout << "\n\n## WARNING: minimum VaR or EL non-negative (min(VaR) = " << min(v_values) << ", min(EL) = " << min(e_values) << ")";
    return(obj * 0 - 1.0e20 / dim_n);
  }
  if (smoothing <= 0) {
    // hard step function
    for (int i1=0; i1<dim_n;++i1) if (obj(i1) > 0) obj(i1) = 0;
  } else {
    // smoothed step function
    obj /= 1 + exp( (ydata - v_values) / smoothing );
  }
  obj =  obj / (e_values * alpha) + v_values / e_values + log(-e_values) - 1;
  return obj;
}




//[[Rcpp::depends(RcppArmadillo)]]
//[[Rcpp::export]]
List PZC_filter_cpp(const vec &ydata, const vec &omega, const mat &A_mat,
                const mat &B_mat, const vec &f0, const double alpha,
                double smoothing = 0) {

  unsigned int dim_n = ydata.size();
  mat f_values = zeros(dim_n, 2);
  vec f_current = f0 + 0;
  vec score = zeros(2);
  
  // run the filter
  for (int i1=0;i1<dim_n;++i1) {
    // store current value of time varying parameters
    f_values(i1,0) = f_current(0);
    f_values(i1,1) = f_current(1);

    // smooth the step function
    double stepper = (ydata(i1) - f_current(0));
    if (smoothing <= 0) stepper = (stepper <= 0) ? 1 : 0; else
      stepper = 1 / (1 + exp( stepper / smoothing));

    // update using (11)+(12) followed by (9)+(16)
    score(0) = -f_current(0) * (stepper - alpha);
    score(1) = stepper * ydata(i1) / alpha - f_current(1);
    f_current = omega + B_mat * f_current + A_mat * score;
  }

  // return results
  vec obj_t = LFZ0_obj_cpp(ydata, f_values.col(0), f_values.col(1), alpha, smoothing);
  double obj = mean(obj_t);
  if (!is_finite(obj)) obj = 1.0e20;
  return(List::create(
      Named("obj") = obj,
      _["obj_t"] = obj_t,
      _["ydata"] = ydata,
      _["VaR"] = f_values.col(0),
      _["EL"] = f_values.col(1)
  ));
}
