#include <RcppArmadillo.h>

//[[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

//[[Rcpp::depends(RcppArmadillo)]]
//[[Rcpp::export]]
List dynamic_tau_filter_cpp(
    const vec &ydata, const double omega, const double beta, 
    const double alpha1, const double alpha2, double f1, 
    const double kappa = 0.1, const double steep = 3, 
    const bool sharp = false) {
  // dynamic.tau.filter() filters the dynamic (1-kappa) quantile:
  // [note: so for kappa = 0.1 it is the UPPER tail quantile]
  //
  // INPUTS:
  // =======
  // ydata: vector of length T holding the data
  // omega: scalar intercept of the transition equation
  // beta: scalar loading of lagged tau
  // alpha1: scalar loading of exceedance indicator
  // alpha2: scalar loading of exceedance size
  // f1: initial value of the threshold
  // kappa: the right tail quantile percentage (0.1 for a 10% right tail prob)
  // steep: if > 0, then the steepness of the sigmoid logistic curve as an approximation of the step function
  // sharp: if TRUE, then use the sharp indicator function rather than the smooth sigmoid approximation
  //
  // OUTPUTS:
  // =======
  // llik: the tick-loss objective function
  // f_out: vector of lenght T holding the dynamic thresholds
  // 
  
  // initialize
  vec f_out = ydata + 0;
  double llik = 0;
  unsigned int dim_n = ydata.size();
  
  // check arguments
  if ((!sharp) && (steep <= 0)) {
    Rcerr << "\n\n*** ERROR IN dynamic_tau_filter_cpp(): sigmoid " <<
      "approximation to the indicator is flat or decreasing (steep = " << 
        steep << "\n";
  }
  
  // loop over the filter
  for (int i1=0; i1<dim_n; ++i1) {
    
    // store the filtered value
    f_out(i1) = f1;
    
    // update the objective
    double PoT = ydata(i1) - f1;
    double sigma = sharp ? ((PoT > 0) ? 1 : 0) : (1 / (1 + exp(-steep * PoT)));
    llik += - kappa * (1 - sigma) * PoT + (1 - kappa) * sigma * PoT;
    
    // update the dynamic threshold
    f1 = omega + beta * f1 -
      alpha1 * ( 
        - (1 - kappa) * sigma 
        + kappa * (1 - sigma)  
        - sigma * (1 - sigma) * PoT * steep
      ) +
      alpha2 * PoT;
  }
  
  // check on integrity of objective value
  if (!is_finite(llik)) llik = 1.0e20;

  return(List::create(
      Named("llik") = llik / dim_n,
      _["f_out"] = f_out
  ));
}
