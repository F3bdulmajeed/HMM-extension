
////////////////////////////////////////////////////////////
//
//Purpose: New HMM method for behavioural state identification from high frequency movement data
//
//Author:  Abdulmajeed Alharbi
//
//Date created: December 02 2025
//
////////////////////////////////////////////////////////////



// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rmath.h>
#include <limits>
using namespace Rcpp;
using namespace arma;

// ---------- Emission helpers (log-densities) ----------


// log Exp(rate = lambda)
// [[Rcpp::export]]
arma::vec dexp_cpp(const arma::vec& y, double lambda) {
  arma::vec safe_y = arma::clamp(y, 1e-10, DBL_MAX);
  return -lambda * safe_y + std::log(lambda);
}

// log Gamma(shape, rate)
// [[Rcpp::export]]
arma::vec dgamma_cpp(const arma::vec& x, double shape, double rate) {
  arma::vec safe_x = arma::clamp(x, 1e-10, DBL_MAX);
  arma::vec  log_x  = arma::log(safe_x);
  return shape * std::log(rate) - R::lgammafn(shape) + (shape - 1.0) * log_x - rate * safe_x;
}

// log soft-truncated von Mises
// [[Rcpp::export]]
arma::vec dsoft_tvm_cpp(const arma::vec& x, double k, double w) {
  arma::vec cosx   = arma::cos(x);
  arma::vec kernel = k * cosx + arma::log1p(-arma::exp(-w * (1.0 - cosx)));
  double I0k  = std::exp(k)     * R::bessel_i(k,     0, 2);
  double I0kw = std::exp(k + w) * R::bessel_i(k + w, 0, 2);
  double diff = I0k - std::exp(-w) * I0kw;
  if (diff <= 0.0) diff = std::numeric_limits<double>::min();
  double normal = std::log(diff);
  return kernel - normal - std::log(2.0 * M_PI);
}

// log von Mises, mean 0
// [[Rcpp::export]]
arma::vec dvm_cpp(const arma::vec& x, double kappa) {
  arma::vec cosx = arma::cos(x);
  double I0 = R::bessel_i(kappa, 0, 2);
  double logI0 = std::log(std::exp(kappa) * I0);
  return kappa * cosx - logI0 - std::log(2.0 * M_PI);
}



// Build log transition cube: for each t,
// [[Rcpp::export]]
arma::cube log_transition_probabilities(
    const arma::mat& X,
    const Rcpp::List& param_list
) {
  
  
  using arma::uword;
  const uword n = X.n_rows, p = X.n_cols;
  arma::mat first = Rcpp::as<arma::mat>(param_list[0]);
  const uword m  = first.n_rows, mm = m * m;
  
  if (param_list.size() != p)
    Rcpp::stop("Size mismatch: parameter_matrix should have %d elements (got %d).", (int)p, (int)param_list.size());
  
  arma::mat P(mm, p);
  for (uword j = 0; j < p; ++j) {
    arma::mat Mj = Rcpp::as<arma::mat>(param_list[j]);
    if (Mj.n_rows != m || Mj.n_cols != m)
      Rcpp::stop("parameter_matrix[%d] must be %dx%d.", (int)j, (int)m, (int)m);
    P.col(j) = arma::vectorise(Mj);
  }
  
  arma::mat S = P * X.t();
  arma::cube out(m, m, n);
  
  for (uword t = 0; t < n; ++t) {
    arma::mat T = arma::reshape(S.col(t), m, m);
    for (uword r = 0; r < m; ++r) {
      arma::rowvec row = T.row(r);
      const double rmax = row.max();
      row -= rmax;
      const double denom = arma::accu(arma::exp(row));
      if (denom > 0.0) T.row(r) = row - std::log(denom);
      else             T.row(r).fill(-std::numeric_limits<double>::infinity());
    }
    out.slice(t) = T;
  }
  return out;
}



// [[Rcpp::export]]
Rcpp::List forward_algorithm(
    const arma::cube& logA,
    const arma::vec&  logpi,
    const arma::mat&  logB
) {
  using arma::uword;
  const uword n = logB.n_rows, m = logB.n_cols;
  arma::mat log_alpha(n, m);
  log_alpha.row(0) = logpi.t() + logB.row(0);
  
  for (uword t = 1; t < n; ++t) {
    arma::mat M = logA.slice(t);
    M.each_col() += log_alpha.row(t-1).t();
    arma::vec cmax = arma::max(M, 0).t();
    M.each_row() -= cmax.t();
    arma::vec denom = arma::sum(arma::exp(M), 0).t();
    arma::rowvec lse = (cmax + arma::log(denom)).t();
    
    log_alpha.row(t) = logB.row(t) + lse;
  }
  
  const arma::rowvec last = log_alpha.row(n - 1);
  const double mx = last.max();
  const double loglik = mx + std::log(arma::accu(arma::exp(last - mx)));
  
  return Rcpp::List::create(
    Rcpp::Named("log_alpha")       = log_alpha,
    Rcpp::Named("log_likelihood")  = loglik
  );
}






// [[Rcpp::export]]
arma::uvec viterbi_cpp(const arma::cube& logA,
                       const arma::vec&  logpi,
                       const arma::mat&  logB) {
  using arma::uword;
  const uword n = logB.n_rows, m = logB.n_cols;
  
  arma::mat delta(n, m);
  arma::umat psi(n, m, fill::zeros);
  
  delta.row(0) = logpi.t() + logB.row(0);
  
  for (uword t = 1; t < n; ++t) {
    const arma::mat& A = logA.slice(t); 
    for (uword j = 0; j < m; ++j) {
      arma::vec cand = delta.row(t-1).t() + A.col(j);
      uword arg;
      double val = cand.max(arg);
      delta(t, j) = val + logB(t, j);
      psi(t, j)   = arg;
    }
  }
  
  arma::uvec path(n);
  uword last;
  delta.row(n-1).max(last);
  path(n-1) = last;
  
  for (int t = (int)n - 2; t >= 0; --t) {
    path(t) = psi(t+1, path(t+1));
  }
  return path + 1u;
}
