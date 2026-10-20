// src/arpls_solver.cpp
//
// Compiled ARPLS fitting engine. This is the sole source of truth for the
// IRLS reweighting loop; R/fit-arpls.R's `arpls_engine_r()` is a pure-R
// reference implementation kept for the benchmarking vignette and as a
// fallback, and must implement the identical convergence criterion below
// if it is ever modified -- the two drifted apart once already (L1 vs L2
// relative weight change), which this file now matches exactly to R's
// `sum(abs(w - w_old)) / sum(abs(w_old))`.
//
// NOTE on sparsity: the system matrix (W + lambda * D'D + ridge * I) is
// built here as a dense arma::mat, not kept sparse via arma::sp_mat /
// arma::spsolve(). An earlier version of this file attempted the sparse
// path; arma::spsolve() requires SuperLU to be compiled in and linked,
// which is not available in this build environment (Windows/MinGW/Rtools,
// no system SuperLU installed), and stood up as a new build dependency
// this close to submission was judged too risky to pursue. This is a
// disclosed, deliberate limitation: the compiled engine is faster than
// the interpreted R reference implementation (same dense algorithm,
// compiled arithmetic + BLAS dense solve vs. R's interpreted loop), but
// does not exploit the pentadiagonal structure of D'D the way a sparse or
// banded solve would. Revisiting this would require either a SuperLU-
// enabled build or a hand-written banded solver for the specific
// second-order difference structure -- documented here as future work,
// not claimed as delivered.
//
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <cmath>

using namespace Rcpp;

// [[Rcpp::export]]
List arpls_fit_cpp(const arma::vec& signal,
                   double lambda,
                   double ratio,
                   int max_iter,
                   double ridge) {
  int n = signal.n_elem;

  // Second-order sparse difference operator, D: (n-2) x n.
  arma::sp_mat D(n - 2, n);
  for (int i = 0; i < n - 2; ++i) {
    D(i, i)     =  1.0;
    D(i, i + 1) = -2.0;
    D(i, i + 2) =  1.0;
  }

  // Densified here -- see file header note on sparsity.
  arma::mat penalty = lambda * arma::mat(D.t() * D);

  arma::vec weights(n, arma::fill::ones);
  arma::vec old_weights(n, arma::fill::ones);
  arma::vec baseline = signal;
  bool converged = false;
  int iter = 0;

  for (iter = 0; iter < max_iter; ++iter) {
    old_weights = weights;

    arma::mat lhs = penalty;
    lhs.diag() += weights;
    lhs.diag() += ridge;

    arma::vec rhs = weights % signal;

    // Exception-throwing form (2-argument arma::solve): on a singular or
    // near-singular system this throws std::runtime_error, which Rcpp
    // converts into a catchable R condition. Do not switch this to the
    // status-returning overload (bool solve(vec&, mat, vec)) -- that form
    // fails silently and defeats the tryCatch() robustness built around
    // gcv_select_lambda().
    baseline = arma::solve(lhs, rhs);

    arma::vec residual = signal - baseline;

    arma::uvec neg_idx = find(residual < 0);
    if (neg_idx.n_elem < 2) {
      weights.fill(1.0);
    } else {
      double location = mean(residual(neg_idx));

      double sum_sq = 0.0;
      for (unsigned int k = 0; k < neg_idx.n_elem; ++k) {
        sum_sq += std::pow(residual(neg_idx(k)) - location, 2);
      }
      double variance = sum_sq / (neg_idx.n_elem - 1);
      double scale = std::sqrt(variance);

      // Threshold matches the R reference implementation's
      // `arpls_update_weights()` exactly (1e-10, not a bare `<= 0`): a
      // near-zero (not just exactly-zero) scale still drives the exponent
      // below toward overflow, which a plain `<= 0` guard would miss.
      if (!std::isfinite(scale) || scale <= 1e-10) {
        weights.fill(1.0);
      } else {
        arma::vec exponent = 2.0 * (residual - (2.0 * scale - location)) / scale;
        for (int j = 0; j < n; ++j) {
          if (exponent(j) > 500.0) exponent(j) = 500.0;
          weights(j) = 1.0 / (1.0 + std::exp(exponent(j)));
        }
      }
    }

    // Matches the R reference implementation's criterion exactly:
    // L1 relative change, normalised by the OLD weights (not the new
    // ones). A prior version of this file used an L2 norm normalised by
    // the new weights, which silently diverged from the R engine's
    // iteration count and final weights near the convergence boundary.
    double denominator = arma::accu(arma::abs(old_weights));
    double relative_change = 0.0;
    if (denominator > 0.0) {
      relative_change = arma::accu(arma::abs(weights - old_weights)) / denominator;
    }

    if (relative_change < ratio) {
      converged = true;
      break;
    }
  }

  // Explicit NumericVector conversion: an arma::vec returned directly
  // inside a List is auto-wrapped by Rcpp as an n x 1 R matrix, not a
  // plain atomic vector. That silently breaks the "identical in shape to
  // the R reference implementation" contract documented in
  // arpls_engine()/arpls_engine_r() -- it happened to go unnoticed because
  // every current call site wraps the result in as.numeric() regardless,
  // but a caller should not have to defensively coerce a documented
  // vector return into the vector it was documented to already be.
  return List::create(
    Named("baseline")  = Rcpp::NumericVector(baseline.begin(), baseline.end()),
    Named("weights")   = Rcpp::NumericVector(weights.begin(), weights.end()),
    Named("n_iter")    = iter + 1,
    Named("converged") = converged
  );
}
