// ==============================================================================
// File: pb_pmf.cpp
// Purpose: The Poisson-binomial recursion of cfisher in C++: the probability mass
//   function of a sum of independent Bernoulli(mu_i) variables at 0..kmax, by the
//   forward recursion over subjects
//     f(k) <- q_j f(k) + mu_j f(k - 1),   q_j = 1 - mu_j,
//   run downward in k so that f(k - 1) is still the previous subject's value. The caller
//   passes q_j itself, computed from the log-odds, because 1 - mu_j loses the digits of
//   q_j when mu_j is near 1. Two versions:
//   - pb_scaled_cpp() works in probability scale. The working vector is rescaled against
//     underflow, and the log of the scale is returned so that the caller can restore it.
//     Probabilities far below the largest one still underflow to zero or lose precision.
//   - pb_log_cpp() works in log scale (log-sum-exp in every update), so that no
//     probability underflows, at about 20 times the cost.
//   R/cfisher_null.R (.pb_log()) runs the first and falls back to the second when an
//   underflow bound does not cover the probabilities it needs. Both take O(n * kmax)
//   time and O(kmax) memory.
#include <Rcpp.h>
#include <algorithm>
#include <cmath>
using namespace Rcpp;

// [[Rcpp::export(.pb_scaled_cpp)]]
List pb_scaled_cpp(NumericVector mu, NumericVector q, int kmax) {
  int n = mu.size();
  if (q.size() != n) stop("mu and q must have the same length");
  if (kmax > n) kmax = n;
  if (kmax < 0) stop("kmax must be non-negative");
  NumericVector f(kmax + 1);
  f[0] = 1.0;
  double logscale = 0.0;                      // the pmf equals f * exp(-logscale)
  const double tiny = 1e-20;                  // normalize early: entries far below the running
                                              // maximum feed the later ones
  for (int j = 0; j < n; ++j) {
    double p = mu[j], qj = q[j];
    int top = std::min(j + 1, kmax);          // counts above j + 1 are impossible so far
    for (int k = top; k >= 1; --k) f[k] = qj * f[k] + p * f[k - 1];
    f[0] *= qj;
    // Every 16 subjects, rescale when every retained entry is tiny, so that a long
    // recursion does not underflow. Without truncation, an update keeps at least
    // max(mu_j, q_j) >= 1/2 of the largest entry, so between checks the maximum falls by
    // at most 2^-16. A truncated recursion can fall faster; .pb_log() then detects any
    // lost probability by its underflow bound.
    if ((j & 15) == 15 || j == n - 1) {
      double mx = 0.0;
      for (int k = 0; k <= top; ++k) if (f[k] > mx) mx = f[k];
      if (mx < tiny && mx > 0.0) {
        for (int k = 0; k <= top; ++k) f[k] /= mx;
        logscale -= std::log(mx);
      }
    }
  }
  // The pmf is f * exp(-logscale).
  return List::create(Named("f") = f, Named("logscale") = logscale);
}

// The same recursion in log scale: lf(k) = log f(k), and each update is
//   lf(k) <- log( exp(log(q_j) + lf(k)) + exp(log(mu_j) + lf(k - 1)) ),
// computed as a + log1p(exp(b - a)) with a the larger of the two terms, so that nothing
// overflows or underflows. The caller passes log(mu_j) and log(q_j) themselves, so that
// neither underflows. Impossible counts have lf = -Inf. Returns log f at 0..kmax.
// [[Rcpp::export(.pb_log_cpp)]]
NumericVector pb_log_cpp(NumericVector logmu, NumericVector logq, int kmax) {
  int n = logmu.size();
  if (logq.size() != n) stop("logmu and logq must have the same length");
  if (kmax > n) kmax = n;
  if (kmax < 0) stop("kmax must be non-negative");
  NumericVector lf(kmax + 1, R_NegInf);
  lf[0] = 0.0;                                // log P(empty sum = 0)
  for (int j = 0; j < n; ++j) {
    double lp = logmu[j], lq = logq[j];
    int top = std::min(j + 1, kmax);          // counts above j + 1 are impossible so far
    for (int k = top; k >= 1; --k) {
      double a = lq + lf[k], b = lp + lf[k - 1];
      if (a < b) std::swap(a, b);             // a is the larger term
      // b = -Inf only at k = j + 1, where lf[k] is -Inf and the update is lp + lf[k - 1]
      lf[k] = (b == R_NegInf) ? a : a + std::log1p(std::exp(b - a));
    }
    lf[0] += lq;
  }
  return lf;
}
