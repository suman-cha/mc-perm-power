// pow_mc.cpp
// ─────────────────────────────────────────────────────────────
// Fast computation of Monte Carlo permutation test power curves
// from Equation (4) of "More Permutations Do Not Always
// Increase Power: Non-monotonicity in MC Permutation Tests"
//
// Compile in R:  Rcpp::sourceCpp("pow_mc.cpp")
// ─────────────────────────────────────────────────────────────

#include <Rcpp.h>
#include <cmath>
#include <vector>
#include <algorithm>
using namespace Rcpp;

// ── log C(n,k) via lgamma ────────────────────────────────────
inline double log_binom(int n, int k) {
    if (k < 0 || k > n) return R_NegInf;
    return lgamma(n + 1.0) - lgamma(k + 1.0) - lgamma(n - k + 1.0);
}

// ── Binomial PMF (log-space for numerical stability) ─────────
inline double binom_pmf(int k, int n, double q) {
    if (k < 0 || k > n) return 0.0;
    if (q <= 0.0) return (k == 0) ? 1.0 : 0.0;
    if (q >= 1.0) return (k == n) ? 1.0 : 0.0;
    return exp(log_binom(n, k) + k * log(q) + (n - k) * log(1.0 - q));
}

// ── Binomial CDF: P(Bin(n,q) <= k) ──────────────────────────
inline double binom_cdf(int k, int n, double q) {
    if (k < 0)  return 0.0;
    if (k >= n) return 1.0;
    // Use R's pbinom for accuracy (regularized incomplete beta)
    return R::pbinom((double)k, (double)n, q, 1, 0);
}

// ─────────────────────────────────────────────────────────────
// q(s) = C(2n-s, n-s) / C(2n, n)
// Conditional exceedance probability given S = s successes
// ─────────────────────────────────────────────────────────────
// [[Rcpp::export]]
NumericVector compute_q(int n) {
    double log_den = log_binom(2 * n, n);
    NumericVector q(n + 1);
    for (int s = 0; s <= n; s++) {
        q[s] = exp(log_binom(2 * n - s, n - s) - log_den);
    }
    return q;
}

// ─────────────────────────────────────────────────────────────
// k_B = floor(alpha * (B+1)) - 1
// Integer critical count
// ─────────────────────────────────────────────────────────────
// [[Rcpp::export]]
IntegerVector compute_kB(double alpha, int Bmax) {
    IntegerVector kB(Bmax);
    for (int B = 1; B <= Bmax; B++) {
        kB[B - 1] = (int)floor(alpha * (B + 1)) - 1;
    }
    return kB;
}

// ─────────────────────────────────────────────────────────────
// Pow(B) = Σ_s P(S=s) · P(Bin(B, q(s)) <= k_B)
//
// Equation (4) — the full unconditional power
// ─────────────────────────────────────────────────────────────
// [[Rcpp::export]]
NumericVector compute_power_curve(int n, double p1, double alpha, int Bmax) {
    // q(s)
    double log_den = log_binom(2 * n, n);
    std::vector<double> q(n + 1);
    for (int s = 0; s <= n; s++) {
        q[s] = exp(log_binom(2 * n - s, n - s) - log_den);
    }

    // P(S = s), S ~ Bin(n, p1)
    std::vector<double> pS(n + 1);
    for (int s = 0; s <= n; s++) {
        pS[s] = R::dbinom((double)s, (double)n, p1, 0);
    }

    // Pow(B) for B = 1 ... Bmax
    NumericVector pow(Bmax);
    for (int B = 1; B <= Bmax; B++) {
        int kB = (int)floor(alpha * (B + 1)) - 1;
        double power = 0.0;
        for (int s = 0; s <= n; s++) {
            power += pS[s] * binom_cdf(kB, B, q[s]);
        }
        pow[B - 1] = power;
    }
    return pow;
}

// ─────────────────────────────────────────────────────────────
// Exact permutation test power
// Reject iff q(S) < alpha, i.e. S >= s_min where
// s_min = min{s : q(s) < alpha}
// ─────────────────────────────────────────────────────────────
// [[Rcpp::export]]
double compute_exact_power(int n, double p1, double alpha) {
    double log_den = log_binom(2 * n, n);

    // Find threshold
    int s_min = n + 1;
    for (int s = 0; s <= n; s++) {
        double qs = exp(log_binom(2 * n - s, n - s) - log_den);
        if (qs < alpha) { s_min = s; break; }
    }

    // Pow_exact = P(S >= s_min) = 1 - P(S <= s_min - 1)
    if (s_min > n) return 0.0;
    return R::pbinom((double)(s_min - 1), (double)n, p1, 0, 0);  // lower.tail=FALSE
}

// ─────────────────────────────────────────────────────────────
// Proposition 3.3 upper bound on one-step power decrease
//   Pow(B) - Pow(B+1) <= 1/sqrt(2*pi*(B+1)) * sqrt(alpha/(1-alpha))
// ─────────────────────────────────────────────────────────────
// [[Rcpp::export]]
NumericVector compute_prop33_bound(double alpha, int Bmax) {
    NumericVector bound(Bmax);
    double coeff = sqrt(alpha / (1.0 - alpha));
    for (int B = 1; B <= Bmax; B++) {
        bound[B - 1] = coeff / sqrt(2.0 * M_PI * (B + 1));
    }
    return bound;
}

// ─────────────────────────────────────────────────────────────
// Conditional rejection probability for a single dataset
//   phi_B(q) = P(Bin(B, q) <= k_B)
// Useful for visualizing the conditional power surface
// ─────────────────────────────────────────────────────────────
// [[Rcpp::export]]
NumericMatrix compute_conditional_power(NumericVector q_grid,
                                         double alpha, int Bmax) {
    int nq = q_grid.size();
    NumericMatrix phi(Bmax, nq);
    for (int B = 1; B <= Bmax; B++) {
        int kB = (int)floor(alpha * (B + 1)) - 1;
        for (int j = 0; j < nq; j++) {
            phi(B - 1, j) = binom_cdf(kB, B, q_grid[j]);
        }
    }
    return phi;
}
