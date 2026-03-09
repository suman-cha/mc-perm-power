#include <Rcpp.h>
using namespace Rcpp;

// ═══════════════════════════════════════════════════════════════
// appendixB_core.cpp
//
// Rcpp back-end for Appendix B: sawtooth numerical illustration.
//
// Statistics:
//   (a) Two-sample mean difference   [Pitman, 1937]
//   (b) MMD^2 (Gaussian kernel)      [Gretton et al., 2012]
//   (c) HSIC  (Gaussian kernel)      [Gretton et al., 2007]
//   (d) Energy distance              [Székely & Rizzo, 2004]
//
// Each perm_* function takes pre-computed matrices/vectors for a
// SINGLE dataset realisation and returns Bmax permuted statistics.
//
// power_curve_from_stats() aggregates observed + permuted stats
// across Nsim replications into the unconditional power curve.
// ═══════════════════════════════════════════════════════════════

// ── Fisher–Yates shuffle (R RNG) ─────────────────────────────
inline void shuffle_inplace(IntegerVector& v) {
    int n = v.size();
    for (int i = n - 1; i > 0; i--) {
        int j = (int)(R::runif(0.0, 1.0) * (i + 1));
        if (j > i) j = i;
        std::swap(v[i], v[j]);
    }
}

// ═══════════════════════════════════════════════════════════════
// Power curve from pre-computed statistics
//
//   T_obs  : Nsim-vector   (observed statistic for each rep)
//   T_perm : Nsim x Bmax   (permuted statistics)
//   alpha  : nominal level
//
// For each rep j and B = 1,...,Bmax:
//   R_B = #{b <= B : T_perm(j,b) >= T_obs(j)}
//   k_B = floor(alpha*(B+1)) - 1
//   reject iff R_B <= k_B
//
// Returns Pow(B) = mean rejection rate over replications.
// ═══════════════════════════════════════════════════════════════

// [[Rcpp::export]]
NumericVector power_curve_from_stats(
        NumericVector T_obs,
        NumericMatrix T_perm,
        double alpha) {

    int Nsim = T_obs.size();
    int Bmax = T_perm.ncol();
    NumericVector power(Bmax, 0.0);

    for (int j = 0; j < Nsim; j++) {
        double tobs = T_obs[j];
        int cum_exceed = 0;
        for (int b = 0; b < Bmax; b++) {
            if (T_perm(j, b) >= tobs) cum_exceed++;
            int B  = b + 1;
            int kB = (int)std::floor(alpha * (double)(B + 1)) - 1;
            if (cum_exceed <= kB) {
                power[b] += 1.0;
            }
        }
    }
    for (int b = 0; b < Bmax; b++) power[b] /= (double)Nsim;
    return power;
}


// ═══════════════════════════════════════════════════════════════
// (a) Two-sample mean difference
//     T = mean(z[G1]) - mean(z[G2])
//     Input: pooled z (length n), first n1 = group 1.
// ═══════════════════════════════════════════════════════════════

// [[Rcpp::export]]
NumericVector perm_meandiff(NumericVector z, int n1, int Bmax) {
    int n  = z.size();
    int n2 = n - n1;
    NumericVector stats(Bmax);
    IntegerVector idx = seq_len(n) - 1;          // 0-based

    for (int b = 0; b < Bmax; b++) {
        IntegerVector pidx = clone(idx);
        shuffle_inplace(pidx);
        double s1 = 0.0, s2 = 0.0;
        for (int i = 0; i < n1; i++) s1 += z[pidx[i]];
        for (int i = n1; i < n;  i++) s2 += z[pidx[i]];
        stats[b] = s1 / (double)n1 - s2 / (double)n2;
    }
    return stats;
}


// ═══════════════════════════════════════════════════════════════
// (b) MMD^2  (U-statistic form)
//     Input: n x n kernel matrix K, first n1 = group 1.
//
//     MMD^2 = 2/(n1(n1-1)) sum_{i<j in G1} K_{ij}
//           + 2/(n2(n2-1)) sum_{i<j in G2} K_{ij}
//           - 2/(n1*n2)    sum_{i in G1, j in G2} K_{ij}
// ═══════════════════════════════════════════════════════════════

// [[Rcpp::export]]
NumericVector perm_mmd(NumericMatrix K, int n1, int Bmax) {
    int n  = K.nrow();
    int n2 = n - n1;
    double d11 = (double)n1 * (n1 - 1);
    double d22 = (double)n2 * (n2 - 1);
    double d12 = (double)n1 * n2;
    NumericVector stats(Bmax);
    IntegerVector idx = seq_len(n) - 1;

    for (int b = 0; b < Bmax; b++) {
        IntegerVector pidx = clone(idx);
        shuffle_inplace(pidx);

        double s11 = 0.0, s22 = 0.0, s12 = 0.0;
        for (int i = 0; i < n; i++) {
            int pi_i = pidx[i];
            bool g1_i = (i < n1);
            for (int j = i + 1; j < n; j++) {
                int pi_j = pidx[j];
                bool g1_j = (j < n1);
                double kij = K(pi_i, pi_j);
                if      ( g1_i &&  g1_j) s11 += kij;
                else if (!g1_i && !g1_j) s22 += kij;
                else                     s12 += kij;
            }
        }
        stats[b] = 2.0 * s11 / d11 + 2.0 * s22 / d22 - 2.0 * s12 / d12;
    }
    return stats;
}


// ═══════════════════════════════════════════════════════════════
// (c) HSIC  (biased V-statistic)
//     Ktilde = H K_x H   (pre-centered X-kernel, computed in R)
//     L      = K_y        (raw Y-kernel)
//
//     Under permutation pi of Y-labels:
//       HSIC_pi = (1/n^2) sum_{i,j} Ktilde_{ij} * L_{pi(i),pi(j)}
// ═══════════════════════════════════════════════════════════════

// [[Rcpp::export]]
NumericVector perm_hsic(NumericMatrix Ktilde, NumericMatrix L,
                        int Bmax) {
    int n = Ktilde.nrow();
    double nn = (double)n * n;
    NumericVector stats(Bmax);
    IntegerVector idx = seq_len(n) - 1;

    for (int b = 0; b < Bmax; b++) {
        IntegerVector pidx = clone(idx);
        shuffle_inplace(pidx);

        double h = 0.0;
        for (int i = 0; i < n; i++) {
            int pi_i = pidx[i];
            for (int j = 0; j < n; j++) {
                h += Ktilde(i, j) * L(pi_i, pidx[j]);
            }
        }
        stats[b] = h / nn;
    }
    return stats;
}


// ═══════════════════════════════════════════════════════════════
// (d) Energy distance
//     Input: n x n Euclidean distance matrix D, first n1 = group 1.
//
//     E = 2/(n1*n2)    sum_{cross} D_{ij}
//       - 2/(n1(n1-1)) sum_{within G1} D_{ij}
//       - 2/(n2(n2-1)) sum_{within G2} D_{ij}
// ═══════════════════════════════════════════════════════════════

// [[Rcpp::export]]
NumericVector perm_energy(NumericMatrix D, int n1, int Bmax) {
    int n  = D.nrow();
    int n2 = n - n1;
    double d11 = (double)n1 * (n1 - 1);
    double d22 = (double)n2 * (n2 - 1);
    double d12 = (double)n1 * n2;
    NumericVector stats(Bmax);
    IntegerVector idx = seq_len(n) - 1;

    for (int b = 0; b < Bmax; b++) {
        IntegerVector pidx = clone(idx);
        shuffle_inplace(pidx);

        double s11 = 0.0, s22 = 0.0, s12 = 0.0;
        for (int i = 0; i < n; i++) {
            int pi_i = pidx[i];
            bool g1_i = (i < n1);
            for (int j = i + 1; j < n; j++) {
                int pi_j = pidx[j];
                bool g1_j = (j < n1);
                double dij = D(pi_i, pi_j);
                if      ( g1_i &&  g1_j) s11 += dij;
                else if (!g1_i && !g1_j) s22 += dij;
                else                     s12 += dij;
            }
        }
        stats[b] = 2.0 * s12 / d12 - 2.0 * s11 / d11 - 2.0 * s22 / d22;
    }
    return stats;
}
