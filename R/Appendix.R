# ═══════════════════════════════════════════════════════════════
# appendixB.R
#
# Appendix B: Numerical illustration of the sawtooth phenomenon
# for four canonical permutation-based test statistics.
#
#   (a) Two-sample mean difference   [Pitman, 1937]
#   (b) MMD^2  (Gaussian kernel)     [Gretton et al., 2012]
#   (c) HSIC   (Gaussian kernel)     [Gretton et al., 2007]
#   (d) Energy distance              [Székely & Rizzo, 2004]
#
# Requirements:  Rcpp  (install.packages("Rcpp"))
#
# Usage:
#   source("appendixB.R")          # in R / RStudio
#   Rscript appendixB.R            # command line
#
# Outputs:
#   figure_appendixB.pdf   (2 x 2 panel, style matches Figure 1)
# ═══════════════════════════════════════════════════════════════

library(Rcpp)

# ── Locate and compile Rcpp source ────────────────────────────
.this_dir <- tryCatch(dirname(sys.frame(1)$ofile),
                      error = function(e) {
                        args <- commandArgs(trailingOnly = FALSE)
                        fa   <- grep("^--file=", args, value = TRUE)
                        if (length(fa)) dirname(normalizePath(sub("^--file=", "", fa[1])))
                        else getwd()
                      })
.cpp_path <- file.path(.this_dir, "..", "src", "appendixB_core.cpp")
if (!file.exists(.cpp_path)) .cpp_path <- file.path(getwd(), "src", "appendixB_core.cpp")
if (!file.exists(.cpp_path)) .cpp_path <- "appendixB_core.cpp"
cat("Compiling", .cpp_path, "...\n")
sourceCpp(.cpp_path)

set.seed(2026)

# ── Global simulation parameters ─────────────────────────────
alpha <- 0.05
Bmax  <- 500
Nsim  <- 10000

# ═══════════════════════════════════════════════════════════════
#  Helper: Gaussian kernel matrix (median heuristic bandwidth)
# ═══════════════════════════════════════════════════════════════
gauss_kernel <- function(X) {
  # X: n x d matrix (or vector -> treated as n x 1)
  if (is.null(dim(X))) X <- matrix(X, ncol = 1)
  D2 <- as.matrix(dist(X))^2
  sigma2 <- median(D2[upper.tri(D2)])
  if (sigma2 == 0) sigma2 <- 1
  exp(-D2 / (2 * sigma2))
}

center_kernel <- function(K) {
  n <- nrow(K)
  H <- diag(n) - matrix(1/n, n, n)
  H %*% K %*% H
}


# ═══════════════════════════════════════════════════════════════
#  (a) Two-sample mean difference
#      X ~ N(delta, 1)  vs  Y ~ N(0, 1),   n1 = n2 = 20
# ═══════════════════════════════════════════════════════════════
cat("Panel (a): Mean difference ...\n")
n1_a <- 20; n2_a <- 20; delta_a <- 0.55
T_obs_a  <- numeric(Nsim)
T_perm_a <- matrix(0, Nsim, Bmax)

for (j in seq_len(Nsim)) {
  x <- rnorm(n1_a, mean = delta_a)
  y <- rnorm(n2_a, mean = 0)
  z <- c(x, y)
  T_obs_a[j]    <- mean(x) - mean(y)
  T_perm_a[j, ] <- perm_meandiff(z, n1_a, Bmax)
}
pow_a <- power_curve_from_stats(T_obs_a, T_perm_a, alpha)
cat("  done.\n")


# ═══════════════════════════════════════════════════════════════
#  (b) MMD^2  (Gaussian kernel two-sample test)
#      X ~ N(delta*1_d, I_d)  vs  Y ~ N(0, I_d)
#      d = 5,  n1 = n2 = 25
# ═══════════════════════════════════════════════════════════════
cat("Panel (b): MMD^2 ...\n")
n1_b <- 25; n2_b <- 25; d_b <- 5; delta_b <- 0.35
n_b  <- n1_b + n2_b

# Compute observed MMD^2 from kernel matrix
obs_mmd2 <- function(K, n1) {
  n2 <- nrow(K) - n1
  K11 <- K[1:n1, 1:n1];           s11 <- (sum(K11) - sum(diag(K11))) / (n1*(n1-1))
  K22 <- K[(n1+1):nrow(K), (n1+1):nrow(K)]; s22 <- (sum(K22) - sum(diag(K22))) / (n2*(n2-1))
  s12 <- sum(K[1:n1, (n1+1):nrow(K)]) / (n1 * n2)
  s11 + s22 - 2 * s12
}

T_obs_b  <- numeric(Nsim)
T_perm_b <- matrix(0, Nsim, Bmax)

for (j in seq_len(Nsim)) {
  X <- matrix(rnorm(n1_b * d_b, mean = delta_b), n1_b, d_b)
  Y <- matrix(rnorm(n2_b * d_b, mean = 0),       n2_b, d_b)
  Z <- rbind(X, Y)
  K <- gauss_kernel(Z)
  T_obs_b[j]    <- obs_mmd2(K, n1_b)
  T_perm_b[j, ] <- perm_mmd(K, n1_b, Bmax)
}
pow_b <- power_curve_from_stats(T_obs_b, T_perm_b, alpha)
cat("  done.\n")


# ═══════════════════════════════════════════════════════════════
#  (c) HSIC  (Gaussian kernel independence test)
#      X ~ N(0,1),  Y = X^2 + epsilon,  epsilon ~ N(0, sigma^2)
#      n = 50   (nonlinear dependence)
# ═══════════════════════════════════════════════════════════════
cat("Panel (c): HSIC ...\n")
n_c <- 50; sigma_c <- 2.5

obs_hsic <- function(Ktilde, L) {
  n <- nrow(Ktilde)
  sum(Ktilde * L) / (n^2)
}

T_obs_c  <- numeric(Nsim)
T_perm_c <- matrix(0, Nsim, Bmax)

for (j in seq_len(Nsim)) {
  xc <- rnorm(n_c)
  yc <- xc^2 + sigma_c * rnorm(n_c)
  Kx     <- gauss_kernel(xc)
  Ly     <- gauss_kernel(yc)
  Ktilde <- center_kernel(Kx)
  T_obs_c[j]    <- obs_hsic(Ktilde, Ly)
  T_perm_c[j, ] <- perm_hsic(Ktilde, Ly, Bmax)
}
pow_c <- power_curve_from_stats(T_obs_c, T_perm_c, alpha)
cat("  done.\n")


# ═══════════════════════════════════════════════════════════════
#  (d) Energy distance
#      X ~ N(0, sigma_alt^2)  vs  Y ~ N(0, 1)   (scale shift)
#      n1 = n2 = 20
# ═══════════════════════════════════════════════════════════════
cat("Panel (d): Energy distance ...\n")
n1_d <- 20; n2_d <- 20; sigma_d <- 1.8
n_d  <- n1_d + n2_d

obs_energy <- function(D, n1) {
  n2 <- nrow(D) - n1
  D11 <- D[1:n1, 1:n1];           s11 <- (sum(D11) - sum(diag(D11))) / (n1*(n1-1))
  D22 <- D[(n1+1):nrow(D), (n1+1):nrow(D)]; s22 <- (sum(D22) - sum(diag(D22))) / (n2*(n2-1))
  s12 <- sum(D[1:n1, (n1+1):nrow(D)]) / (n1 * n2)
  2*s12 - s11 - s22
}

T_obs_d  <- numeric(Nsim)
T_perm_d <- matrix(0, Nsim, Bmax)

for (j in seq_len(Nsim)) {
  x <- rnorm(n1_d, sd = sigma_d)
  y <- rnorm(n2_d, sd = 1)
  z <- c(x, y)
  D <- as.matrix(dist(z))
  T_obs_d[j]    <- obs_energy(D, n1_d)
  T_perm_d[j, ] <- perm_energy(D, n1_d, Bmax)
}
pow_d <- power_curve_from_stats(T_obs_d, T_perm_d, alpha)
cat("  done.\n")


# ═══════════════════════════════════════════════════════════════
#  Plotting — 2 x 2 panel  (matches Figure 1 style)
# ═══════════════════════════════════════════════════════════════

Bseq    <- 1:Bmax
aligned <- Bseq[abs(alpha * (Bseq + 1) - round(alpha * (Bseq + 1))) < 1e-9]

# Reference power = power at the largest aligned B (proxy for B → ∞)
ref_pow <- function(pow) {
  big <- aligned[aligned <= Bmax]
  pow[max(big)]
}

.fig_dir <- file.path(.this_dir, "..", "figures")
if (!dir.exists(.fig_dir)) dir.create(.fig_dir, recursive = TRUE)
pdf(file.path(.fig_dir, "figure_appendixB.pdf"), width = 10, height = 8.4)
par(mfrow = c(2, 2), mar = c(4.5, 4.8, 2.8, 1.2), family = "serif")

plot_panel <- function(pow, title_expr) {
  rp   <- ref_pow(pow)
  ymin <- max(min(pow) - 0.02, 0)
  ymax <- max(pow) + 0.03
  
  plot(Bseq, pow, type = "l", col = "steelblue", lwd = 1.4,
       xlim = c(1, Bmax), ylim = c(ymin, ymax),
       xaxs = "i", yaxs = "i",
       xlab = expression(italic(B)),
       ylab = "Power",
       main = title_expr,
       cex.main = 1.05, cex.lab = 1.05,
       panel.first = grid(col = "grey90"))
  abline(h = rp, lty = 2, col = "orangered", lwd = 1.4)
  points(aligned, pow[aligned], pch = 16, cex = 0.65, col = "forestgreen")
}

# ── Panel (a) ─────────────────────────────────────────────────
plot_panel(pow_a,
           title_expr = bquote("(a) Mean difference: " *
                                 italic(n)[1] == italic(n)[2] ~ "=" ~ 20 * ", " ~ delta == 0.55))

# ── Panel (b) ─────────────────────────────────────────────────
plot_panel(pow_b,
           title_expr = bquote("(b) " * MMD^2 * ": " *
                                 italic(n)[1] == italic(n)[2] ~ "=" ~ 25 * ", " ~
                                 italic(d) == 5 * ", " ~ delta == 0.35))

# ── Panel (c) ─────────────────────────────────────────────────
plot_panel(pow_c,
           title_expr = bquote("(c) HSIC: " *
                                 italic(Y) == italic(X)^2 + epsilon * ", " ~
                                 italic(n) == 50 * ", " ~ sigma[epsilon] == 2.5))

# ── Panel (d) ─────────────────────────────────────────────────
plot_panel(pow_d,
           title_expr = bquote("(d) Energy: " *
                                 italic(n)[1] == italic(n)[2] ~ "=" ~ 20 * ", " ~
                                 sigma[alt] == 1.8))

# ── Shared legend (bottom-right of last panel) ───────────────
# legend("topright",
#        legend = c(expression(Pow(italic(B))),
#                   expression(widehat(Pow)[infinity]),
#                   expression(alpha(italic(B)+1) %in% scriptstyle(N))),
#        col    = c("steelblue", "orangered", "forestgreen"),
#        lty    = c(1, 2, NA), pch = c(NA, NA, 16),
#        lwd    = c(1.4, 1.4, NA), pt.cex = 0.65,
#        bg = "white", cex = 0.85)

dev.off()
cat("\nSaved: figure_appendixB.pdf\n")