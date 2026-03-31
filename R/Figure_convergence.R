# ═══════════════════════════════════════════════════════════════
# Figure_convergence.R
#
# Extended version of Figure 1: shows that the sawtooth amplitude
# decays as O(B^{-1/2}) and Pow(B) converges to Pow_exact.
#
# Usage:
#   source("R/Figure_convergence.R")       # from repo root
#
# Requirements:
#   install.packages("Rcpp")
# ═══════════════════════════════════════════════════════════════

library(Rcpp)

# ── Locate and compile Rcpp source ────────────────────────────
.this_dir <- tryCatch(dirname(sys.frame(1)$ofile),
                      error = function(e) {
                        args <- commandArgs(trailingOnly = FALSE)
                        fa <- grep("^--file=", args, value = TRUE)
                        if (length(fa)) dirname(normalizePath(sub("^--file=", "", fa[1])))
                        else getwd()
                      })

.cpp_path <- file.path(.this_dir, "..", "src", "pow_mc.cpp")
if (!file.exists(.cpp_path)) .cpp_path <- file.path(getwd(), "src", "pow_mc.cpp")
if (!file.exists(.cpp_path)) .cpp_path <- "pow_mc.cpp"
cat("Compiling", .cpp_path, "...\n")
sourceCpp(.cpp_path)
cat("Done.\n\n")

# ═══════════════════════════════════════════════════════════════
#  CONVERGENCE FIGURE
# ═══════════════════════════════════════════════════════════════

alpha <- 0.05
p1    <- 0.16
Bmax1 <- 5000000     # Case 1: extended range (q(4) ≈ 0.0498 ≈ alpha, slow convergence)
Bmax2 <- 10000     # Case 2: extended range

cat("Computing Case 1 (n=15), Bmax =", Bmax1, "...\n")
pow1 <- compute_power_curve(n = 15, p1 = p1, alpha = alpha, Bmax = Bmax1)
pe1  <- compute_exact_power(n = 15, p1 = p1, alpha = alpha)
cat("  done.\n")

cat("Computing Case 2 (n=25), Bmax =", Bmax2, "...\n")
pow2 <- compute_power_curve(n = 25, p1 = p1, alpha = alpha, Bmax = Bmax2)
pe2  <- compute_exact_power(n = 25, p1 = p1, alpha = alpha)
cat("  done.\n")

Bseq1    <- 1:Bmax1*1.05
Bseq2    <- 1:Bmax2*1.05
aligned1 <- Bseq1[abs(alpha * (Bseq1 + 1) - round(alpha * (Bseq1 + 1))) < 1e-9]
aligned2 <- Bseq2[abs(alpha * (Bseq2 + 1) - round(alpha * (Bseq2 + 1))) < 1e-9]

cat(sprintf("\nCase 1: Pow_exact = %.6f\n", pe1))
cat(sprintf("Case 2: Pow_exact = %.6f\n\n", pe2))

# ── Plot ──────────────────────────────────────────────────────
.fig_dir <- file.path(.this_dir, "..", "figures")
if (!dir.exists(.fig_dir)) dir.create(.fig_dir, recursive = TRUE)

pdf(file.path(.fig_dir, "figure_convergence.pdf"), width = 12, height = 7)

par(mfrow = c(2, 1), mar = c(4, 5, 2.5, 1), family = "serif")

# ─── Case 1 (n = 15, B up to 1M) ─────────────────────────────
plot(Bseq1, pow1, type = "l", col = "steelblue", lwd = 1.5,
     xlim = c(1, Bmax1*1.05), ylim = c(0, 0.25),
     xaxs = "i", yaxs = "i", yaxt = "n",
     xlab = expression(italic(B)), ylab = "Power",
     main = expression("Case 1  (" * italic(n) == 15 * ")"),
     cex.main = 1.8, cex.lab = 1.6, cex.axis = 1.4, font.main = 1,
     panel.first = grid(col = "grey90"))
axis(2, at = seq(0, 0.3, by = 0.1), cex.axis = 1.4)
abline(h = pe1, lty = 2, col = "orangered", lwd = 1.5)

# ─── Case 2 (n = 25, B up to 9999) ──────────────────────────
plot(Bseq2, pow2, type = "l", col = "steelblue", lwd = 1.5,
     xlim = c(1, Bmax2*1.05), ylim = c(0, 0.49),
     xaxs = "i", yaxs = "i", yaxt = "n",
     xlab = expression(italic(B)), ylab = "Power",
     main = expression("Case 2  (" * italic(n) == 25 * ")"),
     cex.main = 1.8, cex.lab = 1.6, cex.axis = 1.4, font.main = 1,
     panel.first = grid(col = "grey90"))
axis(2, at = seq(0, 0.5, by = 0.1), cex.axis = 1.4)
abline(h = pe2, lty = 2, col = "orangered", lwd = 1.5)

dev.off()
cat("Saved: figures/figure_convergence.pdf\n")