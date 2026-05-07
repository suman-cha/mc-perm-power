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

# ── Subsample for plotting (keeps PDF small + fast to load) ──
# 5,000,000 points per panel produces a huge PDF (~2 MB) that loads
# slowly in Overleaf. We subsample to N_PLOT evenly spaced points,
# which is visually indistinguishable since the line is continuous.
N_PLOT <- 4000
sub1 <- unique(round(seq(1, length(Bseq1), length.out = N_PLOT)))
sub2 <- unique(round(seq(1, length(Bseq2), length.out = N_PLOT)))
Bseq1_p <- Bseq1[sub1];  pow1_p <- pow1[sub1]
Bseq2_p <- Bseq2[sub2];  pow2_p <- pow2[sub2]
cat(sprintf("Plot subsample  Case 1 %d -> %d points,  Case 2 %d -> %d points\n",
            length(Bseq1), length(sub1),
            length(Bseq2), length(sub2)))

# ── Plot helper ───────────────────────────────────────────────
.fig_dir <- file.path(.this_dir, "..", "figures")
if (!dir.exists(.fig_dir)) dir.create(.fig_dir, recursive = TRUE)

draw_convergence <- function(out_path, pdf_w, pdf_h,
                             mar_vec, cex_main, cex_lab, cex_axis,
                             lwd_line, lwd_ref) {
  pdf(out_path, width = pdf_w, height = pdf_h)
  par(mfrow = c(2, 1), mar = mar_vec, family = "serif")

  # Case 1
  plot(Bseq1_p, pow1_p, type = "l", col = "steelblue", lwd = lwd_line,
       xlim = c(1, Bmax1 * 1.05), ylim = c(0, 0.25),
       xaxs = "i", yaxs = "i", yaxt = "n",
       xlab = expression(italic(B)), ylab = "Power",
       main = expression("Case 1  (" * italic(n) == 15 * ")"),
       cex.main = cex_main, cex.lab = cex_lab, cex.axis = cex_axis,
       font.main = 1, panel.first = grid(col = "grey90"))
  axis(2, at = seq(0, 0.3, by = 0.1), cex.axis = cex_axis)
  abline(h = pe1, lty = 2, col = "orangered", lwd = lwd_ref)

  # Case 2
  plot(Bseq2_p, pow2_p, type = "l", col = "steelblue", lwd = lwd_line,
       xlim = c(1, Bmax2 * 1.05), ylim = c(0, 0.49),
       xaxs = "i", yaxs = "i", yaxt = "n",
       xlab = expression(italic(B)), ylab = "Power",
       main = expression("Case 2  (" * italic(n) == 25 * ")"),
       cex.main = cex_main, cex.lab = cex_lab, cex.axis = cex_axis,
       font.main = 1, panel.first = grid(col = "grey90"))
  axis(2, at = seq(0, 0.5, by = 0.1), cex.axis = cex_axis)
  abline(h = pe2, lty = 2, col = "orangered", lwd = lwd_ref)

  dev.off()
  cat(sprintf("Saved %s\n", basename(out_path)))
}

# 2-column version (current size, keeps current cex)
draw_convergence(
  out_path = file.path(.fig_dir, "figure_convergence_2col.pdf"),
  pdf_w = 8, pdf_h = 5,
  mar_vec = c(4, 5, 2.5, 1),
  cex_main = 1.6, cex_lab = 1.4, cex_axis = 1.2,
  lwd_line = 1.5, lwd_ref = 1.5
)

# 1-column version
# Larger PDF dimensions intentional. LaTeX scales the PDF down to
# \columnwidth (~3.15 in), so a larger source PDF means a smaller
# scaling factor and proportionally smaller text/marks in print.
# Keeping the same cex as 2-col preserves panel-relative proportions.
draw_convergence(
  out_path = file.path(.fig_dir, "figure_convergence_1col.pdf"),
  pdf_w = 10, pdf_h = 6,
  mar_vec = c(4, 5, 2.5, 1),
  cex_main = 1.6, cex_lab = 1.4, cex_axis = 1.2,
  lwd_line = 1.5, lwd_ref = 1.5
)