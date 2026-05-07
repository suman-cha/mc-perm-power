# ═══════════════════════════════════════════════════════════════
# Figure1_clean.R
#
# Re-draws Figure 1 (figure1_bmax150.pdf) using the same code
# pattern as Figure_Bnew.R: a single draw_panel() helper that
# always overlays both filled green dots (at aligned B values
# where alpha (B+1) is integer) and open green circles (scatter
# subset on the descending plateau between aligned values).
#
# Differs from Figure1.R in that the plotting code is factored
# into a draw_panel() helper and the both-marker overlay is
# applied unconditionally for any Bmax (whereas Figure1.R
# suppresses the open circles when Bmax = 400).
#
# Usage
#   source("R/Figure1_clean.R")
#   reproduce_figure1_clean(Bmax = 150, save_pdf = TRUE)
# ═══════════════════════════════════════════════════════════════

library(Rcpp)

# ── 1. Locate and compile pow_mc.cpp ──────────────────────────
.this_script_dir <- tryCatch({
  dirname(sys.frame(1)$ofile)
}, error = function(e) {
  if (requireNamespace("rstudioapi", quietly = TRUE) &&
      rstudioapi::isAvailable()) {
    dirname(rstudioapi::getSourceEditorContext()$path)
  } else {
    args <- commandArgs(trailingOnly = FALSE)
    file_arg <- grep("^--file=", args, value = TRUE)
    if (length(file_arg) > 0) {
      dirname(normalizePath(sub("^--file=", "", file_arg[1])))
    } else {
      getwd()
    }
  }
})

.cpp_path <- file.path(.this_script_dir, "..", "src", "pow_mc.cpp")
if (!file.exists(.cpp_path)) {
  candidates <- c(
    file.path(.this_script_dir, "pow_mc.cpp"),
    file.path(getwd(), "src", "pow_mc.cpp"),
    file.path(getwd(), "pow_mc.cpp")
  )
  found <- candidates[file.exists(candidates)]
  if (length(found) > 0) {
    .cpp_path <- found[1]
  } else {
    stop("Cannot find pow_mc.cpp. Run from repository root.")
  }
}
cat(sprintf("Compiling Rcpp code from %s ...\n", .cpp_path))
sourceCpp(.cpp_path)
cat("Done.\n\n")

# ── 2. Reproduction function (Figure_Bnew.R-style draw_panel) ──
reproduce_figure1_clean <- function(Bmax = 150,
                                    out_path,
                                    pdf_width, pdf_height) {
  alpha <- 0.05
  p1    <- 0.16

  # Compute power curves from Eq. (4)
  cat(sprintf("Computing Case 1 (n = 15, Bmax = %d) ...\n", Bmax))
  pow1 <- compute_power_curve(n = 15, p1 = p1, alpha = alpha, Bmax = Bmax)
  pe1  <- compute_exact_power(n = 15, p1 = p1, alpha = alpha)

  cat(sprintf("Computing Case 2 (n = 25, Bmax = %d) ...\n", Bmax))
  pow2 <- compute_power_curve(n = 25, p1 = p1, alpha = alpha, Bmax = Bmax)
  pe2  <- compute_exact_power(n = 25, p1 = p1, alpha = alpha)

  Bseq <- 1:Bmax
  aligned <- Bseq[abs(alpha * (Bseq + 1) - round(alpha * (Bseq + 1))) < 1e-9]

  cat(sprintf("\n  Case 1  Pow_exact = %.6f\n", pe1))
  cat(sprintf("  Case 2  Pow_exact = %.6f\n\n", pe2))

  pdf(out_path, width = pdf_width, height = pdf_height)
  par(mfrow = c(1, 2), mar = c(5, 5.5, 3, 1.5), family = "serif")

  # Open-circle subset (about 5 per sawtooth segment, excluding aligned neighbourhood)
  .seg_len <- if (length(aligned) >= 2) aligned[2] - aligned[1] else 20
  .by      <- max(1, round(.seg_len / 5))
  .near    <- unique(unlist(lapply(aligned, function(a) (a - 1):(a + 1))))
  .scatter_idx <- setdiff(seq(1, Bmax, by = .by), .near)

  draw_panel <- function(y_data, pe_val, main_expr) {
    plot(Bseq, y_data, type = "l", col = "steelblue", lwd = 1.5,
         xlim = c(0, Bmax), ylim = c(0, 0.55),
         xaxs = "i", yaxs = "i", xaxt = "n", yaxt = "n",
         xlab = expression(italic(B)), ylab = "Power",
         main = main_expr,
         cex.main = 1.8, cex.lab = 1.6, cex.axis = 1.4, font.main = 1,
         panel.first = grid(col = "grey90"))
    axis(1, at = seq(0, Bmax, by = 50), cex.axis = 1.4)
    axis(2, at = seq(0, 0.5, by = 0.1),
         labels = sprintf("%.1f", seq(0, 0.5, by = 0.1)),
         cex.axis = 1.2, las = 1)
    abline(h = pe_val, lty = 2, col = "orangered", lwd = 1.5)
    points(.scatter_idx, y_data[.scatter_idx], pch = 1, cex = 0.9,
           col = adjustcolor("forestgreen", alpha.f = 0.6))
    points(aligned, y_data[aligned], pch = 16, cex = 1.1, col = "forestgreen")
  }

  draw_panel(pow1, pe1,
             main_expr = expression("Case 1  (" * italic(n) == 15 * ")"))
  draw_panel(pow2, pe2,
             main_expr = expression("Case 2  (" * italic(n) == 25 * ")"))

  dev.off()
  cat(sprintf("Saved %s\n", basename(out_path)))
}


# ═══════════════════════════════════════════════════════════════
#  RUN  generate both 2-col and 1-col versions at Bmax = 150
# ═══════════════════════════════════════════════════════════════

cat("══════════════════════════════════════════════════\n")
cat("  Reproducing figure1 (Figure_Bnew.R style) at Bmax = 150 ...\n")
cat("══════════════════════════════════════════════════\n\n")

.fig_dir <- file.path(.this_script_dir, "..", "figures")
if (!dir.exists(.fig_dir)) dir.create(.fig_dir, recursive = TRUE)

# 2-column version (current size)
reproduce_figure1_clean(
  Bmax = 150,
  out_path = file.path(.fig_dir, "figure1_bmax150_2col.pdf"),
  pdf_width = 10, pdf_height = 4
)

# 1-column version
# Larger PDF dimensions intentional. LaTeX scales the PDF down to
# \columnwidth (~3.15 in), so a larger source PDF means a smaller
# scaling factor and proportionally smaller text/marks in print.
# Keeping the same cex as 2-col preserves panel-relative proportions.
reproduce_figure1_clean(
  Bmax = 150,
  out_path = file.path(.fig_dir, "figure1_bmax150_1col.pdf"),
  pdf_width = 12, pdf_height = 5
)

cat("\nAvailable function\n")
cat("  reproduce_figure1_clean(Bmax, out_path, pdf_width, pdf_height)\n")
