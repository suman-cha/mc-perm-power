# ═══════════════════════════════════════════════════════════════
# Figure_Bnew.R
#
# "More Permutations Do Not Always Increase Power
#  Non-monotonicity in Monte Carlo Permutation Tests"
#
# Illustration figure for the equation B_new section.
# Single panel with exactly the visual style and parameters of
# Figure1.R (Case 2). Filled green dots mark the strict local
# maxima alpha*(B+1) in N (the values apply_Bnew recommends);
# open green circles mark plateau values that the modification
# rule recommends moving away from.
#
# Usage
#   source("R/Figure_Bnew.R")
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

# ── 2. apply_Bnew rule (one-line modification, R side) ────────
# A value B is a "jump" iff k_B > k_{B-1}. If the user-specified
# B is a jump, return it unchanged. Otherwise, step back to the
# largest jump value strictly smaller than B. This always returns
# a strict local maximizer of the unconditional power curve.
is_jump <- function(B, alpha) {
  if (B <= 0) return(FALSE)
  kB   <- floor(alpha * (B + 1)) - 1
  kBm1 <- floor(alpha *  B     ) - 1
  kB > kBm1
}

apply_Bnew <- function(B, alpha) {
  if (is_jump(B, alpha)) return(B)
  if (B <= 1) return(B)
  for (Bp in seq.int(B - 1L, 1L, by = -1L)) {
    if (is_jump(Bp, alpha)) return(Bp)
  }
  B
}

# Unit checks
stopifnot(apply_Bnew(50,  0.05) == 39)
stopifnot(apply_Bnew(78,  0.05) == 59)
stopifnot(apply_Bnew(19,  0.05) == 19)
stopifnot(apply_Bnew(50,  0.10) == 49)
stopifnot(apply_Bnew(29,  0.10) == 29)
stopifnot(apply_Bnew(30,  0.90) == 29)
stopifnot(apply_Bnew(31,  0.90) == 31)
cat("[unit] apply_Bnew passed all checks\n")

# ── 3. Panel parameters (identical to Figure1.R Case 2) ───────
alpha <- 0.05
p1    <- 0.16
n     <- 25
Bmax  <- 150

# ── 4. Compute power curve and Pow_exact ──────────────────────
cat(sprintf("Computing power curve (alpha = %.2f, n = %d, p1 = %.2f, Bmax = %d) ...\n",
            alpha, n, p1, Bmax))
pow <- compute_power_curve(n = n, p1 = p1, alpha = alpha, Bmax = Bmax)
pe  <- compute_exact_power (n = n, p1 = p1, alpha = alpha)
cat(sprintf("  Pow_exact = %.4f\n\n", pe))

# Aligned B values where alpha * (B+1) is an integer
Bseq    <- 1:Bmax
aligned <- Bseq[abs(alpha * (Bseq + 1) - round(alpha * (Bseq + 1))) < 1e-9]

# Modified curve: Pow_mod(B) = Pow(apply_Bnew(B, alpha))
# This is the curve obtained after applying the equation B_new
# modification, which replaces the descending plateau by a constant
# value equal to the previous strict local maximum.
pow_mod <- numeric(Bmax)
for (B in 1:Bmax) {
  pow_mod[B] <- pow[apply_Bnew(B, alpha)]
}

# ── 5. Verification asserts (monotonicity of apply_Bnew) ──────
violations <- 0L
for (B in 1:Bmax) {
  Bp <- apply_Bnew(B, alpha)
  if (is_jump(B, alpha)) {
    if (Bp != B) violations <- violations + 1L
  } else {
    if (pow[Bp] + 1e-12 < pow[B]) violations <- violations + 1L
  }
}
cat(sprintf("[verify] monotonicity violations = %d / %d\n\n", violations, Bmax))
stopifnot(violations == 0L)

# ── 6. Render and save PDF (1 row x 2 columns: before vs after) ─
.fig_dir <- file.path(.this_script_dir, "..", "figures")
if (!dir.exists(.fig_dir)) dir.create(.fig_dir, recursive = TRUE)

# Open-circle subset (about 5 per sawtooth segment, excluding aligned neighbourhood)
.seg_len <- if (length(aligned) >= 2) aligned[2] - aligned[1] else 20
.by      <- max(1, round(.seg_len / 5))
.near    <- unique(unlist(lapply(aligned, function(a) (a - 1):(a + 1))))
.scatter_idx <- setdiff(seq(1, Bmax, by = .by), .near)

draw_Bnew <- function(out_path, pdf_w, pdf_h) {
  pdf(out_path, width = pdf_w, height = pdf_h)
  par(mfrow = c(1, 2), mar = c(5, 5.5, 3, 1.5), family = "serif")

  draw_panel <- function(y_data, main_expr) {
    plot(Bseq, y_data, type = "l", col = "steelblue", lwd = 1.5,
         xlim = c(0, Bmax), ylim = c(0, 0.55),
         xaxs = "i", yaxs = "i", xaxt = "n", yaxt = "n",
         xlab = expression(italic(B)), ylab = "Power",
         main = main_expr,
         cex.main = 1.6, cex.lab = 1.6, cex.axis = 1.4, font.main = 1,
         panel.first = grid(col = "grey90"))
    axis(1, at = seq(0, Bmax, by = 50), cex.axis = 1.4)
    axis(2, at = seq(0, 0.5, by = 0.1),
         labels = sprintf("%.1f", seq(0, 0.5, by = 0.1)),
         cex.axis = 1.2, las = 1)
    abline(h = pe, lty = 2, col = "orangered", lwd = 1.5)
    points(.scatter_idx, y_data[.scatter_idx], pch = 1, cex = 0.9,
           col = adjustcolor("forestgreen", alpha.f = 0.6))
    points(aligned, y_data[aligned], pch = 16, cex = 1.1, col = "forestgreen")
  }

  # Left  raw Pow(B) sawtooth
  draw_panel(pow, main_expr = expression(Pow(italic(B))))
  # Right  Pow(B_new) staircase
  draw_panel(pow_mod, main_expr = expression(Pow(italic(B)[new])))

  dev.off()
  cat(sprintf("Saved %s\n", basename(out_path)))
}

# 2-column version (current size)
draw_Bnew(
  out_path = file.path(.fig_dir, "figure_Bnew_2col.pdf"),
  pdf_w = 10, pdf_h = 4
)

# 1-column version
# Larger PDF dimensions intentional. LaTeX scales the PDF down to
# \columnwidth (~3.15 in), so a larger source PDF means a smaller
# scaling factor and proportionally smaller text/marks in print.
# Keeping the same cex as 2-col preserves panel-relative proportions.
draw_Bnew(
  out_path = file.path(.fig_dir, "figure_Bnew_1col.pdf"),
  pdf_w = 12, pdf_h = 5
)

# ── 7. Console log of apply_Bnew on representative B values ──
cat("\n[examples] apply_Bnew recommendations\n")
for (B in c(20, 35, 60, 100, 130)) {
  Bp <- apply_Bnew(B, alpha)
  delta <- pow[Bp] - pow[B]
  tag <- if (Bp == B) "(jump, keep)" else sprintf("-> %d  gain = %+0.4f", Bp, delta)
  cat(sprintf("  B = %3d  Pow = %.4f  %s\n", B, pow[B], tag))
}
