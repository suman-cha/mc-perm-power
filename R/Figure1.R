# ═══════════════════════════════════════════════════════════════
# figure1_reproduce.R
#
# "More Permutations Do Not Always Increase Power:
#  Non-monotonicity in Monte Carlo Permutation Tests"
#
# Reproduces Figure 1 and provides tools for experimentation.
# All computations use Equation (4) via Rcpp.
#
# Usage:
#   source("figure1_reproduce.R")
#
# Requirements:
#   install.packages("Rcpp")  # if not already installed
# ═══════════════════════════════════════════════════════════════

library(Rcpp)

# ── 1. Compile Rcpp code ──────────────────────────────────────
# Locate pow_mc.cpp relative to this script's directory,
# so it works regardless of getwd().
cat("Compiling Rcpp code...\n")

.this_script_dir <- tryCatch({
  # Works when called via source()
  dirname(sys.frame(1)$ofile)
}, error = function(e) {
  # Fallback: try rstudioapi (RStudio), commandArgs (Rscript),
  # or finally just use getwd()
  if (requireNamespace("rstudioapi", quietly = TRUE) &&
      rstudioapi::isAvailable()) {
    dirname(rstudioapi::getSourceEditorContext()$path)
  } else {
    # Rscript --file=... 
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
  # Last resort: search in common locations
  candidates <- c(
    file.path(.this_script_dir, "pow_mc.cpp"),
    file.path(getwd(), "src", "pow_mc.cpp"),
    file.path(getwd(), "pow_mc.cpp")
  )
  found <- candidates[file.exists(candidates)]
  if (length(found) > 0) {
    .cpp_path <- found[1]
  } else {
    stop(paste0(
      "Cannot find 'pow_mc.cpp'.\n",
      "  Searched:\n",
      "    - ", file.path(.this_script_dir, "..", "src", "pow_mc.cpp"), "\n",
      "    - ", file.path(getwd(), "src", "pow_mc.cpp"), "\n\n",
      "  Please run from the repository root:  source('R/Figure1.R')"
    ))
  }
}

cat(sprintf("  Found: %s\n", .cpp_path))
sourceCpp(.cpp_path)
cat("Done.\n\n")

# ═══════════════════════════════════════════════════════════════
#  FIGURE 1 REPRODUCTION
# ═══════════════════════════════════════════════════════════════

reproduce_figure1 <- function(Bmax = 400, save_pdf = FALSE) {
  alpha <- 0.05
  p1    <- 0.16
  
  # ── Compute power curves from Eq. (4) ─────────────────────
  cat("Computing Case 1 (n=15)...\n")
  pow1 <- compute_power_curve(n = 15, p1 = p1, alpha = alpha, Bmax = Bmax)
  pe1  <- compute_exact_power(n = 15, p1 = p1, alpha = alpha)
  
  cat("Computing Case 2 (n=25)...\n")
  pow2 <- compute_power_curve(n = 25, p1 = p1, alpha = alpha, Bmax = Bmax)
  pe2  <- compute_exact_power(n = 25, p1 = p1, alpha = alpha)
  
  # ── Aligned B: alpha(B+1) in N ────────────────────────────
  Bseq <- 1:Bmax
  aligned <- Bseq[abs(alpha * (Bseq + 1) - round(alpha * (Bseq + 1))) < 1e-9]
  
  cat(sprintf("\nCase 1: Pow_exact = %.6f\n", pe1))
  cat(sprintf("Case 2: Pow_exact = %.6f\n", pe2))
  cat(sprintf("Local maxima at B = %s, ...\n\n",
              paste(head(aligned, 8), collapse = ", ")))
  
  # ── Plot ──────────────────────────────────────────────────
  .fig_dir <- file.path(.this_script_dir, "..", "figures")
  if (!dir.exists(.fig_dir)) dir.create(.fig_dir, recursive = TRUE)
  if (save_pdf) pdf(file.path(.fig_dir, "figure1.pdf"), width = 10, height = 4.2)
  
  par(mfrow = c(1, 2), mar = c(4.5, 4.5, 2.5, 1), family = "serif")
  
  # ─── Case 1 ───
  plot(Bseq, pow1, type = "l", col = "steelblue", lwd = 1.5,
       xlim = c(1, Bmax), ylim = c(0, 0.55),
       xaxs = "i", yaxs = "i",
       xlab = expression(italic(B)), ylab = "Power",
       main = expression("Case 1: " * italic(n) == 15),
       panel.first = grid(col = "grey90"))
  abline(h = pe1, lty = 2, col = "orangered", lwd = 1.5)
  points(aligned, pow1[aligned], pch = 16, cex = 0.8, col = "forestgreen")
  
  # ─── Case 2 ───
  plot(Bseq, pow2, type = "l", col = "steelblue", lwd = 1.5,
       xlim = c(1, Bmax), ylim = c(0, 0.55),
       xaxs = "i", yaxs = "i",
       xlab = expression(italic(B)), ylab = "Power",
       main = expression("Case 2: " * italic(n) == 25),
       panel.first = grid(col = "grey90"))
  abline(h = pe2, lty = 2, col = "orangered", lwd = 1.5)
  points(aligned, pow2[aligned], pch = 16, cex = 0.8, col = "forestgreen")
  
  # ─── Legend (on Case 2 panel) ───
  # legend("topright",
  #        legend = c(expression(Pow(italic(B))),
  #                   expression(Pow[exact]),
  #                   expression(alpha %.% (italic(B)+1) %in% italic(N))),
  #        col = c("steelblue", "orangered", "forestgreen"),
  #        lty = c(1, 2, NA), pch = c(NA, NA, 16),
  #        lwd = c(1.5, 1.5, NA), pt.cex = 0.8,
  #        bg = "white", cex = 0.85)
  # 
  if (save_pdf) {
    dev.off()
    cat("Saved to figure1.pdf\n")
  }
  
  invisible(list(pow1 = pow1, pow2 = pow2, pe1 = pe1, pe2 = pe2,
                 aligned = aligned))
}


# ═══════════════════════════════════════════════════════════════
#  GENERAL EXPERIMENTATION FUNCTIONS
# ═══════════════════════════════════════════════════════════════

# ── Compute and plot for arbitrary (n, p1, alpha) ─────────────
run_experiment <- function(n, p1, alpha = 0.05, Bmax = 400,
                           main_title = NULL) {
  
  pow <- compute_power_curve(n, p1, alpha, Bmax)
  pe  <- compute_exact_power(n, p1, alpha)
  q   <- compute_q(n)
  kB  <- compute_kB(alpha, Bmax)
  
  Bseq    <- 1:Bmax
  aligned <- Bseq[abs(alpha * (Bseq + 1) - round(alpha * (Bseq + 1))) < 1e-9]
  
  # ── Local maxima / minima detection ───────────────────────
  local_max <- c()
  local_min <- c()
  for (i in 2:(Bmax - 1)) {
    if (pow[i] > pow[i-1] && pow[i] > pow[i+1]) local_max <- c(local_max, i)
    if (pow[i] < pow[i-1] && pow[i] < pow[i+1]) local_min <- c(local_min, i)
  }
  
  # ── Count non-monotone decreases ─────────────────────────
  decreases <- which(diff(pow) < -1e-15)
  
  # ── Proposition 3.3 bound check ──────────────────────────
  bound <- compute_prop33_bound(alpha, Bmax)
  actual_drop <- c(pmax(diff(-pow), 0), 0)  # Pow(B) - Pow(B+1)
  # Only check on plateaus
  plateaus <- which(diff(kB) == 0)
  bound_ok <- all(actual_drop[plateaus] <= bound[plateaus] + 1e-12)
  
  # ── Print summary ────────────────────────────────────────
  cat(sprintf("═══ n = %d, p1 = %.3f, alpha = %.4f, Bmax = %d ═══\n",
              n, p1, alpha, Bmax))
  cat(sprintf("  q(0)          = %.6f\n", q[1]))
  for (s in 1:min(n, 8)) {
    cat(sprintf("  q(%d)          = %.6f%s\n", s, q[s+1],
                ifelse(q[s+1] < alpha, "  < alpha  ← threshold", "")))
  }
  cat(sprintf("  Pow_exact     = %.6f\n", pe))
  cat(sprintf("  Non-monotone decreases: %d / %d steps\n",
              length(decreases), Bmax - 1))
  cat(sprintf("  Detected local maxima:  %d\n", length(local_max)))
  cat(sprintf("  Aligned B (first 8):    %s\n",
              paste(head(aligned, 8), collapse = ", ")))
  cat(sprintf("  Local max == aligned?   %s\n",
              ifelse(all(local_max %in% aligned) && all(aligned %in% local_max),
                     "YES", "PARTIALLY")))
  mc_exceeds <- sum(pow > pe + 1e-12)
  cat(sprintf("  MC > exact power:       %d / %d values of B\n",
              mc_exceeds, Bmax))
  cat(sprintf("  Max MC power:           %.6f at B = %d\n",
              max(pow), which.max(pow)))
  cat(sprintf("  Prop 3.3 bound holds:   %s\n", ifelse(bound_ok, "YES", "NO")))
  
  # ── Plot ──────────────────────────────────────────────────
  if (is.null(main_title)) {
    main_title <- paste0("n = ", n, ",  p1 = ", p1, ",  α = ", alpha)
  }
  
  par(mar = c(4.5, 4.5, 3, 1), family = "serif")
  plot(Bseq, pow, type = "l", col = "steelblue", lwd = 1.5,
       xlim = c(1, Bmax), ylim = c(0, max(max(pow), pe) * 1.15),
       xaxs = "i", yaxs = "i",
       xlab = expression(italic(B)), ylab = "Power",
       main = main_title,
       panel.first = grid(col = "grey90"))
  abline(h = pe, lty = 2, col = "orangered", lwd = 1.5)
  points(aligned, pow[aligned], pch = 16, cex = 0.8, col = "forestgreen")
  
  legend("bottomright",
         legend = c(expression(Pow(italic(B))),
                    expression(Pow[exact]),
                    expression(alpha %.% (italic(B)+1) %in% italic(N))),
         col = c("steelblue", "orangered", "forestgreen"),
         lty = c(1, 2, NA), pch = c(NA, NA, 16),
         lwd = c(1.5, 1.5, NA), pt.cex = 0.8,
         bg = "white", cex = 0.8)
  
  invisible(list(pow = pow, pe = pe, q = q, kB = kB,
                 aligned = aligned, local_max = local_max,
                 decreases = decreases))
}


# ── Sawtooth detail zoom ─────────────────────────────────────
plot_sawtooth_zoom <- function(n, p1, alpha = 0.05,
                               Brange = c(15, 120)) {
  
  Bmax <- Brange[2]
  pow  <- compute_power_curve(n, p1, alpha, Bmax)
  pe   <- compute_exact_power(n, p1, alpha)
  kB   <- compute_kB(alpha, Bmax)
  
  Bseq    <- Brange[1]:Brange[2]
  aligned <- Bseq[abs(alpha * (Bseq + 1) - round(alpha * (Bseq + 1))) < 1e-9]
  
  par(mfrow = c(2, 1), mar = c(4, 4.5, 2.5, 1), family = "serif")
  
  # ─── Top: power curve zoom ───
  plot(Bseq, pow[Bseq], type = "o", pch = 20, cex = 0.4,
       col = "steelblue", lwd = 1.2,
       xlab = expression(italic(B)), ylab = "Power",
       main = paste0("Sawtooth zoom: n = ", n))
  abline(h = pe, lty = 2, col = "orangered", lwd = 1.2)
  points(aligned, pow[aligned], pch = 16, cex = 1.2, col = "forestgreen")
  abline(v = aligned, lty = 3, col = "grey70")
  
  # ─── Bottom: k_B staircase ───
  plot(Bseq, kB[Bseq], type = "s", col = "grey30", lwd = 1.5,
       xlab = expression(italic(B)),
       ylab = expression(italic(k)[italic(B)]),
       main = expression("Critical count " * italic(k)[italic(B)] *
                           " = " * group(lfloor, alpha %.% (italic(B)+1), rfloor) - 1))
  abline(v = aligned, lty = 3, col = "grey70")
  points(aligned, kB[aligned], pch = 16, cex = 1.0, col = "forestgreen")
}


# ── Proposition 3.3 bound visualization ──────────────────────
plot_prop33 <- function(n, p1, alpha = 0.05, Bmax = 500) {
  
  pow   <- compute_power_curve(n, p1, alpha, Bmax)
  kB    <- compute_kB(alpha, Bmax)
  bound <- compute_prop33_bound(alpha, Bmax)
  
  # Actual decrease on plateaus
  actual <- -diff(pow)
  actual[actual < 0] <- 0
  
  # Identify plateaus
  plateaus <- which(diff(kB) == 0)
  
  par(mar = c(4.5, 4.5, 3, 1), family = "serif")
  plot(plateaus, actual[plateaus], pch = 20, cex = 0.4,
       col = "steelblue",
       xlab = expression(italic(B)),
       ylab = expression(Pow(italic(B)) - Pow(italic(B)+1)),
       main = paste0("Prop. 3.3 bound: n = ", n),
       ylim = c(0, max(bound[plateaus]) * 1.1),
       yaxs = "i",
       panel.first = grid(col = "grey90"))
  lines(1:Bmax, bound, col = "orangered", lwd = 1.5, lty = 2)
  legend("topright",
         legend = c("Actual decrease", "Prop. 3.3 upper bound"),
         col = c("steelblue", "orangered"),
         pch = c(20, NA), lty = c(NA, 2), lwd = c(NA, 1.5),
         bg = "white", cex = 0.8)
}


# ── Conditional power heatmap ─────────────────────────────────
plot_conditional_power <- function(alpha = 0.05, Bmax = 200,
                                   nq = 200) {
  
  q_grid <- seq(0.001, 0.2, length.out = nq)
  phi    <- compute_conditional_power(q_grid, alpha, Bmax)
  
  Bseq    <- 1:Bmax
  aligned <- Bseq[abs(alpha * (Bseq + 1) - round(alpha * (Bseq + 1))) < 1e-9]
  
  par(mar = c(4.5, 4.5, 3, 5), family = "serif")
  image(Bseq, q_grid, phi,
        col = hcl.colors(64, "YlOrRd", rev = TRUE),
        xlab = expression(italic(B)),
        ylab = expression(italic(q)(X)),
        main = expression("Conditional power " *
                            phi[italic(B)](italic(q)) *
                            " = P(Bin(" * italic(B) * "," *
                            italic(q) * ") " <= italic(k)[italic(B)] * ")"))
  abline(h = alpha, lty = 2, col = "white", lwd = 1.5)
  abline(v = aligned, lty = 3, col = "grey50", lwd = 0.5)
  contour(Bseq, q_grid, phi, levels = c(0.5, 0.9, 0.95, 0.99),
          add = TRUE, col = "grey30", labcex = 0.7)
}


# ── Parameter sweep: vary p1 ─────────────────────────────────
sweep_p1 <- function(n = 15, alpha = 0.05, Bmax = 400,
                     p1_vals = seq(0.05, 0.40, by = 0.05)) {
  
  cols <- hcl.colors(length(p1_vals), "Zissou 1")
  
  par(mar = c(4.5, 4.5, 3, 1), family = "serif")
  plot(NULL, xlim = c(1, Bmax), ylim = c(0, 1),
       xaxs = "i", yaxs = "i",
       xlab = expression(italic(B)), ylab = "Power",
       main = paste0("Power sweep: n = ", n, ", α = ", alpha),
       panel.first = grid(col = "grey90"))
  
  for (i in seq_along(p1_vals)) {
    pow <- compute_power_curve(n, p1_vals[i], alpha, Bmax)
    pe  <- compute_exact_power(n, p1_vals[i], alpha)
    lines(1:Bmax, pow, col = cols[i], lwd = 1.3)
    abline(h = pe, col = cols[i], lty = 3, lwd = 0.6)
  }
  
  legend("topleft",
         legend = sprintf("p1 = %.2f", p1_vals),
         col = cols, lwd = 1.5, cex = 0.7, bg = "white")
}


# ═══════════════════════════════════════════════════════════════
#  RUN: Reproduce Figure 1
# ═══════════════════════════════════════════════════════════════

cat("══════════════════════════════════════════════════\n")
cat("  Reproducing Figure 1 ...\n")
cat("══════════════════════════════════════════════════\n\n")

fig1 <- reproduce_figure1(Bmax = 400, save_pdf = FALSE)

cat("\n══════════════════════════════════════════════════\n")
cat("  Figure 1 reproduced.\n")
cat("══════════════════════════════════════════════════\n")
cat("
Available functions for experimentation:

  reproduce_figure1(Bmax, save_pdf)
    → Reproduce the paper's Figure 1

  run_experiment(n, p1, alpha, Bmax)
    → Power curve for arbitrary parameters
    → Example: run_experiment(n=20, p1=0.2, alpha=0.05)

  plot_sawtooth_zoom(n, p1, alpha, Brange)
    → Zoom into sawtooth + k_B staircase
    → Example: plot_sawtooth_zoom(15, 0.16, Brange=c(15,120))

  plot_prop33(n, p1, alpha, Bmax)
    → Visualize Proposition 3.3 bound
    → Example: plot_prop33(15, 0.16)

  plot_conditional_power(alpha, Bmax, nq)
    → Heatmap of phi_B(q) = P(Bin(B,q) <= k_B)

  sweep_p1(n, alpha, Bmax, p1_vals)
    → Overlay power curves for multiple p1 values
    → Example: sweep_p1(n=15, p1_vals=seq(0.1, 0.4, 0.05))

  compute_power_curve(n, p1, alpha, Bmax)  [Rcpp]
    → Returns numeric vector of Pow(B), B = 1..Bmax

  compute_exact_power(n, p1, alpha)  [Rcpp]
    → Returns scalar Pow_exact

  compute_q(n)  [Rcpp]
    → Returns vector q(0), q(1), ..., q(n)

  compute_kB(alpha, Bmax)  [Rcpp]
    → Returns integer vector k_1, ..., k_Bmax

  compute_conditional_power(q_grid, alpha, Bmax)  [Rcpp]
    → Returns Bmax x length(q_grid) matrix of phi_B(q)
")