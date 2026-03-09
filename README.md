# mc-perm-power

Reproducibility code for the paper on non-monotonicity of power in Monte Carlo permutation tests.

## Overview

This repository contains all code needed to reproduce the figures and numerical results in the paper. The paper shows that the unconditional power of the standard (non-randomized) Monte Carlo permutation test is **not** monotone in the number of permutations $B$, and characterizes the sawtooth behavior driven by the integer-valued threshold $k_B = \lfloor \alpha(B+1) \rfloor - 1$.

## Repository structure

```
mc-perm-power/
├── R/
│   ├── Figure1.R            # Figure 1: closed-form Bernoulli example
│   └── Appendix.R           # Appendix B: four test statistics (mean diff, MMD², HSIC, energy)
├── src/
│   ├── pow_mc.cpp           # Core Rcpp: power curve computation via Eq. (4)
│   └── appendixB_core.cpp   # Rcpp: test statistic implementations for Appendix B
├── figures/                 # Pre-generated figures (PDF)
├── LICENSE
└── README.md
```

## Requirements

- **R** >= 4.0
- **Rcpp** (`install.packages("Rcpp")`)
- A C++ compiler supporting C++11 (e.g., g++ or clang++)

No other R packages are required.

## Reproducing the figures

All scripts are self-contained. Run from the repository root:

```r
# Figure 1 (Section 4: Bernoulli example, Cases 1 & 2)
source("R/Figure1.R")

# Appendix B (four test statistics: mean diff, MMD², HSIC, energy distance)
source("R/Appendix.R")
```

Each script automatically compiles the required C++ source files via `Rcpp::sourceCpp()` and saves the output figures to `figures/`.

## License

[MIT](LICENSE)
