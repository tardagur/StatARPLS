
# StatARPLS

StatARPLS estimates the baseline of a spectrum (FTIR, Raman, UV-Vis)
with asymmetrically reweighted penalised least squares (ARPLS),
implemented from first principles with a compiled Rcpp/RcppArmadillo
solver.

It adds two things the existing R tools do not provide:

- **Automatic smoothing parameter**: lambda is chosen by generalised
  cross-validation (GCV).
- **Baseline stability band**: a moving-block residual bootstrap shows
  how much the fitted baseline moves when the noise is resampled.

Results are returned as one S3 object of class `StatARPLSResult`, with
`print()`, `summary()`, `plot()`, `predict()` and `compare()` methods.

## Installation

The package contains C++ code, so Windows users need
[Rtools](https://cran.r-project.org/bin/windows/Rtools/) installed
first.

``` r
# install.packages("pak")
pak::pak("tardagur/StatARPLS")---
output: github_document
---
```

# StatARPLS

StatARPLS estimates the baseline of a spectrum (FTIR, Raman, UV-Vis)
with asymmetrically reweighted penalised least squares (ARPLS),
implemented from first principles with a compiled Rcpp/RcppArmadillo
solver.

It adds two things the existing R tools do not provide:

- **Automatic smoothing parameter**: lambda is chosen by generalised
  cross-validation (GCV).
- **Baseline stability band**: a moving-block residual bootstrap shows
  how much the fitted baseline moves when the noise is resampled.

Results are returned as one S3 object of class `StatARPLSResult`, with
`print()`, `summary()`, `plot()`, `predict()` and `compare()` methods.

## Installation

The package contains C++ code, so Windows users need
[Rtools](https://cran.r-project.org/bin/windows/Rtools/) installed
first.

``` r
# install.packages("pak")
pak::pak("tardagur/StatARPLS")
```

To read the introduction vignette, clone the repository
(`git clone https://github.com/tardagur/StatARPLS.git`) and, from the
cloned folder, install it with the vignette built (needs knitr,
rmarkdown and pandoc):

``` r
devtools::install(build_vignettes = TRUE, dependencies = TRUE)
browseVignettes("StatARPLS")
```

## Example

``` r
library(StatARPLS)

# a simulated spectrum: curved baseline, two peaks and noise
set.seed(501)
wavenumber <- seq(0, 2000, length.out = 500)
signal <- 50 + 0.02 * wavenumber + 30 * sin(wavenumber / 400) +
  80 * exp(-((wavenumber - 700)^2) / (2 * 25^2)) +
  60 * exp(-((wavenumber - 1300)^2) / (2 * 40^2)) +
  rnorm(length(wavenumber), sd = 1.5)

fit  <- gcv_select_lambda(signal, wavenumber, seed = 1)
boot <- bootstrap_baseline_uncertainty(fit, B = 100, seed = 2)

summary(boot)
plot(boot)
```

`wavenumber` must be strictly increasing. See the introduction vignette
(`vignette("StatARPLS-introduction")`, once installed with vignettes as
above) for a full worked example.

## Limitations

- GCV can pick a very small lambda on some reweighted fits, so the
  baseline follows the peaks. The package warns when lambda is the
  smallest value in the grid; inspect the plot.
- The stability band measures sampling stability at a fixed lambda, not
  the baseline’s bias.
- The linear solve is dense, so large spectra (several thousand points)
  are slow.

<!-- -->


    ## Example


    ``` r
    library(StatARPLS)

    fit  <- gcv_select_lambda(signal, wavenumber, seed = 1)
    boot <- bootstrap_baseline_uncertainty(fit, B = 100, seed = 2)

    summary(boot)
    plot(boot)

`wavenumber` must be strictly increasing. See the introduction vignette
(`vignette("StatARPLS-introduction")`) for a full worked example.

## Limitations

- GCV can pick a very small lambda on some reweighted fits, so the
  baseline follows the peaks. The package warns when lambda is the
  smallest value in the grid; inspect the plot.
- The stability band measures sampling stability at a fixed lambda, not
  the baseline’s bias.
- The linear solve is dense, so large spectra (several thousand points)
  are slow.
