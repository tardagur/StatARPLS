#' Generate a synthetic spectrum with known baseline and peaks
#'
#' @description
#' Internal test fixture. Constructs a synthetic spectrum from an
#' independently-specified true baseline plus Gaussian peaks and additive
#' noise, so that fitted results can be checked against a known ground
#' truth rather than against the algorithm's own output.
#'
#' @param n A single positive integer, the number of points.
#' @param noise_sd A single non-negative numeric, the standard deviation of
#'   additive Gaussian noise.
#' @param seed A single integer, for reproducible generation.
#'
#' @return A list with elements `wavenumber`, `true_baseline`, `signal`.
#'
#' @noRd
make_synthetic_spectrum <- function(n = 300, noise_sd = 1.0, seed = 1) {
  set.seed(seed)

  wavenumber <- seq(0, 1000, length.out = n)
  true_baseline <- 20 + 0.01 * wavenumber + 10 * sin(wavenumber / 250)
  peak_signal <- 50 * exp(-((wavenumber - 500)^2) / (2 * 20^2))
  noise <- stats::rnorm(n, mean = 0, sd = noise_sd)

  list(
    wavenumber    = wavenumber,
    true_baseline = true_baseline,
    signal        = true_baseline + peak_signal + noise
  )
}

#' Generate a flat, peak-free synthetic spectrum
#'
#' @description
#' Internal test fixture. A pure-noise-on-constant-baseline spectrum, used
#' to test degenerate/edge-case behaviour where there is no peak signal to
#' distinguish from baseline at all.
#'
#' @param n A single positive integer, the number of points.
#' @param level A single numeric, the constant baseline level.
#' @param noise_sd A single non-negative numeric, additive noise SD.
#' @param seed A single integer, for reproducible generation.
#'
#' @return A list with elements `wavenumber`, `true_baseline`, `signal`.
#'
#' @noRd
make_flat_spectrum <- function(n = 100, level = 10, noise_sd = 0.5, seed = 2) {
  set.seed(seed)
  wavenumber <- seq_len(n)
  list(
    wavenumber    = wavenumber,
    true_baseline = rep(level, n),
    signal        = level + stats::rnorm(n, sd = noise_sd)
  )
}
