#' @import Matrix
NULL

#' Construct a second-order sparse difference operator
#' @noRd
arpls_difference_matrix <- function(n) {
  if (n < 3) {
    stop("`n` must be at least 3 to fit a second-order ARPLS baseline.", call. = FALSE)
  }
  row_index <- rep(seq_len(n - 2), times = 3)
  col_index <- c(seq_len(n - 2), seq_len(n - 2) + 1, seq_len(n - 2) + 2)
  band_values <- rep(c(1, -2, 1), each = n - 2)

  Matrix::sparseMatrix(
    i = row_index, j = col_index, x = band_values,
    dims = c(n - 2, n)
  )
}

#' Update weight vector based on asymmetric residuals
#' @noRd
arpls_update_weights <- function(residual, ratio = 0.05) {
  negative_residual <- residual[residual < 0]

  if (length(negative_residual) < 2) {
    return(rep(1, length(residual)))
  }

  location <- mean(negative_residual)
  scale <- stats::sd(negative_residual)

  if (!is.finite(scale) || scale <= 1e-10) {
    return(rep(1, length(residual)))
  }

  exponent <- 2 * (residual - (2 * scale - location)) / scale
  exponent <- pmin(exponent, 500)
  w <- 1 / (1 + exp(exponent))

  if (any(!is.finite(w))) {
    stop(
      "arpls_update_weights() produced non-finite weights from a ",
      "degenerate residual distribution; check the input spectrum for ",
      "constant or near-constant regions.",
      call. = FALSE
    )
  }

  w
}

#' Validate the wavenumber axis supplied to a fitting routine
#'
#' @description
#' ARPLS penalises the second difference of the baseline, which implicitly
#' assumes the points are evenly spaced along a strictly increasing axis.
#' This helper rejects axes that violate the first assumption outright
#' (repeated values, e.g. from overlapping step-scan segments, or a
#' descending axis, as is common for FTIR files) and warns about axes that
#' violate the second.
#'
#' @param wavenumber A numeric vector of at least 3 points.
#' @param spacing_tolerance A single positive number. A warning is raised
#'   when `(max(step) - min(step)) / median(step)` exceeds this value.
#'
#' @return `wavenumber`, invisibly, if the axis is strictly increasing.
#' Otherwise an informative error (non-finite, repeated, descending, or
#' non-monotone values) or warning (uneven spacing).
#'
#' @noRd
validate_wavenumber_axis <- function(wavenumber, spacing_tolerance = 0.05) {
  bad_wavenumber <- which(!is.finite(wavenumber))
  if (length(bad_wavenumber) > 0) {
    stop(
      "`wavenumber` contains non-finite values at position(s): ",
      paste(utils::head(bad_wavenumber, 5), collapse = ", "),
      if (length(bad_wavenumber) > 5) ", ...",
      call. = FALSE
    )
  }

  step <- diff(wavenumber)
  non_increasing <- which(step <= 0)
  if (length(non_increasing) > 0) {
    stop(
      "`wavenumber` must be strictly increasing, but ", length(non_increasing),
      " step(s) are zero or negative (first after position(s): ",
      paste(utils::head(non_increasing, 5), collapse = ", "),
      if (length(non_increasing) > 5) ", ...",
      "). Sort `signal` and `wavenumber` together first (e.g. with ",
      "`order(wavenumber)`), and average repeated values, such as ",
      "overlapping scan segments, by interpolating onto a regular grid ",
      "with `stats::approx()`.",
      call. = FALSE
    )
  }

  spread <- (max(step) - min(step)) / stats::median(step)
  if (spread > spacing_tolerance) {
    warning(
      "`wavenumber` is unevenly spaced (step sizes range from ",
      signif(min(step), 3), " to ", signif(max(step), 3), "). ARPLS assumes ",
      "evenly spaced points, so the fitted baseline may be distorted; ",
      "consider interpolating onto a regular grid first (e.g. with ",
      "`stats::approx()`).",
      call. = FALSE
    )
  }

  invisible(wavenumber)
}

#' Extract the noise-level residuals from a baseline-corrected signal
#'
#' @description
#' The corrected signal (`signal - baseline`) contains every spectral peak,
#' so its variance reflects peak heights, not measurement noise. This helper
#' keeps only the points that the asymmetric weight function classifies as
#' baseline (weight above 0.5, using the same logistic rule as the fitting
#' loop), so that what remains is a noise estimate. It is the single
#' definition of "noise residuals" shared by the reported residual
#' variance, the provisional stability band, and the bootstrap residual
#' pool.
#'
#' @param corrected A numeric vector, `signal - baseline`.
#'
#' @return A numeric vector of the corrected values at baseline points. If
#' fewer than 3 points are classified as baseline, a warning is raised and
#' the full `corrected` vector is returned instead.
#'
#' @noRd
arpls_noise_residuals <- function(corrected) {
  weights <- arpls_update_weights(corrected)
  is_baseline_point <- weights > 0.5

  if (sum(is_baseline_point) < 3L) {
    warning(
      "Fewer than 3 points were classified as baseline when estimating the ",
      "noise level; using all points, which overstates the noise if the ",
      "spectrum contains peaks.",
      call. = FALSE
    )
    return(corrected)
  }

  corrected[is_baseline_point]
}

#' Noise variance and provisional stability band for a fitted baseline
#'
#' @param baseline A numeric vector, the fitted baseline.
#' @param corrected A numeric vector, `signal - baseline`.
#'
#' @return A list with `residual_variance` (numeric scalar, the variance of
#' the noise residuals) and `stability_bounds` (numeric matrix, n x 2,
#' `baseline` plus or minus 1.96 noise standard deviations).
#'
#' @noRd
arpls_residual_summary <- function(baseline, corrected) {
  residual_variance <- stats::var(arpls_noise_residuals(corrected))
  noise_sd <- sqrt(residual_variance)

  list(
    residual_variance = residual_variance,
    stability_bounds = cbind(
      baseline - 1.96 * noise_sd,
      baseline + 1.96 * noise_sd
    )
  )
}

#' Pure-R reference implementation of the ARPLS IRLS loop
#'
#' @description
#' A pure-R, dense-matrix implementation of the same asymmetrically
#' reweighted penalized least squares fitting loop as `arpls_fit_cpp()`
#' (compiled; see `src/arpls_solver.cpp`). This function is **not** the
#' engine used by `fit_arpls()`, `gcv_select_lambda()`, or
#' `bootstrap_baseline_uncertainty()` -- `arpls_engine()` dispatches to the
#' compiled implementation. This reference implementation is retained for:
#' * the introduction vignette's R-vs-compiled benchmark
#' * a pure-R fallback path, should one ever be wired in for environments
#'   without a working compiler toolchain
#'
#' Both engines currently assemble and solve the system densely: this one
#' densifies the sparse difference-operator penalty matrix, and the
#' compiled engine builds a dense `arma::mat` (see the note at the top of
#' `src/arpls_solver.cpp` on why a sparse solve was not shipped). Its
#' convergence criterion (L1 relative weight change, normalised by the old
#' weights) must be kept identical to `arpls_fit_cpp()`'s if either is ever
#' modified; the two implementations drifted apart once already and are
#' now covered by a parity test (`test-fit-arpls.R`) specifically to catch
#' that again.
#'
#' @param signal A numeric vector of raw spectral intensities.
#' @param lambda A single positive numeric regularisation parameter.
#' @param ratio A single positive numeric convergence threshold.
#' @param max_iter A single positive integer, maximum reweighting iterations.
#' @param ridge A single small positive numeric diagonal stabiliser.
#'
#' @return A list with elements `baseline`, `weights`, `n_iter`, `converged`
#' -- identical in shape to `arpls_fit_cpp()`'s return value.
#'
#' @noRd
arpls_engine_r <- function(signal, lambda, ratio, max_iter, ridge) {
  n <- length(signal)
  w <- rep(1, n)
  D <- arpls_difference_matrix(n)

  # Dense by construction -- this is the reference implementation, not the
  # performance path. See `arpls_fit_cpp()` for the sparse compiled engine.
  P <- lambda * as.matrix(Matrix::crossprod(D))
  diag(P) <- diag(P) + ridge

  converged <- FALSE
  n_iter <- 0L
  baseline <- signal

  for (i in seq_len(max_iter)) {
    n_iter <- i
    w_old <- w

    # Solve the system: (W + P) * z = W * y
    lhs <- diag(w) + P
    rhs <- w * signal
    baseline <- as.numeric(solve(lhs, rhs))

    residual <- signal - baseline
    w <- arpls_update_weights(residual)

    # Check convergence criteria
    weight_diff <- sum(abs(w - w_old)) / sum(abs(w_old))
    if (weight_diff < ratio) {
      converged <- TRUE
      break
    }
  }

  list(
    baseline = baseline,
    weights = w,
    n_iter = n_iter,
    converged = converged
  )
}

#' ARPLS fitting engine dispatcher
#'
#' @description
#' Thin dispatcher used by every fitting, selection, and bootstrap routine
#' in this package. Currently forwards directly to the compiled
#' `arpls_fit_cpp()` (see `src/arpls_solver.cpp`), which solves the dense
#' ARPLS system with Armadillo. The pure-R reference implementation lives
#' separately as `arpls_engine_r()`.
#'
#' @param signal A numeric vector of raw spectral intensities.
#' @param lambda A single positive numeric regularisation parameter.
#' @param ratio A single positive numeric convergence threshold.
#' @param max_iter A single positive integer, maximum reweighting iterations.
#' @param ridge A single small positive numeric diagonal stabiliser.
#'
#' @return A list with elements `baseline`, `weights`, `n_iter`, `converged`.
#'
#' @noRd
arpls_engine <- function(signal, lambda, ratio, max_iter, ridge) {
  arpls_fit_cpp(signal, lambda, ratio, max_iter, ridge)
}

#' Fit an asymmetrically reweighted penalized least squares baseline
#'
#' @description
#' Fits a smooth ARPLS baseline to a spectral signal via iterative
#' reweighted least squares, using the compiled `arpls_fit_cpp()` engine
#' (dispatched via `arpls_engine()`) for the linear solve at each
#' iteration. Implemented entirely from first principles. The reported
#' residual variance and the provisional stability band are based on the
#' noise level estimated from baseline points only, not on the peaks.
#'
#' @param signal A numeric vector of raw spectral intensities.
#' @param wavenumber A numeric vector the same length as `signal`, giving
#'   the independent axis. It must be strictly increasing and should be
#'   evenly spaced; repeated or descending values are an error, and uneven
#'   spacing raises a warning. Defaults to `seq_along(signal)` when not
#'   supplied.
#' @param lambda A single positive numeric value controlling the strength
#'   of the curvature penalty; larger values produce smoother baselines.
#' @param ratio A single positive numeric convergence threshold on the
#'   relative change in the weight vector between iterations.
#' @param max_iter A single positive integer, the maximum number of
#'   reweighting iterations permitted before the routine stops regardless
#'   of convergence.
#' @param ridge A single small positive numeric value added to the
#'   diagonal of the linear system for numerical stability; does not
#'   materially alter the fitted baseline at default magnitude.
#'
#' @return A validated object of class `StatARPLSResult`, as constructed by
#' `StatARPLSResult()`
#'
#' @export
fit_arpls <- function(signal,
                      wavenumber = seq_along(signal),
                      lambda = 1e5,
                      ratio = 0.05,
                      max_iter = 50L,
                      ridge = 1e-6) {
  n <- length(signal)

  if (!is.numeric(signal)) {
    stop("`signal` must be a numeric vector.", call. = FALSE)
  }
  if (n < 3) {
    stop("`signal` must have at least 3 points.", call. = FALSE)
  }
  if (!is.numeric(wavenumber)) {
    stop("`wavenumber` must be a numeric vector.", call. = FALSE)
  }
  if (length(wavenumber) != n) {
    stop("`wavenumber` must be the same length as `signal`.", call. = FALSE)
  }

  bad_signal <- which(!is.finite(signal))
  if (length(bad_signal) > 0) {
    stop(
      "`signal` contains non-finite values (NA/NaN/Inf) at position(s): ",
      paste(utils::head(bad_signal, 5), collapse = ", "),
      if (length(bad_signal) > 5) ", ...",
      call. = FALSE
    )
  }
  bad_wavenumber <- which(!is.finite(wavenumber))
  if (length(bad_wavenumber) > 0) {
    stop(
      "`wavenumber` contains non-finite values at position(s): ",
      paste(utils::head(bad_wavenumber, 5), collapse = ", "),
      if (length(bad_wavenumber) > 5) ", ...",
      call. = FALSE
    )
  }

  validate_wavenumber_axis(wavenumber)

  if (!is.numeric(lambda) || lambda <= 0) {
    stop("`lambda` must be a positive numeric value.", call. = FALSE)
  }
  if (!is.numeric(ratio) || ratio <= 0) {
    stop("`ratio` must be a positive numeric value.", call. = FALSE)
  }
  if (max_iter <= 0) {
    stop("`max_iter` must be an integer greater than 0.", call. = FALSE)
  }

  engine_fit <- arpls_engine(signal, lambda, ratio, max_iter, ridge)

  baseline <- as.numeric(engine_fit$baseline)
  corrected <- signal - baseline
  residual_summary <- arpls_residual_summary(baseline, corrected)

  StatARPLSResult(
    signal            = signal,
    wavenumber        = wavenumber,
    baseline          = baseline,
    corrected         = corrected,
    residual_variance = residual_summary$residual_variance,
    stability_bounds  = residual_summary$stability_bounds,
    lambda            = lambda,
    n_iter            = as.integer(engine_fit$n_iter),
    converged         = as.logical(engine_fit$converged),
    call              = match.call()
  )
}
