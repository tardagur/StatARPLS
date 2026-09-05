#' Low-level constructor for a StatARPLSResult object
#'
#' @description
#' `new_StatARPLSResult()` is the internal, unchecked constructor for the
#' `StatARPLSResult` class. It performs no validation and exists purely to
#' assemble the object's fields cheaply; it is called by the exported
#' `StatARPLSResult()` helper and by internal fitting routines (e.g. during
#' bootstrap resampling) where repeated validation would be wasteful.
#'
#' @param signal A numeric vector of raw, uncorrected spectral intensities.
#' @param wavenumber A numeric vector of the same length as `signal`, giving
#'   the independent axis (wavenumber or index) at which each signal value
#'   was recorded.
#' @param baseline A numeric vector of the same length as `signal`, giving
#'   the fitted ARPLS baseline estimate.
#' @param corrected A numeric vector of the same length as `signal`, giving
#'   the baseline-corrected signal (`signal - baseline`).
#' @param residual_variance A single numeric value giving the estimated
#'   variance of the corrected residuals.
#' @param stability_bounds A numeric matrix with `length(signal)` rows and
#'   exactly 2 columns (`lower`, `upper`), giving bootstrap-derived pointwise
#'   stability bounds on the fitted baseline.
#' @param lambda A single positive numeric value giving the optimal
#'   regularisation parameter selected for this fit.
#' @param n_iter A single positive integer giving the number of iterations
#'   the reweighting scheme required to converge.
#' @param converged A single logical value indicating whether the fitting
#'   algorithm converged within the iteration limit.
#' @param call The matched call of the fitting function that produced this
#'   object, as returned by `match.call()`.
#'
#' @return An object of class `StatARPLSResult`, structured as:
#' * `signal` — numeric vector, length n
#' * `wavenumber` — numeric vector, length n
#' * `baseline` — numeric vector, length n
#' * `corrected` — numeric vector, length n
#' * `residual_variance` — numeric scalar
#' * `stability_bounds` — numeric matrix, n x 2
#' * `lambda` — numeric scalar
#' * `n_iter` — integer scalar
#' * `converged` — logical scalar
#' * `call` — language object
#'
#' @noRd
new_StatARPLSResult <- function(signal,
                                wavenumber,
                                baseline,
                                corrected,
                                residual_variance,
                                stability_bounds,
                                lambda,
                                n_iter,
                                converged,
                                call) {
  structure(
    list(
      signal            = signal,
      wavenumber        = wavenumber,
      baseline          = baseline,
      corrected         = corrected,
      residual_variance = residual_variance,
      stability_bounds  = stability_bounds,
      lambda            = lambda,
      n_iter            = n_iter,
      converged         = converged,
      call              = call
    ),
    class = "StatARPLSResult"
  )
}

#' Validate a StatARPLSResult object's structural invariants
#'
#' @description
#' Enforces the shape, type, and length constraints required for a
#' `StatARPLSResult` object to be internally consistent. Called by
#' `StatARPLSResult()` immediately after construction; never bypass this
#' step when building an object that will be handed to the user.
#'
#' @param x An object of class `StatARPLSResult`, typically produced by
#'   `new_StatARPLSResult()`.
#'
#' @return The input `x`, invisibly, if all checks pass. Throws an
#' informative error via `stop()` on the first violated invariant:
#' * length consistency across `signal`, `wavenumber`, `baseline`,
#'   `corrected`
#' * `stability_bounds` has exactly 2 columns and `length(signal)` rows
#' * `lambda` is a strictly positive finite scalar
#' * `n_iter` is a positive integer scalar
#' * `converged` is a single logical value
#'
#' @noRd
validate_StatARPLSResult <- function(x) {
  n <- length(x$signal)

  if (length(x$wavenumber) != n || length(x$baseline) != n ||
      length(x$corrected) != n) {
    stop(
      "`signal`, `wavenumber`, `baseline`, and `corrected` must all be ",
      "the same length.",
      call. = FALSE
    )
  }

  if (!is.matrix(x$stability_bounds) ||
      nrow(x$stability_bounds) != n ||
      ncol(x$stability_bounds) != 2) {
    stop(
      "`stability_bounds` must be a numeric matrix with ", n,
      " rows and 2 columns (lower, upper).",
      call. = FALSE
    )
  }

  if (!is.numeric(x$lambda) || length(x$lambda) != 1L ||
      !is.finite(x$lambda) || x$lambda <= 0) {
    stop("`lambda` must be a single strictly positive finite number.", call. = FALSE)
  }

  if (!is.numeric(x$n_iter) || length(x$n_iter) != 1L || x$n_iter < 1) {
    stop("`n_iter` must be a single positive integer.", call. = FALSE)
  }

  if (!is.logical(x$converged) || length(x$converged) != 1L) {
    stop("`converged` must be a single logical value.", call. = FALSE)
  }

  invisible(x)
}

#' Construct a StatARPLSResult object
#'
#' @description
#' User-facing constructor for the `StatARPLSResult` class, the single
#' return type produced by every major fitting routine in the StatARPLS
#' package. Combines the raw signal, fitted ARPLS baseline, bootstrap
#' stability diagnostics, and selected regularisation parameter into one
#' validated S3 object.
#'
#' @param signal A numeric vector of raw, uncorrected spectral intensities.
#' @param wavenumber A numeric vector of the same length as `signal`, giving
#'   the independent axis at which each signal value was recorded.
#' @param baseline A numeric vector of the same length as `signal`, giving
#'   the fitted ARPLS baseline estimate.
#' @param corrected A numeric vector of the same length as `signal`, giving
#'   the baseline-corrected signal.
#' @param residual_variance A single numeric value giving the estimated
#'   variance of the corrected residuals.
#' @param stability_bounds A numeric matrix with `length(signal)` rows and
#'   2 columns (`lower`, `upper`), giving bootstrap stability bounds.
#' @param lambda A single positive numeric value giving the optimal
#'   regularisation parameter.
#' @param n_iter A single positive integer giving the number of reweighting
#'   iterations to convergence.
#' @param converged A single logical value indicating convergence status.
#' @param call The matched call of the fitting function; defaults to the
#'   caller's `match.call()` when not supplied.
#'
#' @return A validated object of class `StatARPLSResult` containing:
#' * `signal`, `wavenumber`, `baseline`, `corrected` — numeric vectors,
#'   each of length n
#' * `residual_variance` — numeric scalar
#' * `stability_bounds` — numeric matrix, n x 2
#' * `lambda` — numeric scalar
#' * `n_iter` — integer scalar
#' * `converged` — logical scalar
#' * `call` — language object recording how the object was produced
#'
#' @export
StatARPLSResult <- function(signal,
                            wavenumber,
                            baseline,
                            corrected,
                            residual_variance,
                            stability_bounds,
                            lambda,
                            n_iter,
                            converged,
                            call = sys.call(-1)) {
  result <- new_StatARPLSResult(
    signal            = signal,
    wavenumber        = wavenumber,
    baseline          = baseline,
    corrected         = corrected,
    residual_variance = residual_variance,
    stability_bounds  = stability_bounds,
    lambda            = lambda,
    n_iter            = n_iter,
    converged         = converged,
    call              = call
  )

  validate_StatARPLSResult(result)
}
