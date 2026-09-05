#' Print a StatARPLSResult object
#'
#' @description
#' Prints a concise, human-readable summary of a fitted `StatARPLSResult`
#' object: signal length, convergence status, iteration count, and the
#' selected regularisation parameter. If the object carries GCV grid
#' search or bootstrap metadata (attached by `gcv_select_lambda()` or
#' `bootstrap_baseline_uncertainty()`), a brief summary of each is
#' appended. Full residual diagnostics remain available via
#' `summary.StatARPLSResult()`.
#'
#' @param x An object of class `StatARPLSResult`.
#' @param ... Additional arguments passed to or from other methods; unused
#'   here and retained only for S3 signature compatibility with
#'   `print.default()`.
#'
#' @return `x`, invisibly. Called for its side effect of printing to the
#' console:
#' * package identity line and fitted call
#' * signal length and convergence status
#' * selected lambda and iteration count
#' * GCV grid summary, if present as an attribute
#' * bootstrap replicate summary, if present as an attribute
#'
#' @export
print.StatARPLSResult <- function(x, ...) {
  cat("<StatARPLSResult>\n")
  cat("Call:", deparse(x$call), "\n\n")
  cat("Signal length:   ", length(x$signal), "\n")
  cat("Converged:       ", x$converged, "in", x$n_iter, "iteration(s)\n")
  cat("Lambda (optimal):", format(x$lambda, digits = 4), "\n")
  cat("Residual variance:", format(x$residual_variance, digits = 4), "\n")

  gcv_info <- attr(x, "gcv")
  if (!is.null(gcv_info)) {
    cat("\nGCV selection: ", length(gcv_info$lambda_grid), "candidates evaluated,",
        "effective df =", format(gcv_info$df_effective, digits = 4), "\n")
  }

  bootstrap_info <- attr(x, "bootstrap")
  if (!is.null(bootstrap_info)) {
    cat("Bootstrap:      ", bootstrap_info$B_effective, "of",
        bootstrap_info$B_requested, "replicates converged",
        paste0("(", bootstrap_info$n_failed, " failed)"), "\n")
  }

  invisible(x)
}

#' Summarise a StatARPLSResult object
#'
#' @description
#' Computes residual and stability-band diagnostics for a fitted
#' `StatARPLSResult` object and returns them as a `summary.StatARPLSResult`
#' object with its own `print()` method. When available, GCV grid search
#' and bootstrap replicate metadata are carried through into the summary.
#'
#' @param object An object of class `StatARPLSResult`.
#' @param ... Additional arguments passed to or from other methods; unused
#'   here and retained only for S3 signature compatibility with
#'   `summary.default()`.
#'
#' @return An object of class `summary.StatARPLSResult`, a list containing:
#' * `lambda` — the selected regularisation parameter
#' * `n_iter`, `converged` — fitting diagnostics
#' * `residual_variance` — as stored on the original object
#' * `residual_quantiles` — named numeric vector, the 0/25/50/75/100th
#'   percentiles of `corrected`
#' * `mean_band_width` — numeric scalar, the mean pointwise width of
#'   `stability_bounds`
#' * `gcv` — the object's `"gcv"` attribute, or `NULL` if absent
#' * `bootstrap` — the object's `"bootstrap"` attribute, or `NULL` if
#'   absent
#'
#' @export
summary.StatARPLSResult <- function(object, ...) {
  residual_quantiles <- stats::quantile(
    object$corrected,
    probs = c(0, 0.25, 0.5, 0.75, 1)
  )

  mean_band_width <- mean(
    object$stability_bounds[, 2] - object$stability_bounds[, 1]
  )

  structure(
    list(
      lambda             = object$lambda,
      n_iter             = object$n_iter,
      converged          = object$converged,
      residual_variance  = object$residual_variance,
      residual_quantiles = residual_quantiles,
      mean_band_width    = mean_band_width,
      gcv                = attr(object, "gcv"),
      bootstrap          = attr(object, "bootstrap")
    ),
    class = "summary.StatARPLSResult"
  )
}

#' Print a summary.StatARPLSResult object
#'
#' @param x An object of class `summary.StatARPLSResult`.
#' @param ... Additional arguments passed to or from other methods; unused
#'   here and retained only for S3 signature compatibility.
#'
#' @return `x`, invisibly. Called for its side effect of printing a
#' formatted diagnostics table to the console, including GCV and
#' bootstrap sections when the underlying data is present.
#'
#' @export
print.summary.StatARPLSResult <- function(x, ...) {
  cat("<summary.StatARPLSResult>\n\n")
  cat("Lambda:             ", format(x$lambda, digits = 4), "\n")
  cat("Converged:          ", x$converged, "in", x$n_iter, "iteration(s)\n")
  cat("Residual variance:  ", format(x$residual_variance, digits = 4), "\n")
  cat("Mean stability-band width:", format(x$mean_band_width, digits = 4), "\n\n")
  cat("Corrected-signal quantiles:\n")
  print(round(x$residual_quantiles, 4))

  if (!is.null(x$gcv)) {
    cat("\n--- GCV Parameter Selection ---\n")
    cat("Candidates evaluated:", length(x$gcv$lambda_grid), "\n")
    cat("Optimal lambda:      ", format(x$gcv$lambda_optimal, digits = 4), "\n")
    cat("Effective df:        ", format(x$gcv$df_effective, digits = 4), "\n")
    cat("Min GCV score:       ", format(min(x$gcv$gcv_scores), digits = 4), "\n")
  }

  if (!is.null(x$bootstrap)) {
    cat("\n--- Residual Bootstrap ---\n")
    cat("Replicates requested:", x$bootstrap$B_requested, "\n")
    cat("Replicates converged:", x$bootstrap$B_effective,
        paste0("(", x$bootstrap$n_failed, " failed)"), "\n")
    cat("Mean pointwise SE:   ",
        format(mean(x$bootstrap$boot_se), digits = 4), "\n")
  }

  invisible(x)
}

#' Plot a StatARPLSResult object
#'
#' @description
#' Produces a base R graphics overlay showing the raw signal, the fitted
#' ARPLS baseline, and the bootstrap stability band, on a single plotting
#' surface. No third-party plotting engine is used; this method depends
#' only on the `graphics` package already declared in `Imports`.
#'
#' @param x An object of class `StatARPLSResult`.
#' @param show_band A single logical value; if `TRUE` (the default), the
#'   bootstrap stability band from `x$stability_bounds` is drawn as a
#'   shaded polygon beneath the baseline.
#' @param ... Additional arguments passed through to the underlying
#'   `graphics::plot()` call (e.g. `main`, `xlab`, `col`).
#'
#' @return `NULL`, invisibly. Called for its side effect of drawing to the
#' active graphics device:
#' * raw signal as a line trace
#' * fitted baseline overlaid in a contrasting colour
#' * optional shaded stability band (`show_band = TRUE`)
#' * legend identifying each series
#'
#' @export
plot.StatARPLSResult <- function(x, show_band = TRUE, ...) {
  wavenumber <- x$wavenumber

  graphics::plot(
    wavenumber, x$signal,
    type = "l",
    xlab = "Wavenumber",
    ylab = "Intensity",
    ...
  )

  if (isTRUE(show_band)) {
    graphics::polygon(
      x   = c(wavenumber, rev(wavenumber)),
      y   = c(x$stability_bounds[, 1], rev(x$stability_bounds[, 2])),
      col = grDevices::adjustcolor("grey", alpha.f = 0.4),
      border = NA
    )
  }

  graphics::lines(wavenumber, x$baseline, col = "red", lwd = 2)

  graphics::legend(
    "topright",
    legend = c("Raw signal", "Fitted baseline", "Stability band"),
    col    = c("black", "red", "grey"),
    lty    = c(1, 1, NA),
    pch    = c(NA, NA, 15),
    bty    = "n"
  )

  invisible(NULL)
}
#' Predict baseline values on a new wavenumber grid
#'
#' @description
#' Interpolates (and, if explicitly requested, extrapolates) a fitted
#' ARPLS baseline onto an arbitrary target wavenumber grid, using
#' `stats::approx()` for linear interpolation or `stats::spline()` for
#' smooth cubic interpolation. By default, target points falling outside
#' the range of the originally fitted `wavenumber` vector return `NA`
#' rather than a guessed value, since the fitted ARPLS baseline carries no
#' statistical information beyond the observed spectrum.
#'
#' @param object An object of class `StatARPLSResult`.
#' @param newdata A numeric vector of new wavenumber coordinates at which
#'   to predict the fitted baseline.
#' @param method A single character string, either `"linear"` (the
#'   default, dispatches to `stats::approx()`) or `"spline"` (dispatches
#'   to `stats::spline()` for smooth cubic interpolation).
#' @param rule A single integer, `1` or `2`, controlling behaviour for
#'   `newdata` points outside the fitted wavenumber range. `1` (the
#'   default) returns `NA` for out-of-range points; `2` extends the
#'   boundary value as a constant. Passed through to `stats::approx()`;
#'   ignored (with a warning) for `method = "spline"`, which always
#'   extrapolates.
#' @param ... Additional arguments passed to or from other methods; unused
#'   here and retained only for S3 signature compatibility with
#'   `predict.default()`.
#'
#' @return A numeric vector the same length as `newdata`, giving the
#' predicted baseline value at each requested wavenumber:
#' * in-range points — interpolated from the fitted `wavenumber` /
#'   `baseline` pair
#' * out-of-range points — `NA` under `rule = 1`, or the nearest fitted
#'   boundary value under `rule = 2`
#'
#' @export
predict.StatARPLSResult <- function(object,
                                    newdata,
                                    method = c("linear", "spline"),
                                    rule = 1L,
                                    ...) {
  method <- match.arg(method)

  if (!is.numeric(newdata)) {
    stop("`newdata` must be a numeric vector of wavenumber coordinates.", call. = FALSE)
  }

  if (!rule %in% c(1L, 2L)) {
    stop("`rule` must be 1 (NA outside range) or 2 (constant extrapolation).", call. = FALSE)
  }

  if (method == "linear") {
    prediction <- stats::approx(
      x    = object$wavenumber,
      y    = object$baseline,
      xout = newdata,
      rule = rule
    )$y
  } else {
    if (rule == 1L) {
      warning(
        "`rule = 1` (NA outside range) is not supported by ",
        "`method = \"spline\"`, which always extrapolates. Out-of-range ",
        "predictions should be treated with caution.",
        call. = FALSE
      )
    }

    fitted_range <- range(object$wavenumber)
    prediction <- stats::spline(
      x    = object$wavenumber,
      y    = object$baseline,
      xout = newdata
    )$y

    if (rule == 1L) {
      prediction[newdata < fitted_range[1] | newdata > fitted_range[2]] <- NA_real_
    }
  }

  prediction
}
#' Compare multiple StatARPLSResult fits
#'
#' @description
#' Generic function for side-by-side comparison of two or more fitted
#' model objects. `StatARPLSResult` provides the sole method in this
#' package; see `compare.StatARPLSResult()`.
#'
#' @param x The primary object to compare.
#' @param ... Additional objects to compare against `x`, along with any
#'   further arguments passed to the dispatched method.
#'
#' @return A comparison object; the exact class and structure depend on
#' the method invoked. See `compare.StatARPLSResult()` for the
#' `StatARPLSResult` case.
#'
#' @export
compare <- function(x, ...) {
  UseMethod("compare")
}

#' Compare multiple ARPLS baseline fits side by side
#'
#' @description
#' Assembles a comparative diagnostics table across two or more
#' `StatARPLSResult` objects, contrasting their selected regularisation
#' parameter, effective degrees of freedom (where GCV metadata is
#' present), iterations to convergence, and mean stability-band width.
#' Useful for contrasting fits at different lambda values, on different
#' spectra, or before/after applying `bootstrap_baseline_uncertainty()`.
#'
#' @param x An object of class `StatARPLSResult`, used as the first row
#'   of the comparison table.
#' @param ... One or more additional objects of class `StatARPLSResult`
#'   to compare against `x`. Objects of any other class trigger an error.
#' @param labels An optional character vector of row labels, the same
#'   length as the total number of objects being compared (`x` plus all
#'   `...` arguments). Defaults to the deparsed argument expressions
#'   (e.g. `"fit_low_lambda"`) when not supplied.
#'
#' @return An object of class `compare.StatARPLSResult`, a data frame
#' with one row per compared fit and columns:
#' * `label` — character, the fit's identifying label
#' * `lambda` — numeric, the fitted regularisation parameter
#' * `df_effective` — numeric, effective degrees of freedom from GCV
#'   metadata, or `NA` if the fit was not produced via
#'   `gcv_select_lambda()`
#' * `n_iter` — integer, reweighting iterations to convergence
#' * `converged` — logical, convergence status
#' * `mean_band_width` — numeric, mean pointwise width of
#'   `stability_bounds`
#'
#' @export
compare.StatARPLSResult <- function(x, ..., labels = NULL) {
  dots <- list(...)
  dot_names <- names(dots)

  # Separate trailing `labels =` from genuine StatARPLSResult objects,
  # in case it arrived positionally through `...` rather than by name.
  is_fit <- vapply(dots, inherits, logical(1), what = "StatARPLSResult")
  extra_fits <- dots[is_fit]

  all_fits <- c(list(x), extra_fits)

  if (length(all_fits) < 2) {
    stop(
      "`compare()` requires at least two StatARPLSResult objects: `x` ",
      "plus one or more via `...`.",
      call. = FALSE
    )
  }

  if (!all(vapply(all_fits, inherits, logical(1), what = "StatARPLSResult"))) {
    stop("All objects passed to `compare()` must be of class StatARPLSResult.", call. = FALSE)
  }

  if (is.null(labels)) {
    call_args <- match.call(expand.dots = TRUE)
    arg_names <- as.character(call_args)[-1]
    labels <- arg_names[seq_along(all_fits)]
  } else if (length(labels) != length(all_fits)) {
    stop(
      "`labels` must have length ", length(all_fits),
      " (one per compared fit); got ", length(labels), ".",
      call. = FALSE
    )
  }

  extract_row <- function(fit) {
    gcv_info <- attr(fit, "gcv")
    band_width <- mean(fit$stability_bounds[, 2] - fit$stability_bounds[, 1])

    data.frame(
      lambda           = fit$lambda,
      df_effective      = if (is.null(gcv_info)) NA_real_ else gcv_info$df_effective,
      n_iter           = fit$n_iter,
      converged        = fit$converged,
      mean_band_width  = band_width,
      stringsAsFactors = FALSE
    )
  }

  comparison_rows <- do.call(rbind, lapply(all_fits, extract_row))
  comparison_table <- cbind(label = labels, comparison_rows, stringsAsFactors = FALSE)

  structure(comparison_table, class = c("compare.StatARPLSResult", "data.frame"))
}

#' Print a compare.StatARPLSResult object
#'
#' @param x An object of class `compare.StatARPLSResult`.
#' @param ... Additional arguments passed to or from other methods; unused
#'   here and retained only for S3 signature compatibility.
#'
#' @return `x`, invisibly. Called for its side effect of printing the
#' comparison table to the console as a formatted data frame.
#'
#' @export
print.compare.StatARPLSResult <- function(x, ...) {
  cat("<compare.StatARPLSResult>\n\n")
  print.data.frame(x, row.names = FALSE)
  invisible(x)
}
