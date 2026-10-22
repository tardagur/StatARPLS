#' Estimate the trace of the ARPLS hat matrix via Hutchinson's estimator
#'
#' @description
#' Estimates the trace of the effective smoother matrix using Hutchinson's
#' randomized trace estimator. This avoids the cubic computational cost of
#' explicit matrix inversion by running a small number of random probe vectors
#' through the existing sparse solver pipeline.
#'
#' @param weights A numeric vector of converged ARPLS fitting weights.
#' @param penalty_matrix A sparse numeric matrix representing the regularisation penalty.
#' @param ridge A small positive numeric value added to the diagonal for stability.
#' @param n_probe A single positive integer, the number of random probe vectors.
#' @param seed An optional integer used to seed the random number generator.
#'
#' @return A numeric scalar, the estimated trace (effective degrees of freedom).
#' @noRd
arpls_hat_trace <- function(weights, penalty_matrix, ridge, n_probe = 30L, seed = NULL) {
  n <- length(weights)

  if (!is.null(seed)) {
    set.seed(seed)
  }

  probe_matrix <- matrix(
    sample(c(-1, 1), size = n * n_probe, replace = TRUE),
    nrow = n, ncol = n_probe
  )

  weight_matrix <- Matrix::Diagonal(x = weights)
  system_lhs <- weight_matrix + penalty_matrix + Matrix::Diagonal(n, ridge)

  system_rhs <- as.matrix(sweep(probe_matrix, 1, weights, `*`))

  probe_solutions <- as.matrix(Matrix::solve(system_lhs, system_rhs))

  mean(colSums(probe_matrix * probe_solutions))
}

#' Warn when the selected lambda is the smallest value in the search grid
#'
#' @description
#' GCV rewards fitting the data closely, so for reweighted fits it can select
#' the smallest candidate even though the resulting baseline follows the
#' peaks. A selection at the lower edge of the grid is therefore a warning
#' sign. A selection at the upper edge is not flagged: a very smooth baseline
#' is often correct, and larger values change the fit very little. Only
#' meaningful for grids with at least 3 distinct candidates.
#'
#' @param lambda_grid A numeric vector of candidate lambda values.
#' @param best_lambda A single numeric value, the selected lambda.
#'
#' @return `NULL`, invisibly. Called for its side effect of signalling a
#' warning.
#' @noRd
warn_if_lambda_at_grid_edge <- function(lambda_grid, best_lambda) {
  distinct_grid <- sort(unique(lambda_grid))
  if (length(distinct_grid) < 3L || best_lambda != distinct_grid[1L]) {
    return(invisible(NULL))
  }

  warning(
    "The selected `lambda` (", format(best_lambda, scientific = TRUE),
    ") is the smallest value in `lambda_grid`. GCV can favour ",
    "under-smoothing for reweighted fits, so the baseline may be ",
    "following the peaks. Inspect `plot()`, and consider a larger ",
    "minimum `lambda`.",
    call. = FALSE
  )

  invisible(NULL)
}

#' Select the optimal ARPLS regularisation parameter via GCV
#'
#' @description
#' Evaluates a grid of candidate lambda regularisation values and selects the
#' one that minimises the Generalized Cross-Validation (GCV) score, using
#' Hutchinson's randomized trace estimator to determine the effective degrees
#' of freedom efficiently. A grid point whose linear system is singular or
#' near-singular at the candidate lambda is skipped (scored `Inf`) rather
#' than aborting the entire search.
#'
#' @param signal A numeric vector of raw spectral intensities.
#' @param wavenumber A numeric vector the same length as `signal`.
#' @param lambda_grid A numeric vector of candidate regularisation parameters.
#'   Every value must be strictly positive and finite. The default is
#'   `10^seq(2, 8, length.out = 13)`. GCV rewards fitting the data closely
#'   and, for reweighted fits, can select a very small `lambda` whose
#'   baseline follows the peaks. A warning is issued when the selected value
#'   is the smallest candidate; consider a larger minimum `lambda` then.
#' @param ratio A positive numeric convergence threshold.
#' @param max_iter A positive integer specifying the maximum allowed iterations.
#' @param ridge A small numeric value added to the diagonal for stability.
#' @param n_probe An integer, the number of random probes for trace estimation.
#' @param seed An optional integer used to seed the trace estimator random sweeps.
#'   The caller's global RNG state is saved before seeding and restored on
#'   exit, so supplying `seed` never affects random draws elsewhere in the
#'   caller's session.
#'
#' @return A validated object of class `StatARPLSResult` fitted at the optimal
#' lambda, with GCV diagnostics attached via `attr(result, "gcv")`:
#' * `lambda_grid` — numeric vector, the candidate values searched
#' * `gcv_scores` — numeric vector, the GCV score at each candidate
#'   (`Inf` where the fit failed or degrees of freedom collapsed)
#' * `lambda_optimal` — numeric scalar, the selected lambda
#' * `df_effective` — numeric scalar, the effective degrees of freedom at
#'   the selected lambda
#'
#' @examples
#' set.seed(1)
#' wavenumber <- seq(400, 1800, length.out = 300)
#' peaks <- exp(-((wavenumber - 1000) / 25)^2) +
#'   0.6 * exp(-((wavenumber - 1400) / 30)^2)
#' spectrum <- peaks + 0.001 * (wavenumber - 400) + 0.2 +
#'   rnorm(300, sd = 0.01)
#'
#' gcv_fit <- gcv_select_lambda(
#'   spectrum,
#'   wavenumber = wavenumber,
#'   lambda_grid = 10^seq(4, 8, length.out = 5),
#'   seed = 1
#' )
#' attr(gcv_fit, "gcv")$lambda_optimal
#'
#' @export
gcv_select_lambda <- function(signal,
                              wavenumber = seq_along(signal),
                              lambda_grid = 10^seq(2, 8, length.out = 13),
                              ratio = 0.05,
                              max_iter = 50L,
                              ridge = 1e-6,
                              n_probe = 30L,
                              seed = NULL) {
  n <- length(signal)

  if (!is.numeric(signal)) {
    stop("`signal` must be a numeric vector.", call. = FALSE)
  }
  if (n < 3) {
    stop("`signal` must have at least 3 points.", call. = FALSE)
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
  if (!is.numeric(wavenumber) || length(wavenumber) != n) {
    stop("`wavenumber` must be a numeric vector the same length as `signal`.", call. = FALSE)
  }
  validate_wavenumber_axis(wavenumber)
  if (!is.numeric(lambda_grid) || length(lambda_grid) < 1 ||
      any(!is.finite(lambda_grid)) || any(lambda_grid <= 0)) {
    stop(
      "`lambda_grid` must contain at least 1 strictly positive, finite ",
      "numeric value.",
      call. = FALSE
    )
  }
  if (!is.numeric(ratio) || length(ratio) != 1L || ratio <= 0) {
    stop("`ratio` must be a single positive numeric value.", call. = FALSE)
  }
  if (!is.numeric(max_iter) || length(max_iter) != 1L || max_iter <= 0) {
    stop("`max_iter` must be a single integer greater than 0.", call. = FALSE)
  }
  if (!is.numeric(ridge) || length(ridge) != 1L || ridge < 0) {
    stop("`ridge` must be a single non-negative numeric value.", call. = FALSE)
  }
  if (!is.numeric(n_probe) || length(n_probe) != 1L || n_probe < 1) {
    stop("`n_probe` must be a single positive integer.", call. = FALSE)
  }
  if (!is.null(seed) && (!is.numeric(seed) || length(seed) != 1L)) {
    stop("`seed` must be a single integer, or NULL.", call. = FALSE)
  }

  if (!is.null(seed)) {
    has_old_seed <- exists(".Random.seed", envir = .GlobalEnv)
    old_seed <- if (has_old_seed) get(".Random.seed", envir = .GlobalEnv) else NULL
    on.exit({
      if (has_old_seed) {
        assign(".Random.seed", old_seed, envir = .GlobalEnv)
      } else if (exists(".Random.seed", envir = .GlobalEnv)) {
        rm(".Random.seed", envir = .GlobalEnv)
      }
    }, add = TRUE)
  }

  difference_operator <- arpls_difference_matrix(n)
  dtd <- Matrix::crossprod(difference_operator)

  gcv_scores <- numeric(length(lambda_grid))
  engine_fits <- vector("list", length(lambda_grid))
  df_estimates <- rep(NA_real_, length(lambda_grid))

  for (i in seq_along(lambda_grid)) {
    lambda <- lambda_grid[i]

    grid_point <- tryCatch({
      engine_fit <- arpls_engine(signal, lambda, ratio, max_iter, ridge)
      penalty_matrix <- lambda * dtd
      df_effective <- arpls_hat_trace(
        as.numeric(engine_fit$weights), penalty_matrix, ridge,
        n_probe = n_probe, seed = seed
      )
      df_effective <- min(max(df_effective, 0), n - 1)

      residual <- signal - as.numeric(engine_fit$baseline)
      weighted_rss <- sum(engine_fit$weights * residual^2)
      denominator <- n - df_effective

      list(
        engine_fit = engine_fit,
        df_effective = df_effective,
        gcv = if (denominator <= 0.5) Inf else (n * weighted_rss) / denominator^2
      )
    }, error = function(e) NULL)

    if (is.null(grid_point)) {
      gcv_scores[i] <- Inf
      engine_fits[[i]] <- NULL
      df_estimates[i] <- NA_real_
    } else {
      gcv_scores[i] <- grid_point$gcv
      engine_fits[[i]] <- grid_point$engine_fit
      df_estimates[i] <- grid_point$df_effective
    }
  }

  finite_gcv <- is.finite(gcv_scores)
  if (!any(finite_gcv)) {
    stop(
      "No candidate in `lambda_grid` produced a usable GCV score; every ",
      "grid point either failed to fit or drove effective degrees of ",
      "freedom to within 0.5 of n. Try a different `lambda_grid`.",
      call. = FALSE
    )
  }

  best_index <- which(finite_gcv)[which.min(gcv_scores[finite_gcv])]
  best_lambda <- lambda_grid[best_index]

  warn_if_lambda_at_grid_edge(lambda_grid, best_lambda)
  best_fit <- engine_fits[[best_index]]

  baseline_vector <- as.numeric(best_fit$baseline)

  corrected <- signal - baseline_vector
  residual_summary <- arpls_residual_summary(baseline_vector, corrected)

  result <- StatARPLSResult(
    signal            = signal,
    wavenumber        = wavenumber,
    baseline          = baseline_vector,
    corrected         = corrected,
    residual_variance = residual_summary$residual_variance,
    stability_bounds  = residual_summary$stability_bounds,
    lambda    = best_lambda,
    n_iter    = best_fit$n_iter,
    converged = best_fit$converged,
    call      = match.call()
  )

  attr(result, "gcv") <- list(
    lambda_grid    = lambda_grid,
    gcv_scores     = gcv_scores,
    lambda_optimal = best_lambda,
    df_effective   = df_estimates[best_index]
  )

  result
}

#' Draw one moving block bootstrap resample of a residual series
#'
#' @description
#' Constructs a single pseudo-residual series of length `n_out` by
#' concatenating overlapping blocks of length `block_length`, each drawn with
#' replacement from `residual_pool`, then truncating to length `n_out`. The
#' pool may be shorter than `n_out` (it contains only baseline-classified
#' points), so blocks are simply drawn repeatedly. Used by
#' `bootstrap_baseline_uncertainty()` in place of i.i.d. residual resampling:
#' spectroscopic residuals are serially autocorrelated, and i.i.d. resampling
#' destroys that dependence structure, which understates baseline
#' uncertainty relative to a dependence-preserving scheme.
#'
#' @param residual_pool A numeric vector, the residual series to resample
#'   from.
#' @param block_length A single positive integer, the block length. Must be
#'   at least 1 and no greater than `length(residual_pool)`.
#' @param n_out A single positive integer, the length of the series to
#'   return. Defaults to `length(residual_pool)`.
#'
#' @return A numeric vector of length `n_out`.
#' @noRd
moving_block_bootstrap_sample <- function(residual_pool,
                                          block_length,
                                          n_out = length(residual_pool)) {
  n <- length(residual_pool)
  n_blocks_needed <- ceiling(n_out / block_length)
  max_start <- n - block_length + 1L

  block_starts <- sample(seq_len(max_start), size = n_blocks_needed, replace = TRUE)

  resampled <- unlist(lapply(block_starts, function(start) {
    residual_pool[start:(start + block_length - 1L)]
  }))

  resampled[seq_len(n_out)]
}

#' Quantify baseline uncertainty via moving block residual bootstrap
#'
#' @description
#' Resamples the residuals of a converged fit using a moving block bootstrap
#' (overlapping blocks of contiguous residuals, drawn with replacement) to
#' construct bootstrap pseudo-signals, refits the ARPLS baseline at the
#' original fit's lambda for each replicate, and derives pointwise standard
#' errors and percentile bounds, centred on the fitted baseline. Only
#' residuals at points classified as baseline (final ARPLS weight above 0.5)
#' are resampled, so peaks do not leak into the pseudo-signals; because those
#' points are not contiguous once peaks are removed, block dependence is
#' preserved only approximately across the gaps. The block scheme is used instead
#' of i.i.d. residual resampling because spectroscopic residuals are
#' serially autocorrelated; i.i.d. resampling would understate uncertainty.
#' Lambda is held fixed across replicates so the bootstrap quantifies
#' baseline instability at the selected smoothing level, not uncertainty in
#' the smoothing-parameter choice itself.
#'
#' @param fit A validated object of class `StatARPLSResult`.
#' @param B An integer specifying the number of bootstrap replicates (default 200).
#' @param block_length A single positive integer, the moving block length.
#'   Defaults to `NULL`, which uses the heuristic
#'   `max(2, round(n^(1/3)))`, clamped to at most `floor(n / 2)` so that at
#'   least two distinct block start positions exist. Supply this explicitly
#'   when you have domain knowledge of the spectrum's autocorrelation length
#'   (e.g. from the instrument's spectral resolution).
#' @param ratio A positive numeric convergence threshold.
#' @param max_iter A positive integer, maximum allowed iterations per replicate.
#' @param ridge A small positive numeric adjustment for linear system stability.
#' @param seed An optional integer used to seed the bootstrap resampling loop.
#'   The caller's global RNG state is saved before seeding and restored on
#'   exit, so supplying `seed` never affects random draws elsewhere in the
#'   caller's session.
#'
#' @return A validated object of class `StatARPLSResult` with moving-block
#' bootstrap stability bounds, and bootstrap diagnostics attached via
#' `attr(result, "bootstrap")`:
#' * `B_requested` — integer scalar, replicates requested
#' * `B_effective` — integer scalar, replicates that converged and were kept
#' * `n_failed` — integer scalar, replicates that errored or failed to converge
#' * `block_length` — integer scalar, the block length used
#' * `boot_se` — numeric vector, pointwise bootstrap standard errors
#'
#' @examples
#' set.seed(1)
#' wavenumber <- seq(400, 1800, length.out = 300)
#' peaks <- exp(-((wavenumber - 1000) / 25)^2) +
#'   0.6 * exp(-((wavenumber - 1400) / 30)^2)
#' spectrum <- peaks + 0.001 * (wavenumber - 400) + 0.2 +
#'   rnorm(300, sd = 0.01)
#'
#' fit <- fit_arpls(spectrum, wavenumber = wavenumber, lambda = 1e6)
#' boot_fit <- bootstrap_baseline_uncertainty(fit, B = 20, seed = 1)
#' attr(boot_fit, "bootstrap")$B_effective
#' plot(boot_fit)
#'
#' @export
bootstrap_baseline_uncertainty <- function(fit,
                                           B = 200L,
                                           block_length = NULL,
                                           ratio = 0.05,
                                           max_iter = 50L,
                                           ridge = 1e-6,
                                           seed = NULL) {
  if (!inherits(fit, "StatARPLSResult")) {
    stop("`fit` must be a StatARPLSResult object.", call. = FALSE)
  }
  if (!is.numeric(B) || length(B) != 1L || B < 2) {
    stop("`B` must be a single integer of at least 2.", call. = FALSE)
  }
  if (!is.numeric(ratio) || length(ratio) != 1L || ratio <= 0) {
    stop("`ratio` must be a single positive numeric value.", call. = FALSE)
  }
  if (!is.numeric(max_iter) || length(max_iter) != 1L || max_iter <= 0) {
    stop("`max_iter` must be a single integer greater than 0.", call. = FALSE)
  }
  if (!is.numeric(ridge) || length(ridge) != 1L || ridge < 0) {
    stop("`ridge` must be a single non-negative numeric value.", call. = FALSE)
  }
  if (!is.null(seed) && (!is.numeric(seed) || length(seed) != 1L)) {
    stop("`seed` must be a single integer, or NULL.", call. = FALSE)
  }

  n <- length(fit$signal)

  if (is.null(block_length)) {
    block_length <- max(2L, round(n^(1 / 3)))
    block_length <- min(block_length, max(2L, floor(n / 2)))
  } else {
    if (!is.numeric(block_length) || length(block_length) != 1L ||
        block_length < 1 || block_length > n) {
      stop(
        "`block_length` must be a single integer between 1 and ",
        "length(fit$signal) (", n, ").",
        call. = FALSE
      )
    }
  }
  block_length <- as.integer(block_length)

  if (!is.null(seed)) {
    has_old_seed <- exists(".Random.seed", envir = .GlobalEnv)
    old_seed <- if (has_old_seed) get(".Random.seed", envir = .GlobalEnv) else NULL
    on.exit({
      if (has_old_seed) {
        assign(".Random.seed", old_seed, envir = .GlobalEnv)
      } else if (exists(".Random.seed", envir = .GlobalEnv)) {
        rm(".Random.seed", envir = .GlobalEnv)
      }
    }, add = TRUE)
    set.seed(seed)
  }

  # Resample noise only: points the final weights classify as baseline.
  # Peak points would otherwise be injected into every pseudo-signal.
  residual_pool <- arpls_noise_residuals(fit$corrected)
  block_length <- min(block_length, length(residual_pool))
  boot_baselines <- matrix(NA_real_, nrow = n, ncol = B)
  n_failed <- 0L

  for (b in seq_len(B)) {
    resampled_residual <- moving_block_bootstrap_sample(
      residual_pool, block_length, n_out = n
    )
    pseudo_signal <- fit$baseline + resampled_residual

    engine_fit <- tryCatch(
      arpls_engine(pseudo_signal, fit$lambda, ratio, max_iter, ridge),
      error = function(e) NULL
    )

    if (is.null(engine_fit) || !isTRUE(engine_fit$converged)) {
      n_failed <- n_failed + 1L
      next
    }

    boot_baselines[, b] <- as.numeric(engine_fit$baseline)
  }

  keep <- !apply(boot_baselines, 2, anyNA)
  boot_baselines <- boot_baselines[, keep, drop = FALSE]

  if (ncol(boot_baselines) < 2) {
    stop(
      "Fewer than 2 bootstrap replicates converged (", n_failed,
      " of ", B, " failed).", call. = FALSE
    )
  }

  boot_se <- apply(boot_baselines, 1, stats::sd)
  bootstrap_quantiles <- t(apply(
    boot_baselines, 1, stats::quantile, probs = c(0.025, 0.975)
  ))
  # Centre the band on the fitted baseline. The asymmetric weighting makes
  # the replicate baselines sit slightly off the original fit; raw
  # percentiles would inherit that offset and often exclude the fit itself.
  bootstrap_ci <- fit$baseline + (bootstrap_quantiles - rowMeans(boot_baselines))

  updated_fit <- StatARPLSResult(
    signal            = fit$signal,
    wavenumber        = fit$wavenumber,
    baseline          = fit$baseline,
    corrected         = fit$corrected,
    residual_variance = fit$residual_variance,
    stability_bounds  = bootstrap_ci,
    lambda            = fit$lambda,
    n_iter            = fit$n_iter,
    converged         = fit$converged,
    call              = fit$call
  )

  attr(updated_fit, "gcv") <- attr(fit, "gcv")
  attr(updated_fit, "bootstrap") <- list(
    B_requested  = B,
    B_effective  = ncol(boot_baselines),
    n_failed     = n_failed,
    block_length = block_length,
    boot_se      = boot_se
  )

  updated_fit
}
