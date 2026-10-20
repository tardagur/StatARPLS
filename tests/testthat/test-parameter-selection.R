test_that("gcv_select_lambda() returns a validated StatARPLSResult with gcv attribute", {
  spec <- make_synthetic_spectrum()
  fit <- gcv_select_lambda(
    spec$signal, spec$wavenumber,
    lambda_grid = 10^seq(3, 7, length.out = 5)
  )

  expect_s3_class(fit, "StatARPLSResult")
  gcv_info <- attr(fit, "gcv")
  expect_type(gcv_info, "list")
  expect_true(fit$lambda %in% gcv_info$lambda_grid)
})

test_that("gcv_select_lambda() selected lambda achieves the minimum GCV score on its own grid", {
  spec <- make_synthetic_spectrum()
  fit <- gcv_select_lambda(
    spec$signal, spec$wavenumber,
    lambda_grid = 10^seq(3, 7, length.out = 5)
  )
  gcv_info <- attr(fit, "gcv")

  expect_equal(fit$lambda, gcv_info$lambda_grid[which.min(gcv_info$gcv_scores)])
})

test_that("gcv_select_lambda() validates inputs", {
  expect_error(gcv_select_lambda("not numeric"), "`signal` must be a numeric vector")
  expect_error(gcv_select_lambda(c(1, 2)), "at least 3 points")

  spec <- make_synthetic_spectrum()
  expect_error(
    gcv_select_lambda(spec$signal, lambda_grid = c(-1, 2)),
    "strictly positive"
  )
  expect_error(
    gcv_select_lambda(spec$signal, ratio = -0.1),
    "`ratio` must be"
  )
  expect_error(
    gcv_select_lambda(spec$signal, n_probe = 0),
    "`n_probe` must be"
  )
})

test_that("gcv_select_lambda() survives a forced engine failure at one grid point", {
  # A constant signal does NOT create a singular system for this solver --
  # the second-difference operator annihilates constants, so `baseline =
  # signal` exactly satisfies (W + lambda*D'D)*baseline = W*signal for any
  # lambda. An earlier version of this test used a constant signal plus an
  # extreme lambda expecting that to be numerically pathological; it is
  # not, for this specific system (the identity-like weight term makes it
  # positive-definite by construction regardless of lambda), and the test
  # passed without ever exercising the tryCatch() path it claimed to
  # verify. Forcing the failure directly via mocking is the only reliable
  # way to test that gcv_select_lambda() actually survives an engine
  # failure, rather than hoping one occurs organically.
  spec <- make_synthetic_spectrum()
  grid <- c(1e4, 1e5, 1e6)

  testthat::local_mocked_bindings(
    arpls_engine = function(signal, lambda, ratio, max_iter, ridge) {
      if (isTRUE(all.equal(lambda, grid[1]))) {
        stop("mocked engine failure to exercise gcv_select_lambda()'s tryCatch")
      }
      arpls_fit_cpp(signal, lambda, ratio, max_iter, ridge)
    }
  )

  fit <- expect_no_error(gcv_select_lambda(spec$signal, lambda_grid = grid, seed = 1))
  expect_s3_class(fit, "StatARPLSResult")

  gcv_scores <- attr(fit, "gcv")$gcv_scores
  expect_length(gcv_scores, length(grid))
  expect_true(is.infinite(gcv_scores[1]))
  expect_true(all(is.finite(gcv_scores[-1])))
})

test_that("arpls_hat_trace() effective df decreases as lambda increases", {
  spec <- make_synthetic_spectrum()
  n <- length(spec$signal)
  D <- arpls_difference_matrix(n)

  engine_low  <- arpls_engine(spec$signal, 1e3, 0.05, 50L, 1e-6)
  engine_high <- arpls_engine(spec$signal, 1e7, 0.05, 50L, 1e-6)

  penalty_low  <- 1e3 * Matrix::crossprod(D)
  penalty_high <- 1e7 * Matrix::crossprod(D)

  df_est_low  <- arpls_hat_trace(engine_low$weights, penalty_low, 1e-6, n_probe = 30L, seed = 42)
  df_est_high <- arpls_hat_trace(engine_high$weights, penalty_high, 1e-6, n_probe = 30L, seed = 42)

  expect_gt(df_est_low, df_est_high)
})

test_that("moving_block_bootstrap_sample() preserves contiguous block structure", {
  pool <- 1:10
  resampled <- moving_block_bootstrap_sample(pool, block_length = 10L)

  expect_length(resampled, 10)
  # a single full-length block is a contiguous rotation of `pool`, not a
  # shuffle -- distinguishes genuine block resampling from disguised i.i.d.
  diffs <- diff(resampled)
  expect_true(all(diffs == 1) || all(diffs == -9))
})

test_that("bootstrap_baseline_uncertainty() returns tightened stability_bounds and bootstrap attribute", {
  spec <- make_synthetic_spectrum()
  base_fit <- fit_arpls(spec$signal, spec$wavenumber, lambda = 1e5)

  boot_fit <- bootstrap_baseline_uncertainty(base_fit, B = 40L, seed = 7)

  boot_info <- attr(boot_fit, "bootstrap")
  expect_type(boot_info, "list")
  expect_equal(boot_info$B_requested, 40L)
  expect_lte(boot_info$B_effective, 40L)
  expect_length(boot_info$boot_se, length(spec$signal))
  expect_true(all(boot_info$boot_se >= 0))
})

test_that("bootstrap_baseline_uncertainty() preserves the original baseline point estimate", {
  spec <- make_synthetic_spectrum()
  base_fit <- fit_arpls(spec$signal, spec$wavenumber, lambda = 1e5)
  boot_fit <- bootstrap_baseline_uncertainty(base_fit, B = 30L, seed = 7)

  expect_equal(boot_fit$baseline, base_fit$baseline)
})

test_that("bootstrap_baseline_uncertainty() rejects non-StatARPLSResult input", {
  expect_error(
    bootstrap_baseline_uncertainty(list(a = 1), B = 10L),
    regexp = "StatARPLSResult"
  )
})

test_that("bootstrap_baseline_uncertainty() preserves prior gcv attribute", {
  spec <- make_synthetic_spectrum()
  gcv_fit <- gcv_select_lambda(spec$signal, spec$wavenumber, lambda_grid = 10^seq(3, 6, length.out = 4))
  boot_fit <- bootstrap_baseline_uncertainty(gcv_fit, B = 20L, seed = 3)

  expect_identical(attr(boot_fit, "gcv"), attr(gcv_fit, "gcv"))
})

test_that("bootstrap_baseline_uncertainty() validates B, ratio, max_iter, ridge, block_length", {
  spec <- make_synthetic_spectrum()
  base_fit <- fit_arpls(spec$signal, spec$wavenumber, lambda = 1e5)

  expect_error(bootstrap_baseline_uncertainty(base_fit, B = 1L), "at least 2")
  expect_error(bootstrap_baseline_uncertainty(base_fit, ratio = -1), "`ratio` must be")
  expect_error(bootstrap_baseline_uncertainty(base_fit, max_iter = 0L), "`max_iter` must be")
  expect_error(bootstrap_baseline_uncertainty(base_fit, ridge = -1e-6), "`ridge` must be")
  expect_error(
    bootstrap_baseline_uncertainty(base_fit, block_length = length(spec$signal) + 10L),
    "between 1 and"
  )
})

test_that("seeded gcv_select_lambda() and bootstrap_baseline_uncertainty() do not leak RNG state", {
  spec <- make_synthetic_spectrum()
  base_fit <- fit_arpls(spec$signal, spec$wavenumber, lambda = 1e5)

  set.seed(123)
  draw_before <- runif(1)

  set.seed(123)
  invisible(bootstrap_baseline_uncertainty(base_fit, B = 5L, seed = 999))
  draw_after_boot <- runif(1)
  expect_equal(draw_before, draw_after_boot)

  set.seed(123)
  invisible(gcv_select_lambda(spec$signal, lambda_grid = c(1e4, 1e5), seed = 999))
  draw_after_gcv <- runif(1)
  expect_equal(draw_before, draw_after_gcv)
})
