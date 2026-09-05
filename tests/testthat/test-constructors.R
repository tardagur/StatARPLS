test_that("StatARPLSResult() constructs a valid object with correct class and fields", {
  n <- 10
  result <- StatARPLSResult(
    signal            = rnorm(n),
    wavenumber        = seq_len(n),
    baseline          = rep(0, n),
    corrected         = rnorm(n),
    residual_variance = 1,
    stability_bounds  = cbind(rep(-1, n), rep(1, n)),
    lambda            = 100,
    n_iter            = 5L,
    converged         = TRUE,
    call              = quote(fit_arpls(x))
  )

  expect_s3_class(result, "StatARPLSResult")
  expect_named(
    result,
    c("signal", "wavenumber", "baseline", "corrected", "residual_variance",
      "stability_bounds", "lambda", "n_iter", "converged", "call")
  )
  expect_length(result$signal, n)
})

test_that("validate_StatARPLSResult() rejects mismatched vector lengths", {
  expect_error(
    StatARPLSResult(
      signal            = rnorm(10),
      wavenumber        = seq_len(9),          # deliberately wrong length
      baseline          = rep(0, 10),
      corrected         = rnorm(10),
      residual_variance = 1,
      stability_bounds  = cbind(rep(-1, 10), rep(1, 10)),
      lambda            = 100,
      n_iter            = 5L,
      converged         = TRUE,
      call              = quote(fit_arpls(x))
    ),
    regexp = "same length"
  )
})

test_that("validate_StatARPLSResult() rejects malformed stability_bounds", {
  n <- 10
  expect_error(
    StatARPLSResult(
      signal            = rnorm(n),
      wavenumber        = seq_len(n),
      baseline          = rep(0, n),
      corrected         = rnorm(n),
      residual_variance = 1,
      stability_bounds  = matrix(0, nrow = n, ncol = 3),  # wrong ncol
      lambda            = 100,
      n_iter            = 5L,
      converged         = TRUE,
      call              = quote(fit_arpls(x))
    ),
    regexp = "stability_bounds"
  )
})

test_that("validate_StatARPLSResult() rejects non-positive lambda", {
  n <- 5
  expect_error(
    StatARPLSResult(
      signal = rnorm(n), wavenumber = seq_len(n), baseline = rep(0, n),
      corrected = rnorm(n), residual_variance = 1,
      stability_bounds = cbind(rep(-1, n), rep(1, n)),
      lambda = -5, n_iter = 1L, converged = TRUE, call = quote(f())
    ),
    regexp = "lambda"
  )
})

test_that("validate_StatARPLSResult() rejects non-scalar converged", {
  n <- 5
  expect_error(
    StatARPLSResult(
      signal = rnorm(n), wavenumber = seq_len(n), baseline = rep(0, n),
      corrected = rnorm(n), residual_variance = 1,
      stability_bounds = cbind(rep(-1, n), rep(1, n)),
      lambda = 10, n_iter = 1L, converged = c(TRUE, FALSE), call = quote(f())
    ),
    regexp = "converged"
  )
})
