test_that("compare() generic dispatches correctly to compare.StatARPLSResult()", {
  spec <- make_synthetic_spectrum()
  fit_a <- fit_arpls(spec$signal, spec$wavenumber, lambda = 1e4)
  fit_b <- fit_arpls(spec$signal, spec$wavenumber, lambda = 1e6)

  result <- compare(fit_a, fit_b)

  expect_s3_class(result, "compare.StatARPLSResult")
  expect_equal(nrow(result), 2)
})

test_that("compare.StatARPLSResult() reports NA df_effective for non-GCV fits", {
  spec <- make_synthetic_spectrum()
  fit_a <- fit_arpls(spec$signal, spec$wavenumber, lambda = 1e4)
  fit_b <- fit_arpls(spec$signal, spec$wavenumber, lambda = 1e6)

  result <- compare(fit_a, fit_b)

  expect_true(all(is.na(result$df_effective)))
})

test_that("compare.StatARPLSResult() reports non-NA df_effective for GCV fits", {
  spec <- make_synthetic_spectrum()
  gcv_fit <- gcv_select_lambda(spec$signal, spec$wavenumber, lambda_grid = 10^seq(3, 6, length.out = 3))
  plain_fit <- fit_arpls(spec$signal, spec$wavenumber, lambda = 1e6)

  result <- compare(gcv_fit, plain_fit)

  expect_false(is.na(result$df_effective[result$label == "gcv_fit"]))
})

test_that("compare.StatARPLSResult() respects explicit labels argument", {
  spec <- make_synthetic_spectrum()
  fit_a <- fit_arpls(spec$signal, spec$wavenumber, lambda = 1e4)
  fit_b <- fit_arpls(spec$signal, spec$wavenumber, lambda = 1e6)

  result <- compare(fit_a, fit_b, labels = c("low", "high"))

  expect_equal(result$label, c("low", "high"))
})

test_that("compare.StatARPLSResult() errors with fewer than two fits", {
  spec <- make_synthetic_spectrum()
  fit_a <- fit_arpls(spec$signal, spec$wavenumber, lambda = 1e4)

  expect_error(compare(fit_a), regexp = "at least two")
})

test_that("compare.StatARPLSResult() errors on non-StatARPLSResult input in ...", {
  spec <- make_synthetic_spectrum()
  fit_a <- fit_arpls(spec$signal, spec$wavenumber, lambda = 1e4)

  # A non-StatARPLSResult object passed via `...` must be rejected for what
  # it actually is -- a type violation -- not silently dropped (which would
  # also corrupt row labels; see methods.R). This asserts on the cause, not
  # the count-based symptom the pre-fix implementation produced.
  expect_error(compare(fit_a, list(not = "a fit")), regexp = "must be of class")
})

test_that("print.compare.StatARPLSResult() runs without error", {
  spec <- make_synthetic_spectrum()
  fit_a <- fit_arpls(spec$signal, spec$wavenumber, lambda = 1e4)
  fit_b <- fit_arpls(spec$signal, spec$wavenumber, lambda = 1e6)

  expect_output(print(compare(fit_a, fit_b)), "compare.StatARPLSResult")
})
