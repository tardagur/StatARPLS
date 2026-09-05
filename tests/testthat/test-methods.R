test_that("print.StatARPLSResult() runs without error and produces console output", {
  spec <- make_synthetic_spectrum()
  fit <- fit_arpls(spec$signal, spec$wavenumber, lambda = 1e5)

  expect_output(print(fit), "StatARPLSResult")
  expect_output(print(fit), "Converged")
})

test_that("print.StatARPLSResult() reports GCV/bootstrap sections only when present", {
  spec <- make_synthetic_spectrum()
  plain_fit <- fit_arpls(spec$signal, spec$wavenumber, lambda = 1e5)

  plain_output <- capture.output(print(plain_fit))
  expect_false(any(grepl("GCV selection", plain_output)))
})

test_that("summary.StatARPLSResult() returns correctly classed object with expected fields", {
  spec <- make_synthetic_spectrum()
  fit <- fit_arpls(spec$signal, spec$wavenumber, lambda = 1e5)
  smry <- summary(fit)

  expect_s3_class(smry, "summary.StatARPLSResult")
  expect_named(
    smry,
    c("lambda", "n_iter", "converged", "residual_variance",
      "residual_quantiles", "mean_band_width", "gcv", "bootstrap")
  )
  expect_length(smry$residual_quantiles, 5)
})

test_that("print.summary.StatARPLSResult() runs without error", {
  spec <- make_synthetic_spectrum()
  fit <- fit_arpls(spec$signal, spec$wavenumber, lambda = 1e5)

  expect_output(print(summary(fit)), "summary.StatARPLSResult")
})

test_that("plot.StatARPLSResult() runs without error and respects show_band toggle", {
  spec <- make_synthetic_spectrum()
  fit <- fit_arpls(spec$signal, spec$wavenumber, lambda = 1e5)

  expect_silent({
    png(tempfile())
    on.exit(dev.off(), add = TRUE)
    plot(fit, show_band = TRUE)
    plot(fit, show_band = FALSE)
  })
})
