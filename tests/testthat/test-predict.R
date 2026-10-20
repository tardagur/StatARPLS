test_that("predict.StatARPLSResult() recovers exact fitted values at original wavenumbers", {
  spec <- make_synthetic_spectrum()
  fit <- fit_arpls(spec$signal, spec$wavenumber, lambda = 1e5)

  predicted <- predict(fit, newdata = spec$wavenumber, method = "linear")

  expect_equal(predicted, fit$baseline, tolerance = 1e-8)
})

test_that("predict.StatARPLSResult() interpolates sensibly at a midpoint", {
  spec <- make_synthetic_spectrum()
  fit <- fit_arpls(spec$signal, spec$wavenumber, lambda = 1e5)

  midpoint <- (spec$wavenumber[1] + spec$wavenumber[2]) / 2
  predicted <- predict(fit, newdata = midpoint, method = "linear")

  expect_true(predicted >= min(fit$baseline[1:2]) && predicted <= max(fit$baseline[1:2]))
})

test_that("predict.StatARPLSResult() returns NA outside range under default rule = 1", {
  spec <- make_synthetic_spectrum()
  fit <- fit_arpls(spec$signal, spec$wavenumber, lambda = 1e5)

  out_of_range <- max(spec$wavenumber) + 1000
  predicted <- predict(fit, newdata = out_of_range)

  expect_true(is.na(predicted))
})

test_that("predict.StatARPLSResult() extrapolates as constant under rule = 2", {
  spec <- make_synthetic_spectrum()
  fit <- fit_arpls(spec$signal, spec$wavenumber, lambda = 1e5)

  out_of_range <- max(spec$wavenumber) + 1000
  predicted <- predict(fit, newdata = out_of_range, rule = 2)

  expect_equal(predicted, fit$baseline[length(fit$baseline)], tolerance = 1e-8)
})

test_that("predict.StatARPLSResult() spline method warns when rule = 1 is requested", {
  spec <- make_synthetic_spectrum()
  fit <- fit_arpls(spec$signal, spec$wavenumber, lambda = 1e5)

  out_of_range <- max(spec$wavenumber) + 1000
  expect_warning(
    predict(fit, newdata = out_of_range, method = "spline", rule = 1),
    regexp = "spline"
  )
})

test_that("predict.StatARPLSResult() rejects non-numeric newdata", {
  spec <- make_synthetic_spectrum()
  fit <- fit_arpls(spec$signal, spec$wavenumber, lambda = 1e5)

  expect_error(predict(fit, newdata = "not numeric"), regexp = "numeric")
})
