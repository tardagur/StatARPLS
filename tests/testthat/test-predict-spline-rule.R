# Regression tests for the peer-review finding that
# predict(method = "spline", rule = 2) ignored `rule`.

make_spline_rule_fit <- function() {
  set.seed(15)
  wavenumber <- seq(400, 1800, length.out = 80)
  signal <- 5 + 0.002 * wavenumber +
    3 * exp(-((wavenumber - 1000) / 40)^2) + stats::rnorm(80, sd = 0.5)
  fit_arpls(signal, wavenumber, lambda = 1e5)
}

test_that("predict() with method = 'spline' and rule = 2 clamps to the boundary values", {
  fit <- make_spline_rule_fit()
  below <- min(fit$wavenumber) - 1
  above <- max(fit$wavenumber) + 1

  predicted <- predict(fit, newdata = c(below, above), method = "spline", rule = 2)

  expect_equal(predicted[1], fit$baseline[1], tolerance = 1e-6)
  expect_equal(predicted[2], fit$baseline[length(fit$baseline)], tolerance = 1e-6)
})

test_that("predict() spline and linear agree on out-of-range points under rule = 2", {
  fit <- make_spline_rule_fit()
  outside <- c(min(fit$wavenumber) - 50, max(fit$wavenumber) + 50)

  expect_equal(
    predict(fit, newdata = outside, method = "spline", rule = 2),
    predict(fit, newdata = outside, method = "linear", rule = 2),
    tolerance = 1e-6
  )
})

test_that("predict() with method = 'spline' leaves in-range points unchanged by rule", {
  fit <- make_spline_rule_fit()
  inside <- fit$wavenumber[c(10, 40, 70)]

  expect_equal(
    predict(fit, newdata = inside, method = "spline", rule = 2),
    suppressWarnings(predict(fit, newdata = inside, method = "spline", rule = 1)),
    tolerance = 1e-8
  )
})

test_that("predict() with method = 'spline' and rule = 1 still warns and returns NA outside the range", {
  fit <- make_spline_rule_fit()
  below <- min(fit$wavenumber) - 1

  expect_warning(
    predicted <- predict(fit, newdata = below, method = "spline", rule = 1),
    "not supported"
  )
  expect_true(is.na(predicted))
})
