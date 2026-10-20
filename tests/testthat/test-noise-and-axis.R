make_peaky_spectrum <- function(n = 400, noise_sd = 2, seed = 11) {
  set.seed(seed)
  wavenumber <- seq(400, 2000, length.out = n)
  baseline <- 40 + 0.02 * (wavenumber - 400)
  peaks <- 400 * exp(-((wavenumber - 900) / 20)^2) +
    300 * exp(-((wavenumber - 1500) / 30)^2)
  list(
    signal = baseline + peaks + rnorm(n, sd = noise_sd),
    wavenumber = wavenumber,
    noise_sd = noise_sd
  )
}

test_that("residual_variance reflects noise, not peaks", {
  s <- make_peaky_spectrum()
  fit <- gcv_select_lambda(s$signal, s$wavenumber, seed = 1)
  expect_lt(abs(sqrt(fit$residual_variance) - s$noise_sd) / s$noise_sd, 0.35)
  expect_lt(sqrt(fit$residual_variance), 0.2 * stats::sd(fit$corrected))
})

test_that("provisional stability band width is set by noise SD", {
  s <- make_peaky_spectrum()
  fit <- gcv_select_lambda(s$signal, s$wavenumber, seed = 1)
  half_width <- (fit$stability_bounds[, 2] - fit$stability_bounds[, 1]) / 2
  expect_equal(half_width, rep(1.96 * sqrt(fit$residual_variance), length(half_width)))
})

test_that("bootstrap band is centred on the fitted baseline and narrow", {
  s <- make_peaky_spectrum()
  fit <- gcv_select_lambda(s$signal, s$wavenumber, seed = 1)
  boot <- bootstrap_baseline_uncertainty(fit, B = 30, seed = 2)
  b <- boot$stability_bounds
  expect_true(all(b[, 1] <= boot$baseline & boot$baseline <= b[, 2]))
  expect_lt(stats::median(b[, 2] - b[, 1]), 4 * 1.96 * sqrt(fit$residual_variance))
})

test_that("moving_block_bootstrap_sample honours n_out with a short pool", {
  set.seed(1)
  out <- moving_block_bootstrap_sample(c(1, 2, 3, 4, 5), block_length = 2L, n_out = 12L)
  expect_length(out, 12L)
  expect_true(all(out %in% 1:5))
})

test_that("arpls_noise_residuals falls back with a warning when few baseline points", {
  expect_warning(
    out <- arpls_noise_residuals(c(-1, -1.2, 500, 600, 700, 800)),
    "Fewer than 3 points"
  )
  expect_length(out, 6L)
})

test_that("wavenumber validation rejects bad axes with clear messages", {
  y <- c(1, 2, 3, 4, 5, 6)
  expect_error(validate_wavenumber_axis(c(1, 2, 2, 3, 4, 5)), "strictly increasing")
  expect_error(validate_wavenumber_axis(c(6, 5, 4, 3, 2, 1)), "Sort `signal`")
  expect_error(validate_wavenumber_axis(c(1, 2, NA, 4, 5, 6)), "non-finite")
  expect_warning(validate_wavenumber_axis(c(1, 2, 3, 10, 11, 12)), "unevenly spaced")
  expect_silent(validate_wavenumber_axis(seq(400, 4000, length.out = 50)))
  expect_error(gcv_select_lambda(y, c(1, 2, 2, 3, 4, 5)), "strictly increasing")
})

test_that("fit_arpls rejects duplicated or descending wavenumbers", {
  s <- make_peaky_spectrum(n = 100)
  expect_error(fit_arpls(s$signal, rev(s$wavenumber), lambda = 1e5), "strictly increasing")
  wn <- s$wavenumber
  wn[50] <- wn[49]
  expect_error(fit_arpls(s$signal, wn, lambda = 1e5), "strictly increasing")
})

test_that("a lambda selected at the lower grid edge triggers a warning", {
  expect_warning(
    warn_if_lambda_at_grid_edge(10^(2:6), 1e2),
    "smallest value in `lambda_grid`"
  )
  expect_no_warning(warn_if_lambda_at_grid_edge(10^(2:6), 1e4))
  expect_no_warning(warn_if_lambda_at_grid_edge(10^(2:6), 1e6))
  expect_no_warning(warn_if_lambda_at_grid_edge(c(1e3, 1e5), 1e3))
})
