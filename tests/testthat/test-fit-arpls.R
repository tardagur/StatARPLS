test_that("fit_arpls produces a structurally valid object on mock spectral data", {
  # 1. Simulate a synthetic spectroscopic signal with peaks
  set.seed(501)
  wavenumbers <- seq(400, 1800, length.out = 300)
  true_bg <- 10 + 0.005 * (wavenumbers - 400)
  peak_profile <- 20 * exp(-((wavenumbers - 900) / 25)^2)
  random_noise <- rnorm(300, mean = 0, sd = 0.4)
  mock_spectrum <- true_bg + peak_profile + random_noise

  # 2. Run the fitting algorithm operation
  res <- fit_arpls(signal = mock_spectrum, wavenumber = wavenumbers, lambda = 1e4, ratio = 1e-4)

  # 3. Structural assertions matching constructor contracts
  expect_s3_class(res, "StatARPLSResult")
  expect_length(res$baseline, 300)
  expect_true(all(is.finite(res$baseline)))
  expect_true(res$converged)
  expect_equal(dim(res$stability_bounds), c(300, 2))

  # 4. Assert parameter error trapping
  expect_error(
    fit_arpls(signal = mock_spectrum, wavenumber = wavenumbers[-1]),
    "wavenumber` must be the same length as `signal"
  )
})
