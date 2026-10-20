# data-raw/real-data-check-sre-ir.R
# Dev-only script (excluded from the build by `^data-raw$` in .Rbuildignore).
# Compares StatARPLS (automatic lambda) with the STATISTICAL baseline methods
# in the `baseline` package (als, irls, modpolyfit) on ChemoSpec::SrE.IR.
# `baseline` and `ChemoSpec` are NOT package dependencies; run interactively.
#
# Run from the package root:  devtools::load_all(); source("data-raw/real-data-check-sre-ir.R")

stopifnot(
  requireNamespace("ChemoSpec", quietly = TRUE),
  requireNamespace("baseline", quietly = TRUE)
)

dir.create("data-raw/output", showWarnings = FALSE, recursive = TRUE)

# ---- 1. Load and prepare the data --------------------------------------------
data("SrE.IR", package = "ChemoSpec")
wavenumber <- SrE.IR$freq
spectra <- SrE.IR$data            # rows = samples, columns = wavenumbers
sample_names <- SrE.IR$names

axis_order <- order(wavenumber)   # StatARPLS requires a strictly increasing axis
wavenumber <- wavenumber[axis_order]
spectra <- spectra[, axis_order, drop = FALSE]
stopifnot(!anyDuplicated(wavenumber))
cat("Spectra:", nrow(spectra), " points:", ncol(spectra), "\n")

# ---- 2. Real data: runtime, negativity and parameter sensitivity -------------
als_lambdas <- c(1e3, 1e5, 1e7)   # the hand-picked values a user must guess

# `baseline::baseline.als()` may take lambda on the log10 scale (it computes
# 10^lambda internally); passing 1e3 then overflows to Inf and fails with
# "NA/NaN/Inf in foreign function call". Detect this from the function body
# so that the lambda values above always mean the same thing.
als_uses_log10 <- any(grepl("10 *\\^", deparse(baseline::baseline.als)))
cat("baseline::baseline.als takes lambda on the log10 scale:", als_uses_log10, "\n")

fit_als <- function(y, lambda) {
  lambda_arg <- if (als_uses_log10) log10(lambda) else lambda
  as.numeric(baseline::getBaseline(
    baseline::baseline(matrix(y, nrow = 1), method = "als",
                       lambda = lambda_arg, p = 0.01, maxit = 50)
  ))
}

time_it <- function(expr) unname(system.time(expr)[["elapsed"]])

summarise_baseline <- function(signal, baseline_vector) {
  corrected <- signal - baseline_vector
  c(
    prop_negative = mean(corrected < 0),
    min_corrected = min(corrected)
  )
}

rows <- list()
for (i in seq_len(nrow(spectra))) {
  y <- spectra[i, ]

  t_stat <- time_it(fit <- gcv_select_lambda(y, wavenumber, seed = 1))
  rows[[length(rows) + 1]] <- data.frame(
    sample = sample_names[i], method = "StatARPLS (auto)", setting = signif(fit$lambda, 3),
    seconds = t_stat, t(summarise_baseline(y, fit$baseline))
  )

  for (lam in als_lambdas) {
    t_als <- time_it(als_fit <- fit_als(y, lam))
    rows[[length(rows) + 1]] <- data.frame(
      sample = sample_names[i], method = "AsLS (fixed)", setting = lam,
      seconds = t_als,
      t(summarise_baseline(y, als_fit))
    )
  }

  t_irls <- time_it(bl <- baseline::baseline(matrix(y, nrow = 1), method = "irls"))
  rows[[length(rows) + 1]] <- data.frame(
    sample = sample_names[i], method = "IRLS (default)", setting = NA_real_,
    seconds = t_irls,
    t(summarise_baseline(y, as.numeric(baseline::getBaseline(bl))))
  )

  t_poly <- time_it(bl <- baseline::baseline(matrix(y, nrow = 1), method = "modpolyfit", degree = 4))
  rows[[length(rows) + 1]] <- data.frame(
    sample = sample_names[i], method = "ModPolyFit (deg 4)", setting = 4,
    seconds = t_poly,
    t(summarise_baseline(y, as.numeric(baseline::getBaseline(bl))))
  )
}
real_results <- do.call(rbind, rows)
utils::write.csv(real_results, "data-raw/output/sre-ir-real-results.csv", row.names = FALSE)

cat("\nMedian runtime and negativity by method/setting (all spectra):\n")
print(aggregate(cbind(seconds, prop_negative) ~ method + setting, data = real_results, FUN = median))

# How much does the baseline move if the user guesses lambda differently?
# For AsLS: largest pointwise gap between lambda = 1e3, 1e5, 1e7, relative to signal range.
spread <- vapply(seq_len(nrow(spectra)), function(i) {
  y <- spectra[i, ]
  fits <- vapply(als_lambdas, function(lam) fit_als(y, lam), numeric(length(y)))
  max(apply(fits, 1, max) - apply(fits, 1, min)) / diff(range(y))
}, numeric(1))
cat("\nAsLS baseline spread across lambda (fraction of signal range): median",
    signif(stats::median(spread), 3), " max", signif(max(spread), 3), "\n")

# ---- 3. Figure and stability band for one spectrum ---------------------------
example_index <- 1L
y <- spectra[example_index, ]
fit <- gcv_select_lambda(y, wavenumber, seed = 1)
fit_boot <- bootstrap_baseline_uncertainty(fit, B = 50, seed = 2)

grDevices::png("data-raw/output/sre-ir-example.png", width = 1400, height = 900, res = 150)
graphics::plot(wavenumber, y, type = "l", col = "grey40",
               xlab = "Wavenumber (cm-1)", ylab = "Absorbance",
               main = paste("SrE.IR:", sample_names[example_index]))
graphics::polygon(c(wavenumber, rev(wavenumber)),
                  c(fit_boot$stability_bounds[, 1], rev(fit_boot$stability_bounds[, 2])),
                  col = grDevices::adjustcolor("steelblue", 0.35), border = NA)
graphics::lines(wavenumber, fit$baseline, col = "steelblue", lwd = 2)
cols <- c("firebrick", "darkorange", "forestgreen")
for (k in seq_along(als_lambdas)) {
  graphics::lines(wavenumber, fit_als(y, als_lambdas[k]), col = cols[k], lty = 2)
}
graphics::legend("topright", bty = "n", cex = 0.8,
                 legend = c("Signal", paste0("StatARPLS (lambda = ", signif(fit$lambda, 3), ", band)"),
                            paste0("AsLS lambda = ", format(als_lambdas, scientific = TRUE))),
                 col = c("grey40", "steelblue", cols), lty = c(1, 1, 2, 2, 2), lwd = c(1, 2, 1, 1, 1))
grDevices::dev.off()

# ---- 4. Simulation with a KNOWN baseline (the only true accuracy test) -------
# Peaks and baseline shapes mimic an FTIR spectrum on the same wavenumber grid.
simulate_known <- function(seed, noise_sd = 0.01) {
  set.seed(seed)
  x <- wavenumber
  scaled <- (x - min(x)) / diff(range(x))
  true_baseline <- 0.2 + 0.4 * scaled + 0.25 * sin(2.5 * pi * scaled + seed %% 3)
  peak_centres <- stats::runif(8, min(x) + 100, max(x) - 100)
  peaks <- Reduce(`+`, lapply(peak_centres, function(centre) {
    stats::runif(1, 0.1, 0.8) * exp(-((x - centre) / stats::runif(1, 15, 60))^2)
  }))
  list(signal = true_baseline + peaks + stats::rnorm(length(x), sd = noise_sd),
       truth = true_baseline)
}

rmse <- function(a, b) sqrt(mean((a - b)^2))
sim_rows <- list()
for (seed in 1:20) {
  sim <- simulate_known(seed)
  stat_fit <- gcv_select_lambda(sim$signal, wavenumber, seed = 1)
  out <- data.frame(seed = seed, method = "StatARPLS (auto)", rmse = rmse(stat_fit$baseline, sim$truth))
  for (lam in als_lambdas) {
    out <- rbind(out, data.frame(seed = seed, method = paste0("AsLS lambda=", format(lam, scientific = TRUE)),
                                 rmse = rmse(fit_als(sim$signal, lam), sim$truth)))
  }
  bl <- baseline::baseline(matrix(sim$signal, nrow = 1), method = "irls")
  out <- rbind(out, data.frame(seed = seed, method = "IRLS (default)",
                               rmse = rmse(as.numeric(baseline::getBaseline(bl)), sim$truth)))
  bl <- baseline::baseline(matrix(sim$signal, nrow = 1), method = "modpolyfit", degree = 4)
  out <- rbind(out, data.frame(seed = seed, method = "ModPolyFit (deg 4)",
                               rmse = rmse(as.numeric(baseline::getBaseline(bl)), sim$truth)))
  sim_rows[[seed]] <- out
}
sim_results <- do.call(rbind, sim_rows)
utils::write.csv(sim_results, "data-raw/output/sre-ir-simulation-rmse.csv", row.names = FALSE)

cat("\nSimulation (known baseline), baseline RMSE over 20 spectra:\n")
print(aggregate(rmse ~ method, data = sim_results,
                FUN = function(v) c(median = stats::median(v), mean = mean(v))))
cat("\nStatARPLS best in this many of 20 simulated spectra:",
    sum(tapply(sim_results$rmse, sim_results$seed, function(v) which.min(v) == 1L)), "\n")
