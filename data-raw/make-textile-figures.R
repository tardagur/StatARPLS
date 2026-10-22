# data-raw/make-textile-figures.R
# Dev-only script (excluded from the build by `^data-raw$` in .Rbuildignore).
# Rebuilds the four-fibre textile example and saves the figures used in the
# presentation. Data: Mendeley textile ATR-FTIR, 02_feature_matrix_raw (CC BY 4.0).
#
# Run from the package root:
#   devtools::load_all()
#   source("data-raw/make-textile-figures.R")
# A file chooser opens: select the CSV in 02_feature_matrix_raw.

raw_path <- file.choose()
raw <- utils::read.csv(raw_path, check.names = FALSE)

# ---- 1. Split metadata and spectra; transmission -> absorbance ----------------
meta <- raw[, 1:6]
wavenumber <- as.numeric(names(raw)[-(1:6)])
transmission <- as.matrix(raw[, -(1:6)])

scale_max <- if (max(transmission, na.rm = TRUE) > 1.5) 100 else 1
absorbance <- -log10(pmax(transmission / scale_max, 1e-4))

axis_order <- order(wavenumber)
wavenumber <- wavenumber[axis_order]
absorbance <- absorbance[, axis_order, drop = FALSE]
keep <- seq(1, length(wavenumber), by = 2)
wavenumber <- wavenumber[keep]
absorbance <- absorbance[, keep, drop = FALSE]
stopifnot(!anyDuplicated(wavenumber))
cat("Axis points:", length(wavenumber), "\n")   # expect 1726

# ---- 2. Four fibres, automatic lambda on the safer (larger-minimum) grid -------
chosen <- c("Protein", "Cellulose", "Polyester", "Acrylic")
idx <- vapply(chosen, function(type) which(meta$Type == type)[1], integer(1))
print(meta[idx, ])   # which spectra are shown (first of each type); cite in the report
safe_grid <- 10^seq(4, 9, length.out = 11)

fits_safe <- lapply(idx, function(i)
  gcv_select_lambda(absorbance[i, ], wavenumber, lambda_grid = safe_grid, seed = 1))
boots_safe <- lapply(fits_safe, function(f)
  bootstrap_baseline_uncertainty(f, B = 50, seed = 2))

print(data.frame(
  type   = chosen,
  lambda = vapply(fits_safe, function(f) f$lambda, numeric(1)),
  df     = vapply(fits_safe, function(f) attr(f, "gcv")$df_effective, numeric(1))
))

# ---- 3. Figures ---------------------------------------------------------------
grDevices::graphics.off()   # avoids "invalid graphics state" from the Plots pane
dir.create("presentation/figures", showWarnings = FALSE, recursive = TRUE)

draw_panels <- function(fit_list, show_band, use_ylim = FALSE) {
  op <- graphics::par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
  on.exit(graphics::par(op))
  for (k in seq_along(fit_list)) {
    panel_title <- paste0(chosen[k], " (lambda = ",
                          format(fit_list[[k]]$lambda, digits = 2, scientific = TRUE), ")")
    if (use_ylim) {
      plot(fit_list[[k]], show_band = show_band, main = panel_title,
           ylab = "Absorbance",
           ylim = c(-0.01, max(absorbance[idx[k], ]) * 1.05))
    } else {
      plot(fit_list[[k]], show_band = show_band, main = panel_title,
           ylab = "Absorbance")
    }
  }
}

grDevices::png("presentation/figures/textile-fits.png", width = 2000, height = 1400, res = 200)
draw_panels(fits_safe, show_band = FALSE)
grDevices::dev.off()

grDevices::png("presentation/figures/textile-bands.png", width = 2000, height = 1400, res = 200)
draw_panels(boots_safe, show_band = TRUE, use_ylim = TRUE)
grDevices::dev.off()

cat("Saved presentation/figures/textile-fits.png and textile-bands.png\n")
