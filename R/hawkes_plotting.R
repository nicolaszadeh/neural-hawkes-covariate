# ============================================================
# Plotting utilities for the Hawkes project
# ============================================================


# ============================================================
# Make sure the parent folder of a file exists
# ============================================================

ensure_parent_dir <- function(filename) {
  folder <- dirname(filename)

  if (!dir.exists(folder)) {
    dir.create(
      folder,
      recursive = TRUE,
      showWarnings = FALSE
    )
  }

  invisible(TRUE)
}


# ============================================================
# Save a plot as PNG
# ============================================================

save_png <- function(
    filename,
    plot_function,
    width = 1000,
    height = 900,
    res = 120
) {
  ensure_parent_dir(filename)

  png(
    filename = filename,
    width = width,
    height = height,
    res = res
  )

  on.exit(dev.off())

  plot_function()

  invisible(filename)
}


# ============================================================
# Save a plot as PDF
# ============================================================

save_pdf <- function(
    filename,
    plot_function,
    width = 8,
    height = 7
) {
  ensure_parent_dir(filename)

  pdf(
    file = filename,
    width = width,
    height = height
  )

  on.exit(dev.off())

  plot_function()

  invisible(filename)
}


# ============================================================
# Save the same plot as both PNG and PDF
# ============================================================

save_plot_both <- function(
    file_stem,
    plot_function,
    png_width = 1000,
    png_height = 900,
    png_res = 120,
    pdf_width = 8,
    pdf_height = 7
) {
  png_file <- paste0(file_stem, ".png")
  pdf_file <- paste0(file_stem, ".pdf")

  save_png(
    filename = png_file,
    plot_function = plot_function,
    width = png_width,
    height = png_height,
    res = png_res
  )

  save_pdf(
    filename = pdf_file,
    plot_function = plot_function,
    width = pdf_width,
    height = pdf_height
  )

  invisible(c(
    png = png_file,
    pdf = pdf_file
  ))
}