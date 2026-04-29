# Auto-sourced by testthat before running any test file.
# Since the package isn't installed (devtools dependencies missing on this
# system), we source the R/ sources directly so test functions can find them.

pkg_root <- normalizePath(file.path(testthat::test_path(), "..", ".."),
                          mustWork = TRUE)
r_files <- list.files(file.path(pkg_root, "R"), pattern = "\\.R$", full.names = TRUE)
for (f in r_files) source(f, local = FALSE)

# Load bundled dataset (not accessible via data() when package isn't installed)
prostate_rda <- file.path(pkg_root, "data", "prostate_data.rda")
if (file.exists(prostate_rda)) load(prostate_rda, envir = globalenv())
