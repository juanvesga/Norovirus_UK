# Check packages ----------------------------------------------------------


required <- c(
  "here",
  "qs2",
  "readxl",
  "data.table",
  "gridExtra",
  "odin2",
  "dust2",
  "monty",  
  "tools",
  "matrixStats",
  "lubridate",
  "ISOweek",
  "tidyr",
  "dplyr",
  "ggplot2")

installed <- rownames(installed.packages())
not_installed <- required[!required %in% installed]

if (length(not_installed)) {
  msg <- "Some required packages are missing. Please install the following:\n"
  pkgs <- paste("\t", not_installed, collapse = "\n")
  msg <- sprintf("%s%s", msg, pkgs)
  stop(msg)
} else {
  # Set the root directory --------------------------------------------------
  root <- here::here()
  
  # List the fils we need to run in order
  files <- file.path(
    root, "R",
    c(
      "01-parameters.R",
      "02-data-calibration.R",
      "03-data-polymod-short.R",
      "04-data-school-uk.R",
      "05-data-comix-short.R",
      "06-data-covid.R",
      "07-generate-input-parameters.R",
      "08-generate-priors.R",
      "09-run-mcmc.R",
      "10-mcmc-dx.R",
      "11-run-mcmc-fits.R",
      "12-plot-mcmc-fits.R",
      "13-model-selection.R",
      "14-prepare-files-DIC.R",
      "15-life-events.R",
      "16-run-r0-estimates.R",
      "17-run-emergence.R",
      "18-plot-emergence.R",
      "19-table-emergence.R"
    )
  )
  
  # Loop over the files
  for (f in files) {
    local(source(f, local = TRUE))
  }
  
  message("Analysis complete. Results are in the 'results/' directory.")
}



