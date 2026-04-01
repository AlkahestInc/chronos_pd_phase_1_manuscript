#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  if (!requireNamespace("lintr", quietly = TRUE)) {
    message("Package 'lintr' is required but not installed.")
    quit(status = 1)
  }
})


# The hook will get called for each changed file passed on the commandline
args <- commandArgs(trailingOnly = TRUE)

# Exit if no files were passed
if (length(args) == 0) {
  quit(status = 0)
}

# Only keep files that exist (deleted filenames may be passed)
files <- args[file.exists(args)]

# Exit if there are no files to process
if (length(files) == 0) {
  quit(status = 0)
}

tryCatch(
  {
    for (f in files) {
      lintr::lint(f)
    }
  },
  error = function(e) {
    message("Linting failed: ", conditionMessage(e))
    quit(status = 1)
  }
)

quit(status = 0)
