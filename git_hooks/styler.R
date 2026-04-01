#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  if (!requireNamespace("styler", quietly = TRUE)) {
    message("Package 'styler' is required but not installed.")
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

# Set the styler cache to be permanant
# https://styler.r-lib.org/reference/caching.html
options(list(styler.cache_root = "styler-perm"))

tryCatch(
  {
    styler::style_file(files)
  },
  error = function(e) {
    message("Styling failed: ", conditionMessage(e))
    quit(status = 1)
  }
)

quit(status = 0)
