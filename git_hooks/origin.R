#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  if (!requireNamespace("origin", quietly = TRUE)) {
    message("Package 'origin' is required but not installed.")
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

# Set the origin options according to the repo's config file
cfg <- yaml::read_yaml(".origin")

cfg[["overwrite"]] <- TRUE
cfg[["ask_before_applying_changes"]] <- FALSE
cfg[["verbose"]] <- FALSE
cfg[["use_markers_for_logging"]] <- FALSE
# cfg[["check_local_conflicts"]] <- FALSE
# cfg[["path_to_local_functions"]] <- "./"

names(cfg) <- paste0("origin.", names(cfg))

options(cfg)

tryCatch(
  {
    for (f in files) {
      suppressMessages(origin::originize_file(f))
    }
  },
  error = function(e) {
    message("Originizing failed: ", conditionMessage(e))
    quit(status = 1)
  }
)

quit(status = 0)
