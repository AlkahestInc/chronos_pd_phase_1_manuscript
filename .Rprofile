source("renv/activate.R")

get_config <- local({
  cfg <- NULL

  function() {
    if (is.null(cfg)) {
      cfg <<- config::get(file = "config.yaml")
    }
    cfg
  }
})

#####
# Specify the default logger layout
#####
logger::layout_glue_generator(format = "{level} [{time}] {fn}: {msg}") |>
  logger::log_layout()
logger::log_threshold(get_config()[["log_level"]])

#####
# Automatically load the packages in the "package" directory
#####
devtools::load_all("packages/manuscript")

#####
# Process the .origin.yml file for it's configuration options
#####
# .local_origin_config <- function(path = ".origin") {
#   if (!file.exists(path)) return(invisible(NULL))
#   if (!requireNamespace("yaml", quietly = TRUE)) return(invisible(NULL))
#
#   cfg <- yaml::read_yaml(path)
#   names(cfg) <- paste0("origin.", names(cfg))
#   options(cfg)
#
#   invisible(TRUE)
# }
# try(.local_origin_config(), silent = TRUE)
# rm(.local_origin_config)
# cfg <- yaml::read_yaml(path)
# names(cfg) <- paste0("origin.", names(cfg))
# options(cfg)
