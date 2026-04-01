here::i_am("scripts/R/.Rprofile")

#####
# Specify the default logger layout
#####
logger::layout_glue_generator(format = "{level} [{time}] {fn}: {msg}") |>
  logger::log_layout()
logger::log_threshold(get_config()[["log_level"]])
