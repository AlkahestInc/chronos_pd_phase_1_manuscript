## SPDX-FileCopyrightText: 2026 Alkahest, Inc. <chronospd_data_access@alkahest.com>
##
## SPDX-License-Identifier: MIT

ukb_rap_load_proteomics_data <- function(rap_path) {
	to_file <- tempfile(fileext = ".rds")

	ukb_rap_download(
		from = rap_path,
		to = to,
		overwrite = TRUE
	)

	load(to, envir = environment())

	list(
		data = dataProt,
		row_metadata = dataOther,
		col_metadata = NA,
		participant_metadata = dons
	)
}
