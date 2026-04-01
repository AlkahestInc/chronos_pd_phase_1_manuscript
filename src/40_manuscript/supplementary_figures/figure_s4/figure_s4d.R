## SPDX-FileCopyrightText: 2026 Alkahest, Inc. <chronospd_data_access@alkahest.com>
##
## SPDX-License-Identifier: MIT

here::i_am("src/40_manuscript/supplementary_figures/figure_s4/figure_s4d.R")

# Define input files
## protein files:
file_proteomics_somalogic <- glue::glue("{get_config()$input_dir}/name=Proteomics_Somalogic/version=none/dataset.parquet")
file_annotations_somalogic <- glue::glue("{get_config()$input_dir}/name=Annotations_Somalogic/version=none/dataset.parquet")
file_proteomics_alamarbio <- glue::glue("{get_config()$input_dir}/name=Proteomics_AlamarBio/version=none/dataset.parquet")
file_annotations_alamarbio <- glue::glue("{get_config()$input_dir}/name=Annotations_AlamarBio/version=none/dataset.parquet")
file_proteomics_olink <- glue::glue("{get_config()$input_dir}/name=Proteomics_Olink/version=none/dataset.parquet")
file_annotations_olink <- glue::glue("{get_config()$input_dir}/name=Annotations_Olink/version=none/dataset.parquet")

## sample and donor information files:
file_grifols_donor_pd <- glue::glue("{get_config()$input_dir}/name=Grifols_Donor_PD/version=none/dataset.parquet")
file_grifols_demographics_donor <- glue::glue("{get_config()$input_dir}/name=Grifols_Demographics_Donor/version=none/dataset.parquet")
file_grifols_sample_level <-  glue::glue("{get_config()$input_dir}/name=Grifols_Sample_Level/version=none/dataset.parquet")

## RWD files:
file_rx_kythera <- glue::glue("{get_config()$input_dir}/name=RX_Kythera/version=none/dataset.parquet")
file_ldopa_medication_ndcs <- glue::glue("{get_config()$input_dir}/name=LDOPA_medication_NDCs/version=none/dataset.parquet")

# Define output file
output_figure_file <- glue::glue("{get_config()$output_dir}/supplementary_figures/figure_s4/figure_s4d.pdf")
dir.create(dirname(output_figure_file), showWarnings = FALSE, recursive = TRUE)

# list of possibly levodopa-related protein_ids (five proteins including SUMF1, OIP5, TUBA1A, GPT, and DDC showed minimal change before onset and marked shifts afterward):
olink_hits <- c("OID44598", "OID43147") # DDC, SUMF1
soma_hits <- c("seq.16015.19", "seq.3709.4", "seq.33162.51", "seq.6941.11", "seq.29387.44", "seq.3538.26") # GPT, GPT, TUBA1A, SUMF1, OIP5, DDC
alamarbio_hits <- c("DDC")
hits <- c(olink_hits, soma_hits, alamarbio_hits)

# medication data of the donors
rx_kythera <- arrow::read_parquet(file_rx_kythera) |>
  dplyr::select(patient_id, date_of_service, national_drug_code) # medication data

# raw-protein input
protein_input_somalogic <- arrow::read_parquet(file_proteomics_somalogic)
protein_input_olink <- arrow::read_parquet(file_proteomics_olink)
protein_input_alamarbio <- arrow::read_parquet(file_proteomics_alamarbio)

# probe annotations
protein_annotation <- dplyr::bind_rows(
  arrow::read_parquet(file_annotations_somalogic) |>
    dplyr::select(ID, EntrezGeneSymbol) |>
    dplyr::rename(gene_symbol = EntrezGeneSymbol) |>
    dplyr::mutate(source = "SomaLogic 11k"),
  arrow::read_parquet(file_annotations_olink) |>
    dplyr::select(ID, ASSAY) |>
    dplyr::rename(gene_symbol = ASSAY) |>
    dplyr::mutate(source = "Olink 5k"),
  arrow::read_parquet(file_annotations_alamarbio) |>
    dplyr::select(ID, Target) |>
    dplyr::rename(gene_symbol = Target) |>
    dplyr::mutate(source = "Alamar CNS 120")) |>
  dplyr::rename(protein_id = ID)

# merging donor-sample level information
grifols_donor_pd <- arrow::read_parquet(file_grifols_donor_pd) |>
  dplyr::select(
    patient_id,
    estimated_pd_onset_date,
    age_at_baseline,
    group
  )
grifols_demographics_donor <- arrow::read_parquet(file_grifols_demographics_donor) |>
  dplyr::select(
    patient_id,
    gender
  ) |>
  dplyr::distinct()
grifols_sample_level <- arrow::read_parquet(file_grifols_sample_level) |>
  dplyr::select(
    patient_id,
    sample_collection_date,
    recent_plasma_donations,
    cluster,
    years_since_pd_onset,
    aliquot_barcode
  )
grifols_sample_input <- grifols_sample_level |>
  dplyr::inner_join(
    grifols_donor_pd,
    by = "patient_id"
  ) |>
  dplyr::inner_join(
    grifols_demographics_donor,
    by = "patient_id"
  )

# list of levodopa medication (NDCs)
pd_meds <- arrow::read_parquet(file_ldopa_medication_ndcs) |>
  dplyr::select(NDCPACKAGECODE_11) |>
  dplyr::rename(national_drug_code = NDCPACKAGECODE_11) |>
  dplyr::distinct() |>
  dplyr::mutate(MED_GROUP = "PD_incl")

# phase1 PD-cohort with onset
phase1_donors <- grifols_sample_input |>
  dplyr::filter(group == "pd") |>
  dplyr::select(patient_id, estimated_pd_onset_date) |>
  dplyr::distinct()

# collecting medication
patient_meds <- phase1_donors |>
  dplyr::inner_join(
    rx_kythera,
    by = dplyr::join_by("patient_id")
  ) |>
  dplyr::select(patient_id, national_drug_code, date_of_service, estimated_pd_onset_date) |>
  dplyr::filter(!(is.na(national_drug_code))) |>
  dplyr::filter(!(is.na(date_of_service))) |>
  dplyr::distinct()

# calculate levodopa effect for the proteomic platforms
lmm_results_somalogic <- levodopa_processing(
  protein_input_somalogic,
  grifols_sample_input |>
    dplyr::filter(aliquot_barcode %in% protein_input_somalogic$aliquot_barcode),
  protein_annotation |>
    dplyr::filter(source == "SomaLogic 11k"),
  patient_meds,
  phase1_donors,
  pd_meds
)
lmm_results_olink <- levodopa_processing(
  protein_input_olink,
  grifols_sample_input |>
    dplyr::filter(aliquot_barcode %in% protein_input_olink$aliquot_barcode),
  protein_annotation |>
    dplyr::filter(source == "Olink 5k"),
  patient_meds,
  phase1_donors,
  pd_meds
)
lmm_results_alamarbio <- levodopa_processing(
  protein_input_alamarbio,
  grifols_sample_input |>
    dplyr::filter(aliquot_barcode %in% protein_input_alamarbio$aliquot_barcode),
  protein_annotation |>
    dplyr::filter(source == "Alamar CNS 120"),
  patient_meds,
  phase1_donors,
  pd_meds
)

# merge LMM tables
data <- dplyr::bind_rows(
  lmm_results_alamarbio |>
    dplyr::mutate(source = "Alamar CNS 120"),
  lmm_results_olink |>
    dplyr::mutate(source = "Olink 5k"),
  lmm_results_somalogic |>
    dplyr::mutate(source = "SomaLogic 11k"))

lmm_results_plot <- data |>
  dplyr::mutate(
    CohensD_group = dplyr::case_when(
      cohensD >= 0.8 ~ "Cohen's D ≥ 0.8",
      cohensD <= -0.8 ~ "Cohen's D ≤ -0.8",
      TRUE ~ "|Cohen's D| < 0.8"
    ),
    is_med_affected = protein_id %in% hits # & p.value <= 0.05
  )

# volcano plot parameters
threshold_005 <- -log10(0.05)
color_values = c(
  "|Cohen's D| < 0.8" = "grey",
  "Cohen's D ≤ -0.8" = "#1f77b4",
  "Cohen's D ≥ 0.8" = "#d62728"
)
x_label = "log2FC between levodopa+ and levodopa- PDs"
y_label = "-log10(p.value)"
point_size = 0.7
text_size = 2.5
axis_text_size = 13
axis_title_size = 13
legend_title_size = 13
legend_text_size = 11
strip_text_size = 15
max_overlaps = 10
p_size = 4
gene_symbol_size = 5


# volcano plot
p <- lmm_results_plot |>
  dplyr::mutate(
    CohensD_group = factor(
      CohensD_group,
      levels = c("|Cohen's D| < 0.8", "Cohen's D ≤ -0.8", "Cohen's D ≥ 0.8")
    )
  ) |>
  dplyr::arrange(CohensD_group) |>
  ggplot2::ggplot(
    ggplot2::aes(
      x = estimate,
      y = neg_log10_p,
      colour = CohensD_group
    )
  ) +
  ggplot2::facet_wrap(~source, nrow = 1) +
  ggplot2::geom_point(size = point_size) +
  ggplot2::scale_colour_manual(values = color_values) +
  ggplot2::geom_hline(
    yintercept = threshold_005,
    linetype = "dashed",
    color = "darkgreen"
  ) +
  ggplot2::annotate(
    "text",
    x = -Inf, y = threshold_005,
    label = "p = 0.05",
    hjust = -0.1, vjust = -0.4,
    color = "darkgreen",
    size = p_size
  ) +
  ggplot2::labs(
    x = x_label,
    y = y_label,
    color = "Effect size (Cohen's D)"
  ) +
  ggplot2::theme_minimal() +
  ggplot2::theme(
    axis.text = ggplot2::element_text(size = axis_text_size),
    axis.title = ggplot2::element_text(size = axis_title_size, face = "bold"),
    strip.text = ggplot2::element_text(size = strip_text_size, face = "bold"),
    legend.title = ggplot2::element_text(size = legend_title_size, face = "bold"),
    legend.text  = ggplot2::element_text(size = legend_text_size)
  ) +
  ggrepel::geom_label_repel(
    data = lmm_results_plot |>
      dplyr::filter(is_med_affected),
    ggplot2::aes(
      label = gene_symbol
    ),
    fill = scales::alpha("white", 0.6),
    color = "black",
    size = gene_symbol_size,
    label.size = NA,
    box.padding = 0.2,
    min.segment.length = 0.0,
    max.overlaps = Inf,
    seed = 1234
  )

ggplot2::ggsave(
  filename = output_figure_file,
  plot = p,
  device = cairo_pdf,
  width = 10,
  height = 4,
  units = "in",
  dpi = 600
)
