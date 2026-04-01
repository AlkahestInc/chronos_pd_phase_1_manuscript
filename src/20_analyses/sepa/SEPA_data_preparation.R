## SPDX-FileCopyrightText: 2026 Alkahest, Inc. <chronospd_data_access@alkahest.com>
##
## SPDX-License-Identifier: MIT

here::i_am("src/20_analyses/sepa/SEPA_data_preparation.R")

# ----------------------------------------------------------------------------
# Load Chronos results
# ----------------------------------------------------------------------------
raw_datasets <- list()

# Somalogic 11k linear mixed model rerun results
raw_datasets$chronos_somalogic_all <- glue::glue("{get_config()$input_dir}/name=chronos_somalogic_all-pd_master_summary/version=29806/dataset.parquet") |>
  arrow::open_dataset() |>
  dplyr::collect() |>
  wide_to_long_probe_map()

raw_datasets$chronos_somalogic_pre <- glue::glue("{get_config()$input_dir}/name=chronos_somalogic_pre_onset-pd_master_summary/version=29807/dataset.parquet") |>
  arrow::open_dataset() |>
  dplyr::collect() |>
  wide_to_long_probe_map()

raw_datasets$chronos_somalogic_post <- glue::glue("{get_config()$input_dir}/name=chronos_somalogic_post_onset-pd_master_summary/version=29808/dataset.parquet") |>
  arrow::open_dataset() |>
  dplyr::collect() |>
  wide_to_long_probe_map()

# Olink 5k linear mixed model rerun results
raw_datasets$chronos_olink_all <- glue::glue("{get_config()$input_dir}/name=chronos_olink_all-pd_master_summary/version=29809/dataset.parquet") |>
  arrow::open_dataset() |>
  dplyr::collect() |>
  wide_to_long_probe_map()

raw_datasets$chronos_olink_pre <- glue::glue("{get_config()$input_dir}/name=chronos_olink_all_pre_onset-pd_master_summary/version=29810/dataset.parquet") |>
  arrow::open_dataset() |>
  dplyr::collect() |>
  wide_to_long_probe_map()

raw_datasets$chronos_olink_post <- glue::glue("{get_config()$input_dir}/name=chronos_olink_all_post_onset-pd_master_summary/version=29811/dataset.parquet") |>
  arrow::open_dataset() |>
  dplyr::collect() |>
  wide_to_long_probe_map()

# Alamar linear mixed model rerun results
# Note: need to reformat to UNIPROT-separated line format
# Note: P69905 will be HBA1 (remove HBA2 association)
alamar_all_raw <- glue::glue("{get_config()$input_dir}/name=chronos_alamar_all-pd_master_summary/version=29812/dataset.parquet") |>
  arrow::open_dataset() |>
  dplyr::collect()
alamar_all_raw$gene_symbol[which(alamar_all_raw$gene_symbol == "HBA1:HBA2")] <- "HBA1"
alamar_all_cleaned <- wide_to_long_probe_map(alamar_all_raw, sep = ":")
raw_datasets$chronos_alamar_all <- alamar_all_cleaned
rm(alamar_all_raw, alamar_all_cleaned)

alamar_pre_raw <- glue::glue("{get_config()$input_dir}/name=chronos_alamar_pre_onset-pd_master_summary/version=29813/dataset.parquet") |>
  arrow::open_dataset() |>
  dplyr::collect()
alamar_pre_raw$gene_symbol[which(alamar_pre_raw$gene_symbol == "HBA1:HBA2")] <- "HBA1"
alamar_pre_cleaned <- wide_to_long_probe_map(alamar_pre_raw, sep = ":")
raw_datasets$chronos_alamar_pre <- alamar_pre_cleaned
rm(alamar_pre_raw, alamar_pre_cleaned)

alamar_post_raw <- glue::glue("{get_config()$input_dir}/name=chronos_alamar_post_onset-pd_master_summary/version=29814/dataset.parquet") |>
  arrow::open_dataset() |>
  dplyr::collect()
alamar_post_raw$gene_symbol[which(alamar_post_raw$gene_symbol == "HBA1:HBA2")] <- "HBA1"
alamar_post_cleaned <- wide_to_long_probe_map(alamar_post_raw, sep = ":")
raw_datasets$chronos_alamar_post <- alamar_post_cleaned
rm(alamar_post_raw, alamar_post_cleaned)

# Joint model results
# Note: need to reformat to UNIPROT-separated line format
# Note: rename platform_id to probe_id
# Note: rename p_value to p_val
# Note: we don't have CohensD for joint model
somalogic_jm_raw <- glue::glue("{get_config()$input_dir}/name=chronos_somalogic_joint_model-pd_master_summary/version=29819/dataset.parquet") |>
  arrow::open_dataset() |>
  dplyr::collect()
somalogic_jm_cleaned <- wide_to_long_probe_map(somalogic_jm_raw)
colnames(somalogic_jm_cleaned)[which(colnames(somalogic_jm_cleaned) == "platform_id")] <- "probe_id"
colnames(somalogic_jm_cleaned)[which(colnames(somalogic_jm_cleaned) == "p_value")] <- "p_val"
raw_datasets$chronos_somalogic_jm <- somalogic_jm_cleaned
rm(somalogic_jm_raw, somalogic_jm_cleaned)

olink_jm_raw <- glue::glue("{get_config()$input_dir}/name=chronos_olink_joint_model-pd_master_summary/version=29820/dataset.parquet") |>
  arrow::open_dataset() |>
  dplyr::collect()
olink_jm_cleaned <- wide_to_long_probe_map(olink_jm_raw)
colnames(olink_jm_cleaned)[which(colnames(olink_jm_cleaned) == "platform_id")] <- "probe_id"
colnames(olink_jm_cleaned)[which(colnames(olink_jm_cleaned) == "p_value")] <- "p_val"
raw_datasets$chronos_olink_jm <- olink_jm_cleaned
rm(olink_jm_raw, olink_jm_cleaned)

alamar_jm_raw <- glue::glue("{get_config()$input_dir}/name=chronos_alamar_joint_model-pd_master_summary/version=29821/dataset.parquet") |>
  arrow::open_dataset() |>
  dplyr::collect()
alamar_jm_cleaned <- wide_to_long_probe_map(alamar_jm_raw)
colnames(alamar_jm_cleaned)[which(colnames(alamar_jm_cleaned) == "platform_id")] <- "probe_id"
colnames(alamar_jm_cleaned)[which(colnames(alamar_jm_cleaned) == "p_value")] <- "p_val"
raw_datasets$chronos_alamar_jm <- alamar_jm_cleaned
rm(alamar_jm_raw, alamar_jm_cleaned)

# Note: we still have the bijectivity problem - need to resolve probe --> SYMBOL
raw_datasets$chronos_somalogic_all$gene_symbol[
  which(nchar(raw_datasets$chronos_somalogic_all$gene_symbol) ==
          max(nchar(raw_datasets$chronos_somalogic_all$gene_symbol), na.rm = TRUE))
]

# Check: redundant probes pointing to the same UNIPROT
max(table(raw_datasets$chronos_somalogic_all$probe_id)) == 5

# ----------------------------------------------------------------------------
# Clean up probe annotations (one probe -> one UNIPROT)
# ----------------------------------------------------------------------------

annotation_maps <- load_annotation_mappings_sparse()

cleaned_datasets <- list()
for (q in 1:length(raw_datasets)) {
  cleaned_datasets[[length(cleaned_datasets) + 1]] <- list()
  
  current_dataset <- raw_datasets[[q]]
  cleaned_probes <- prepare_probe_annotations_sparse(
    current_dataset,
    probe.id = "probe_id",
    annot.col = "uniprot_id",
    maps = annotation_maps
  )$data
  
  if (max(table(cleaned_probes$probe_id)) != 1) {
    stop("redundant probes")
  }
  
  cleaned_datasets[[length(cleaned_datasets)]]$data <- cleaned_probes
  cleaned_datasets[[length(cleaned_datasets)]]$universe <- "UNIPROT"
  cleaned_datasets[[length(cleaned_datasets)]]$background <- sort(unique(cleaned_probes$uniprot_id))
  
  rm(current_dataset, cleaned_probes)
}
names(cleaned_datasets) <- names(raw_datasets)

# ----------------------------------------------------------------------------
# Validation checks: probe-to-UNIPROT projection consistency
# ----------------------------------------------------------------------------

# Check Somalogic pre vs post
soma_pre_post_check <- dplyr::left_join(
  data.frame(
    probe_id = cleaned_datasets$chronos_somalogic_pre$data$probe_id,
    uniprot1 = cleaned_datasets$chronos_somalogic_pre$data$uniprot_id
  ),
  data.frame(
    probe_id = cleaned_datasets$chronos_somalogic_post$data$probe_id,
    uniprot2 = cleaned_datasets$chronos_somalogic_post$data$uniprot_id
  ),
  by = "probe_id"
)
if (!identical(soma_pre_post_check$uniprot1, soma_pre_post_check$uniprot2)) {
  stop("projection error")
}
rm(soma_pre_post_check)

# Check Somalogic pre vs all
soma_pre_all_check <- dplyr::left_join(
  data.frame(
    probe_id = cleaned_datasets$chronos_somalogic_pre$data$probe_id,
    uniprot1 = cleaned_datasets$chronos_somalogic_pre$data$uniprot_id
  ),
  data.frame(
    probe_id = cleaned_datasets$chronos_somalogic_all$data$probe_id,
    uniprot2 = cleaned_datasets$chronos_somalogic_all$data$uniprot_id
  ),
  by = "probe_id"
)
if (!identical(soma_pre_all_check$uniprot1, soma_pre_all_check$uniprot2)) {
  stop("projection error")
}
rm(soma_pre_all_check)

# Check Olink pre vs post
olink_pre_post_check <- dplyr::left_join(
  data.frame(
    probe_id = cleaned_datasets$chronos_olink_pre$data$probe_id,
    uniprot1 = cleaned_datasets$chronos_olink_pre$data$uniprot_id
  ),
  data.frame(
    probe_id = cleaned_datasets$chronos_olink_post$data$probe_id,
    uniprot2 = cleaned_datasets$chronos_olink_post$data$uniprot_id
  ),
  by = "probe_id"
)
if (!identical(olink_pre_post_check$uniprot1, olink_pre_post_check$uniprot2)) {
  stop("projection error")
}
rm(olink_pre_post_check)

# Check Olink pre vs all
olink_pre_all_check <- dplyr::left_join(
  data.frame(
    probe_id = cleaned_datasets$chronos_olink_pre$data$probe_id,
    uniprot1 = cleaned_datasets$chronos_olink_pre$data$uniprot_id
  ),
  data.frame(
    probe_id = cleaned_datasets$chronos_olink_all$data$probe_id,
    uniprot2 = cleaned_datasets$chronos_olink_all$data$uniprot_id
  ),
  by = "probe_id"
)
if (!identical(olink_pre_all_check$uniprot1, olink_pre_all_check$uniprot2)) {
  stop("projection error")
}
rm(olink_pre_all_check)

# Check Alamar pre vs post
alamar_pre_post_check <- dplyr::left_join(
  data.frame(
    probe_id = cleaned_datasets$chronos_alamar_pre$data$probe_id,
    uniprot1 = cleaned_datasets$chronos_alamar_pre$data$uniprot_id
  ),
  data.frame(
    probe_id = cleaned_datasets$chronos_alamar_post$data$probe_id,
    uniprot2 = cleaned_datasets$chronos_alamar_post$data$uniprot_id
  ),
  by = "probe_id"
)
if (!identical(alamar_pre_post_check$uniprot1, alamar_pre_post_check$uniprot2)) {
  stop("projection error")
}
rm(alamar_pre_post_check)

# Check Alamar pre vs all
alamar_pre_all_check <- dplyr::left_join(
  data.frame(
    probe_id = cleaned_datasets$chronos_alamar_pre$data$probe_id,
    uniprot1 = cleaned_datasets$chronos_alamar_pre$data$uniprot_id
  ),
  data.frame(
    probe_id = cleaned_datasets$chronos_alamar_all$data$probe_id,
    uniprot2 = cleaned_datasets$chronos_alamar_all$data$uniprot_id
  ),
  by = "probe_id"
)
if (!identical(alamar_pre_all_check$uniprot1, alamar_pre_all_check$uniprot2)) {
  stop("projection error")
}
rm(alamar_pre_all_check)

# Verify: only 1 measurement per probe
max(table(cleaned_datasets$chronos_somalogic_all$data$probe_id))

# ----------------------------------------------------------------------------
# Get ranked entity lists for pathway analysis
# Entity = biological object (here: UNIPROT ID)
# Ranking by p-value; break ties by effect size
# ----------------------------------------------------------------------------

pathway_ranks <- list()

# Flavors for ORA search: upregulated only, downregulated only, general perturbation
for (flavor in c("up", "down", "mix")) {
  # Cohorts: preDx, postDx, all samples, special case of preDx + joint model
  for (cohort in c("pre", "post", "all", "prejm")) {
    if (cohort == "pre") {
      # 11260 UniProts is our full capacity
      background <- sort(unique(c(
        cleaned_datasets$chronos_somalogic_pre$background,
        cleaned_datasets$chronos_olink_pre$background,
        cleaned_datasets$chronos_alamar_pre$background
      )))
      
      soma_pre <- cleaned_datasets$chronos_somalogic_pre$data
      olink_pre <- cleaned_datasets$chronos_olink_pre$data
      alamar_pre <- cleaned_datasets$chronos_alamar_pre$data
      
      combined_data <- data.frame(
        probe_id = c(soma_pre$probe_id, olink_pre$probe_id, alamar_pre$probe_id),
        entity = c(soma_pre$uniprot_id, olink_pre$uniprot_id, alamar_pre$uniprot_id),
        p_val = c(soma_pre$p_val, olink_pre$p_val, alamar_pre$p_val),
        estimate = c(soma_pre$estimate, olink_pre$estimate, alamar_pre$estimate)
      )
      
      rm(soma_pre, olink_pre, alamar_pre)
    }
    
    if (cohort == "post") {
      background <- sort(unique(c(
        cleaned_datasets$chronos_somalogic_post$background,
        cleaned_datasets$chronos_olink_post$background,
        cleaned_datasets$chronos_alamar_post$background
      )))
      
      soma_post <- cleaned_datasets$chronos_somalogic_post$data
      olink_post <- cleaned_datasets$chronos_olink_post$data
      alamar_post <- cleaned_datasets$chronos_alamar_post$data
      
      combined_data <- data.frame(
        probe_id = c(soma_post$probe_id, olink_post$probe_id, alamar_post$probe_id),
        entity = c(soma_post$uniprot_id, olink_post$uniprot_id, alamar_post$uniprot_id),
        p_val = c(soma_post$p_val, olink_post$p_val, alamar_post$p_val),
        estimate = c(soma_post$estimate, olink_post$estimate, alamar_post$estimate)
      )
      
      rm(soma_post, olink_post, alamar_post)
    }
    
    if (cohort == "all") {
      background <- sort(unique(c(
        cleaned_datasets$chronos_somalogic_all$background,
        cleaned_datasets$chronos_olink_all$background,
        cleaned_datasets$chronos_alamar_all$background
      )))
      
      soma_all <- cleaned_datasets$chronos_somalogic_all$data
      olink_all <- cleaned_datasets$chronos_olink_all$data
      alamar_all <- cleaned_datasets$chronos_alamar_all$data
      
      combined_data <- data.frame(
        probe_id = c(soma_all$probe_id, olink_all$probe_id, alamar_all$probe_id),
        entity = c(soma_all$uniprot_id, olink_all$uniprot_id, alamar_all$uniprot_id),
        p_val = c(soma_all$p_val, olink_all$p_val, alamar_all$p_val),
        estimate = c(soma_all$estimate, olink_all$estimate, alamar_all$estimate)
      )
      
      rm(soma_all, olink_all, alamar_all)
    }
    
    if (cohort == "prejm") {
      background <- sort(unique(c(
        cleaned_datasets$chronos_somalogic_pre$background,
        cleaned_datasets$chronos_olink_pre$background,
        cleaned_datasets$chronos_alamar_pre$background,
        cleaned_datasets$chronos_somalogic_jm$background,
        cleaned_datasets$chronos_olink_jm$background,
        cleaned_datasets$chronos_alamar_jm$background
      )))
      
      soma_pre <- cleaned_datasets$chronos_somalogic_pre$data
      olink_pre <- cleaned_datasets$chronos_olink_pre$data
      alamar_pre <- cleaned_datasets$chronos_alamar_pre$data
      soma_jm <- cleaned_datasets$chronos_somalogic_jm$data
      olink_jm <- cleaned_datasets$chronos_olink_jm$data
      alamar_jm <- cleaned_datasets$chronos_alamar_jm$data
      
      combined_data <- data.frame(
        probe_id = c(
          soma_pre$probe_id, olink_pre$probe_id, alamar_pre$probe_id,
          soma_jm$probe_id, olink_jm$probe_id, alamar_jm$probe_id
        ),
        entity = c(
          soma_pre$uniprot_id, olink_pre$uniprot_id, alamar_pre$uniprot_id,
          soma_jm$uniprot_id, olink_jm$uniprot_id, alamar_jm$uniprot_id
        ),
        p_val = c(
          soma_pre$p_val, olink_pre$p_val, alamar_pre$p_val,
          soma_jm$p_val, olink_jm$p_val, alamar_jm$p_val
        ),
        estimate = c(
          soma_pre$estimate, olink_pre$estimate, alamar_pre$estimate,
          soma_jm$estimate, olink_jm$estimate, alamar_jm$estimate
        )
      )
      
      rm(soma_pre, olink_pre, alamar_pre, soma_jm, olink_jm, alamar_jm)
    }
    
    # JM does not have CohensD - need to break ties based on estimate
    combined_data$rank <- multi_rank(combined_data$p_val, -abs(combined_data$estimate))
    
    # Filter by directionality
    if (flavor == "up") {
      combined_data$rank[which(combined_data$estimate < 0)] <- NA
    }
    if (flavor == "down") {
      combined_data$rank[which(combined_data$estimate > 0)] <- NA
    }
    
    # Get best rank per entity
    entity_ranks <- data.frame(entity = unique(combined_data$entity), v = NA)
    for (k in 1:nrow(entity_ranks)) {
      entity_ranks$v[k] <- collapse::fmin(
        combined_data$rank[which(combined_data$entity == entity_ranks$entity[k])],
        na.rm = TRUE
      )
    }
    entity_ranks <- entity_ranks[order(entity_ranks$v), ]
    entity_ranks$rank <- rank(entity_ranks$v, ties.method = "min", na.last = "keep")
    
    pathway_ranks[[length(pathway_ranks) + 1]] <- list()
    names(pathway_ranks)[length(pathway_ranks)] <- paste0(
      "lme_SomaOlinkAlamar_", cohort, "_", flavor
    )
    pathway_ranks[[length(pathway_ranks)]]$universe <- "UNIPROT"
    pathway_ranks[[length(pathway_ranks)]]$ranks <- dplyr::select(entity_ranks, c("entity", "rank"))
    pathway_ranks[[length(pathway_ranks)]]$background <- background
    
    rm(combined_data, entity_ranks, background)
  }
}

print("data prepared for analysis")
